#include "./MultiOPRF.h"
#include "libOTe/Base/BaseOT.h"
#include "libOTe/TwoChooseOne/Kos/KosOtExtReceiver.h"
#include "libOTe/TwoChooseOne/Kos/KosOtExtSender.h"
#include "cryptoTools/Crypto/PRNG.h"
#include "cryptoTools/Common/BitVector.h"
#include "volePSI/Paxos.h"
#include "../Common/BaxosUtils.h"
#include "../Common/SockUtils.h"
#include <vector>

#define MULTI_OPRF_PAXOS_SSP 40
#define MULTI_OPRF_PAXOS_NTHREADS 1

static constexpr int64_t max_send_size_bytes(2147483648L); // 2 GBs

using KosOtExtSender = osuCrypto::KosOtExtSender;
using KosOtExtReceiver = osuCrypto::KosOtExtReceiver;
using u8 = osuCrypto::u8;
using Baxos = volePSI::Baxos;
using PaxosParam = volePSI::PaxosParam;

static const size_t comp_sec = sparse_comp::multi_oprf::comp_sec_param;
static const size_t ell = sparse_comp::multi_oprf::ell;

sparse_comp::multi_oprf::Sender::~Sender() {
    delete this->randSetupOtMsgs;
    delete this->okvs;
}

Proto sparse_comp::multi_oprf::Sender::setup(coproto::Socket& sock, PRNG& prng, size_t num_instances, std::vector<Sender*>& senders) {
    MC_BEGIN(Proto, &sock, &prng, num_instances, &senders,
             otRecv = KosOtExtReceiver(),
             allRandSetupOtMsgs = (std::vector<block>*) nullptr,
             allRandSetupOtChoices = (BitVector*) nullptr,
             choice_blocks = (block*) nullptr);
            
        senders.resize(num_instances);    
        
        otRecv.mIsMalicious = false;

        allRandSetupOtMsgs = new std::vector<block>(ell*num_instances);
        allRandSetupOtChoices = new BitVector(ell*num_instances);

        allRandSetupOtChoices->randomize(prng); // Sample random setup ot choices

        MC_AWAIT(otRecv.genBaseOts(prng, sock));

        MC_AWAIT(otRecv.receive(*allRandSetupOtChoices, *allRandSetupOtMsgs, prng, sock));

        choice_blocks = allRandSetupOtChoices->blocks();

        for (size_t i=0;i < num_instances;i++) {
            senders[i] = new Sender();
            senders[i]->s = choice_blocks[i];
            senders[i]->randSetupOtMsgs = new std::vector<block>(ell);
            
            for (size_t j=0;j < ell;j++) {
                senders[i]->randSetupOtMsgs->at(j) = allRandSetupOtMsgs->at(ell*i + j);
            }
        }

        delete allRandSetupOtChoices;
        delete allRandSetupOtMsgs;

    MC_END();
}

sparse_comp::multi_oprf::Receiver::~Receiver() {
    delete this->randSetupOtMsgs;
}

Proto sparse_comp::multi_oprf::Receiver::setup(coproto::Socket& sock, PRNG& prng, size_t num_instances, std::vector<Receiver*>& receivers) {
    MC_BEGIN(Proto, &sock, &prng, num_instances, &receivers,
             allRandSetupOtMsgs = (std::vector<std::array<block, 2>>*) nullptr,
             otSender = KosOtExtSender());

        otSender.mIsMalicious = false;

        allRandSetupOtMsgs = new std::vector<std::array<block, 2>>(ell*num_instances);
        prng.get((u8*)allRandSetupOtMsgs->data()->data(), sizeof(block) * 2 * allRandSetupOtMsgs->size());

        MC_AWAIT(otSender.genBaseOts(prng, sock));
        MC_AWAIT(otSender.send(*allRandSetupOtMsgs, prng, sock));

        for (size_t i=0;i < num_instances;i++) {
            receivers[i] = new Receiver();
            receivers[i]->randSetupOtMsgs = new std::vector<std::array<block, 2>>(ell);

            for (size_t j=0;j < ell;j++) {
                receivers[i]->randSetupOtMsgs->at(j)[0] = allRandSetupOtMsgs->at(ell*i + j)[0];
                receivers[i]->randSetupOtMsgs->at(j)[1] = allRandSetupOtMsgs->at(ell*i + j)[1];
            }
        }

        delete allRandSetupOtMsgs;

    MC_END();
}

static void compute_R_blocks(size_t ell,
                              std::vector<std::array<block, 2>>& randSetupOtMsgs,
                              std::vector<block>& ys,
                              std::vector<block>& rs) {
    std::vector<block> ciphers(ys.size());
    std::vector<block> bit1_blocks(128);
    std::vector<AES> t_prfs(ell);
    std::vector<AES> u_prfs(ell);

    for (size_t i=0;i < ell;i++) { //Instantiate PRFs using ot messages as keys.
        t_prfs[i].setKey(randSetupOtMsgs[i][0]);
        u_prfs[i].setKey(randSetupOtMsgs[i][1]);
    }

    for (size_t i=0;i < 128;i++) { // Instantiate blocks with one hot bit at position i.
        bit1_blocks[i] = i < 64 ? block(0,1) : block(1,0);

        bit1_blocks[i] = bit1_blocks[i] << (i % 64);
    }

    for (size_t i=0;i < ell;i++) {
        t_prfs[i].ecbEncBlocks(ys.data(), ys.size(), ciphers.data());
        for (size_t j=0;j < ys.size();j++) {
            rs[j] = rs[j] ^ (ciphers[j] & bit1_blocks[i]);
        }

        u_prfs[i].ecbEncBlocks(ys.data(), ys.size(), ciphers.data());
        for (size_t j=0;j < ys.size();j++) {
            rs[j] = rs[j] ^ (ciphers[j] & bit1_blocks[i]);
        }
    }

}

static void compute_T_blocks(size_t ell,
                              std::vector<std::array<block, 2>>& randSetupOtMsgs,
                              std::vector<block>& ys,
                              std::vector<block>& ts) {
    std::vector<block> ciphers(ys.size());
    std::vector<block> bit1_blocks(128);
    std::vector<AES> t_prfs(ell);

    for (size_t i=0;i < ell;i++) { //Instantiate PRFs using ot messages as keys.
        t_prfs[i].setKey(randSetupOtMsgs[i][0]);
    }

    for (size_t i=0;i < 128;i++) { // Instantiate blocks with one hot bit at position i.
        bit1_blocks[i] = i < 64 ? block(0,1) : block(1,0);

        bit1_blocks[i] = bit1_blocks[i] << (i % 64);
    }

    for (size_t i=0;i < ell;i++) {
        t_prfs[i].ecbEncBlocks(ys.data(), ys.size(), ciphers.data());
        for (size_t j=0;j < ys.size();j++) {
            ts[j] = ts[j] ^ (ciphers[j] & bit1_blocks[i]);
        }
    }

}

static void compute_Q_blocks(size_t ell,
                             std::vector<block>& randSetupOtMsgs,
                             std::vector<block>& xs,
                             std::vector<block>& qs) {
    block* ciphers = new block[xs.size()]; 
    AES q_prf;

  /*  for (size_t i=0;i < ell;i++) { //Instantiate PRFs using ot messages as keys.
        q_prfs[i].setKey(randSetupOtMsgs[i]);
    }

    for (size_t i=0;i < 128;i++) { // Instantiate blocks with one hot bit at position i.
        bit1_blocks[i] = i < 64 ? block(0,1) : block(1,0);

        bit1_blocks[i] = bit1_blocks[i] << (i % 64);
    }
    */


    //std::chrono::high_resolution_clock::time_point start = std::chrono::high_resolution_clock::now();
    for (size_t i=0;i < 64;i++) {
        block bit1_block = block(0,1) << i; 
        q_prf.setKey(randSetupOtMsgs[i]);

        q_prf.ecbEncBlocks(xs.data(), xs.size(), ciphers);
        for (size_t j=0;j < xs.size();j++) {
            qs[j] = qs[j] ^ (ciphers[j] & bit1_block);
        }
    }

    for (size_t i=0;i < 64;i++) {
        block bit1_block = block(1,0) << i; 
        q_prf.setKey(randSetupOtMsgs[i+64]);

        q_prf.ecbEncBlocks(xs.data(), xs.size(), ciphers);
        for (size_t j=0;j < xs.size();j++) {
            qs[j] = qs[j] ^ (ciphers[j] & bit1_block);
        }
    }

    delete[] ciphers;

    //std::chrono::high_resolution_clock::time_point end = std::chrono::high_resolution_clock::now();
    //std::chrono::duration<double> elapsed = end - start;
    //std::cout << "Time taken to execute compute_Q_blocks third loop: " << elapsed.count() << " seconds" << std::endl;

}

Proto sparse_comp::multi_oprf::Sender::send(coproto::Socket& sock, size_t query_num) {
    MC_BEGIN(Proto, this, &sock, query_num,
             paxosBlockCount = size_t(0),
             t = coproto::task<void>{});

        //std::cout << "waiting to receive multioprf okvs (s)" << std::endl;

        this->query_num = query_num;

        paxosBlockCount = sparse_comp::baxosBlockCount(query_num, MULTI_OPRF_PAXOS_SSP);

        t = sparse_comp::receive<block, max_send_size_bytes>(sock, paxosBlockCount, *okvs);
        MC_AWAIT(t);

        //MC_AWAIT(sock.recvResize(*(this->okvs)));

        //std::cout << "received multioprf okvs (s); okvs byte size: " << (this->okvs->size() * sizeof(block)) << std::endl;

    MC_END();
}

static void decode_okvs(std::vector<block>& idxs, std::vector<block>& vals, size_t okvs_num_encoded, std::vector<block>& okvs) {

    Baxos paxos;
    paxos.init(okvs_num_encoded, sparse_comp::baxosBinSize(okvs_num_encoded), 3, MULTI_OPRF_PAXOS_SSP, PaxosParam::GF128, oc::ZeroBlock);
    
    paxos.decode<block>(idxs, vals, okvs);

}

void sparse_comp::multi_oprf::Sender::eval(std::vector<block>& idxs, std::vector<block>& vals) {
    std::vector<block> qs(idxs.size());
    std::vector<block> ps(idxs.size());
    std::vector<block> vs(idxs.size());

    //std::chrono::high_resolution_clock::time_point start = std::chrono::high_resolution_clock::now();
    compute_Q_blocks(ell, *(this->randSetupOtMsgs), idxs, qs);
    //std::chrono::high_resolution_clock::time_point end = std::chrono::high_resolution_clock::now();
    //std::chrono::duration<double> elapsed = end - start;
    //std::cout << "Time taken to execute compute_Q_blocks: " << elapsed.count() << " seconds" << std::endl;

    //start = std::chrono::high_resolution_clock::now();
    decode_okvs(idxs, ps, this->query_num, *(this->okvs));
    //end = std::chrono::high_resolution_clock::now();
    //elapsed = end - start;
    //std::cout << "Time taken to execute decode_okvs: " << elapsed.count() << " seconds" << std::endl;

    //start = std::chrono::high_resolution_clock::now();
    for (size_t i=0;i < idxs.size();i++) {
        vs[i] = qs[i] ^ ((this->s) & (ps[i]));
    }
    //end = std::chrono::high_resolution_clock::now();
    //elapsed = end - start;
    //std::cout << "Time taken to execute vs->at(i) = qs->at(i) ^ ((this->s) & (ps->at(i)));: " << elapsed.count() << " seconds" << std::endl;

    //start = std::chrono::high_resolution_clock::now();
    this->aes.hashBlocks(vs.data(), vs.size(), vals.data());
    //end = std::chrono::high_resolution_clock::now();
    //elapsed = end - start;
    //std::cout << "Time taken to execute >aes.hashBlocks: " << elapsed.count() << " seconds" << std::endl;
}

static void encode_okvs(std::vector<block>& idxs, std::vector<block>& vals, std::vector<block>& okvs) {

    Baxos paxos;
    paxos.init(idxs.size(), sparse_comp::baxosBinSize(idxs.size()), 3, MULTI_OPRF_PAXOS_SSP, PaxosParam::GF128, oc::ZeroBlock);
    okvs.resize(paxos.size());

    paxos.solve<block>(idxs, vals, okvs, nullptr, MULTI_OPRF_PAXOS_NTHREADS);

}

Proto sparse_comp::multi_oprf::Receiver::receive(coproto::Socket& sock, std::vector<block>& idxs, std::vector<block>& vals) {
    MC_BEGIN(Proto, this, &sock, &idxs, &vals,
            rs = (std::vector<block>*) nullptr,
            ts = (std::vector<block>*) nullptr,
            okvs = (std::vector<block>*) nullptr,
            ec = macoro::result<void>{},
            i = size_t(0),
            t = coproto::task<void>{},
            start = std::chrono::high_resolution_clock::time_point{},
            end = std::chrono::high_resolution_clock::time_point{},
            elapsed = std::chrono::duration<double>{});

        rs = new std::vector<block>(idxs.size());
        ts = new std::vector<block>(idxs.size());
        okvs = new std::vector<block>();

        //start = std::chrono::high_resolution_clock::now();
        compute_R_blocks(ell, *(this->randSetupOtMsgs), idxs, *rs);
        //end = std::chrono::high_resolution_clock::now();
        //elapsed = end - start;
        //std::cout << "Time taken to execute compute_R_blocks: " << elapsed.count() << " seconds" << std::endl;

        //compute_R_blocks(ell, *(this->randSetupOtMsgs), idxs, *rs);
        //start = std::chrono::high_resolution_clock::now();
        encode_okvs(idxs, *rs, *okvs);
        //end = std::chrono::high_resolution_clock::now();
        //elapsed = end - start;
        //std::cout << "Time taken to execute encode_okvs: " << elapsed.count() << " seconds" << std::endl;
        //encode_okvs(idxs, *rs, *okvs);

        // std::cout << "sending multioprf okvs (r)" << std::endl;

        //MC_AWAIT(sock.send(*okvs));

        // std::cout << "okvs byte size: " << (okvs->size() * sizeof(block)) << std::endl;

        t = sparse_comp::send<block, max_send_size_bytes>(sock, *okvs);
        MC_AWAIT(t);

        //MC_AWAIT_SET(ec, sock.send(std::move(*okvs)) | macoro::wrap());

        // std::cout << "multioprf okvs sent (r)" << std::endl;

        //start = std::chrono::high_resolution_clock::now();
        compute_T_blocks(ell, *(this->randSetupOtMsgs), idxs, *ts);
        //end = std::chrono::high_resolution_clock::now();
        //elapsed = end - start;
        //std::cout << "Time taken to execute compute_T_blocks: " << elapsed.count() << " seconds" << std::endl;
        //compute_T_blocks(ell, *(this->randSetupOtMsgs),idxs,*ts);

        //start = std::chrono::high_resolution_clock::now();
        this->aes.hashBlocks(ts->data(), ts->size(), vals.data());

        //end = std::chrono::high_resolution_clock::now();
        //elapsed = end - start;
        //std::cout << "Time taken to execute >aes.hashBlocks: " << elapsed.count() << " seconds" << std::endl;


        delete rs;
        delete ts;
        delete okvs;

    MC_END();
}

/*
 MC_BEGIN(Proto, this, &sock, &prng, base_ot_num, max_query_num);

        std::vector<block> recvMsg(max_query_num), baseRecv(base_ot_num);
        this->randSetupOtMsgs = new std::vector<std::array<block, 2>>(max_query_num);
        BitVector baseChoice(base_ot_num);
        
        baseChoice.randomize(prng); // Sample random base choices
        prng.get((u8*)sendMsg->data()->data(), sizeof(block) * 2 * sendMsg->size()); // Sample random OT messages.

        IknpOtExtSender sender;

        sender.setBaseOts(baseRecv, baseChoice);
        MC_AWAIT(sender.send(*sendMsg, prng, sock));



    MC_END();
*/