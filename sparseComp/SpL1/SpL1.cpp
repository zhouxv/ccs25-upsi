#include "./SpL1.h"
#include "../Common/ZN.h"
#include "../SpBSOT/SpBSOT.h"
#include "../BlockSpBSOT/BlockSpBSOT.h"
#include "../CustomOPRF/CustomizedOPRF.h"
#include <array>
#include <cmath>

#define MAX_SSP 128
#define IN_COMP_BIT_LEN 8

template<typename T, size_t N>
using array = std::array<T, N>;

template<size_t tr, size_t ts, size_t k, size_t n>
using BlockSpBSOTSender = sparse_comp::block_sp_bsot::Sender<tr,ts,k,n>; 

template<size_t ts, size_t tr, size_t k, size_t n>
using BlockSpBSOTReceiver = sparse_comp::block_sp_bsot::Receiver<ts,tr,k,n>; 

template<size_t tr, size_t t, size_t k, size_t n, uint64_t M>
using SpBSOTSender = sparse_comp::sp_bsot::Sender<tr,t,k,n,M>;

template<size_t ts, size_t t, size_t k, size_t n, uint64_t M>
using SpBSOTReceiver = sparse_comp::sp_bsot::Receiver<ts,t,k,n,M>;

using OprfSender = sparse_comp::custom_oprf::Sender;
using OprfReceiver = sparse_comp::custom_oprf::Receiver;

template<size_t t, size_t d, uint64_t N>
static void gen_zero_shares(array<array<ZN<N>,d>,t>& zero_shares) {

    for (size_t i=0;i < t;i++) {
        for (size_t j=0;j < d;j++) {
            zero_shares[i][j] = ZN<N>(0);
        }
    }

}
/*
template<size_t t>
void hash_z_shares(AES& aes, array<array<block,1>,t>& z_vec_shares, array<block,t>& hashed_z_shares) {

    for (size_t i=0;i < t;i++) {
        
        hashed_z_shares[i] = aes.hashBlock(z_vec_shares[i][0]);

    }

}
*/

template<size_t t>
void hash_z_shares(AES& aes, array<array<block,1>,t>& z_vec_shares, vector<block>& hashed_z_shares) {

    for (size_t i=0;i < t;i++) {
        
        hashed_z_shares[i] = aes.hashBlock(z_vec_shares[i][0]);

    }

}

template<size_t tr, size_t ts, size_t d, uint16_t twotol, uint8_t delta, uint64_t M>
static Proto sender_compute_h_shares(OprfSender* oprfSender,
                                    Socket& sock, 
                                    PRNG& prng,
                                    array<point,ts>& ordIndexSet,
                                    array<array<ZN<twotol>,d>,ts>& in_vals,
                                    array<array<ZN<M>,d>,ts>& h_shares) {
    
    static_assert(M == d*(delta + 1) + 1,"the following identity must be fulfilled: M = d*(delta + 1) + 1");

    MC_BEGIN(Proto, oprfSender, &sock, &prng, &ordIndexSet,&in_vals,&h_shares,
             zero_shares = (array<array<ZN<twotol>,d>,ts>*) nullptr,
             msg_vecs = (array<array<array<ZN<M>,twotol>,d>,ts>*) nullptr,
             bsotSender = (SpBSOTSender<tr,ts,d,twotol,M>*) nullptr,
             diff = int32_t(0),
             abs_diff = int32_t(0));
        zero_shares = new array<array<ZN<twotol>,d>,ts>();
        msg_vecs = new array<array<array<ZN<M>,twotol>,d>,ts>();
        bsotSender = new SpBSOTSender<tr,ts,d,twotol,M>(prng, oprfSender);

        for (size_t i=0;i < ts;i++) {
            array<array<ZN<M>,twotol>,d>& msg_vec_batch = msg_vecs->at(i);

            for (size_t j=0;j < d;j++) {
                array<ZN<M>,twotol>& msg_vec = msg_vec_batch[j];
                
                for (int32_t h=0;h < twotol;h++) {

                    diff = h-((int32_t) in_vals[i][j].to_uint64_t());
                    abs_diff = (int32_t) abs(diff);
                    msg_vec[h] = (abs_diff <= delta) ? ZN<M>(abs_diff) : ZN<M>(delta+1);
                    
                }

            }

        }

        gen_zero_shares<ts,d,twotol>(*zero_shares);

        MC_AWAIT(bsotSender->send(sock, ordIndexSet, *msg_vecs, *zero_shares, h_shares));

        delete zero_shares;
        delete msg_vecs;
        delete bsotSender;

    MC_END();

}


template<size_t ts, size_t tr, size_t d, uint16_t twotol, uint8_t delta, uint64_t M>
static Proto recvr_compute_h_shares(OprfReceiver* oprfReceiver,
                                   Socket& sock, 
                                   PRNG& prng,
                                   array<point,tr>& ordIndexSet,
                                   array<array<ZN<twotol>,d>,tr>& in_vals,
                                   array<array<ZN<M>,d>,tr>& h_shares) {
    
    static_assert(M == d*(delta + 1) + 1,"the following identity must be fulfilled: M = d*(delta + 1) + 1");
    
    MC_BEGIN(Proto, oprfReceiver, &sock, &prng, &ordIndexSet, &in_vals, &h_shares,
             bsotReceiver = (SpBSOTReceiver<ts, tr,d,twotol,M>*) nullptr);
        
        bsotReceiver = new SpBSOTReceiver<ts, tr,d,twotol,M>(prng, oprfReceiver);
       
        MC_AWAIT(bsotReceiver->receive(sock, ordIndexSet, in_vals, h_shares));

        delete bsotReceiver;

    MC_END();

}

template<size_t tr, size_t ts, size_t d, uint8_t delta, uint64_t M>
static Proto sender_comp_z_shares(OprfSender* oprfSender,
                                  Socket& sock, 
                                  PRNG& prng,
                                  array<point,ts>& ordIndexSet,
                                  array<array<ZN<M>,1>,ts>& g_vec_shares,
                                  array<array<block,1>,ts>& z_vec_shares) {
    static_assert(M == d*(delta + 1) + 1,"the following identity must be fulfilled: M = d*(delta + 1) + 1");

    MC_BEGIN(Proto, oprfSender, &sock, &prng, &ordIndexSet, &g_vec_shares, &z_vec_shares,
             msg_vecs = (array<array<array<block,M>,1>,ts>*) nullptr,
             bsotSender = (BlockSpBSOTSender<tr,ts, 1, M>*) nullptr,
             r = block(0,0));
        msg_vecs = new array<array<array<block,M>,1>,ts>();
        bsotSender = new BlockSpBSOTSender<tr,ts, 1, M>(prng, oprfSender);

        for (size_t i=0;i < ts;i++) {
            r = prng.get<block>();
            for (size_t j=0;j < M;j++) {
                if(j <= delta) {
                    msg_vecs->at(i)[0][j] = block(0,0);
                } else {
                    msg_vecs->at(i)[0][j] = r;
                }
            }
        }

        MC_AWAIT(bsotSender->send(sock, ordIndexSet, *msg_vecs, g_vec_shares, z_vec_shares));

        delete msg_vecs;
        delete bsotSender;

    MC_END();
}

template<size_t ts, size_t tr, size_t d, uint8_t delta, uint64_t M>
static Proto receiver_comp_z_shares(OprfReceiver* oprfReceiver,
                                    Socket& sock, 
                                    PRNG& prng,
                                    array<point,tr>& ordIndexSet,
                                    array<array<ZN<M>,1>,tr>& g_vec_shares,
                                    array<array<block,1>,tr>& z_vec_shares) {
    static_assert(M == d*(delta + 1) + 1,"the following identity must be fulfilled: M = d*(delta + 1) + 1");

    MC_BEGIN(Proto, oprfReceiver, &sock, &prng, &ordIndexSet, &g_vec_shares, &z_vec_shares,
             bsotReceiver = (BlockSpBSOTReceiver<ts,tr, 1, M>*) nullptr);

        bsotReceiver = new BlockSpBSOTReceiver<ts,tr, 1, M>(prng, oprfReceiver);

        MC_AWAIT(bsotReceiver->receive(sock, ordIndexSet, g_vec_shares, z_vec_shares));

        delete bsotReceiver;

    MC_END();
}

template<size_t t, size_t d, uint8_t delta, uint64_t M>
static void comp_g_shares(array<array<ZN<M>,d>,t>& h_vec_shares, array<array<ZN<M>,1>,t>& g_shares) {

    static_assert(M == d*(delta + 1) + 1,"the following identity must be fulfilled: M = d*(delta + 1) + 1");

    for (size_t i=0;i < t;i++) {
        g_shares[i][0] = ZN<M>(0);

        for (size_t j=0;j < d;j++) {
            
            g_shares[i][0] = g_shares[i][0] + h_vec_shares[i][j];

        }

    }

}

template<size_t ts, size_t tr>
static void compute_intersec_hashed_z_shares(AES& aes, 
                                             array<array<block,1>,tr>& rec_z_shares, 
                                             vector<block>& hashed_sen_z_shares,
                                             vector<size_t>& inter_pos) {

    vector<block>* hashed_rec_z_shares = new vector<block>(tr);
    
    hash_z_shares(aes, rec_z_shares, *hashed_rec_z_shares);

    std::unordered_set<block> set;

    for (size_t i=0;i < ts;i++) {
        set.insert(hashed_sen_z_shares[i]);
    }

    for(size_t i=0;i < tr;i++) {
        if(set.contains(hashed_rec_z_shares->at(i))) {
            inter_pos.push_back(i);
        }
    }

    delete hashed_rec_z_shares;

}

template<size_t t, size_t d, uint64_t N>
static void in_values_to_zn(array<array<uint32_t,d>,t>& in_values, array<array<ZN<N>,d>,t>& convted_values) {

    for (size_t i=0;i < t;i++) {
        for (size_t j=0;j < d;j++) {
            convted_values[i][j] = ZN<N>(in_values[i][j]);
        }
    }

}

template<size_t tr, size_t ts, size_t d, uint8_t delta, uint8_t ssp>
Proto sparse_comp::sp_l1::Sender<tr,ts,d,delta,ssp>::send(
                                                    Socket& sock, 
                                                    array<point,ts>& ordIndexSet, 
                                                    array<array<uint32_t,d>,ts>& in_values) {
    static_assert(ssp <= MAX_SSP,"ssp must be less or equal to 128");
    
    constexpr uint8_t l = IN_COMP_BIT_LEN;
    static_assert(l <= 8,"l must be less or equal to 8");
    constexpr const uint16_t twotol = (uint16_t) pow(2, l);
    constexpr const uint64_t M = d*(delta + 1) + 1;

    //constexpr const uint64_t two_to_ssp = std::pow(2,ssp);

    constexpr const size_t oprf_instances = 2;

    MC_BEGIN(Proto, this, &sock, &ordIndexSet, &in_values,
             zn_in_values = (array<array<ZN<twotol>,d>,ts>*) nullptr,
             h_vec_shares = (array<array<ZN<M>,d>,ts>*) nullptr,
             g_vec_shares = (array<array<ZN<M>,1>,ts>*) nullptr,
             z_vec_shares = (array<array<block,1>,ts>*) nullptr,
             hashed_z_shares = vector<block>(),
             oprfSenders = std::vector<OprfSender*>(oprf_instances),
             prt = Proto());
        

        MC_AWAIT(OprfSender::setup(sock, *(this->prng), oprf_instances, oprfSenders)); // Setup OPRFs

        zn_in_values = new array<array<ZN<twotol>,d>,ts>();
        in_values_to_zn<ts,d,twotol>(in_values, *zn_in_values);

        h_vec_shares = new array<array<ZN<M>,d>,ts>();
        prt = sender_compute_h_shares<tr,ts,d,twotol,delta,M>(oprfSenders[0], sock, *(this->prng), ordIndexSet, *zn_in_values, *h_vec_shares);
        MC_AWAIT(prt);

        delete zn_in_values;
        
        //std::cout << "M: " << M << std::endl;
        //std::cout << "sender h vec shares" << std::endl;
        //std::cout << h_vec_shares[1][0].to_uint64_t() << "," << h_vec_shares[1][1].to_uint64_t() << std::endl;

        g_vec_shares = new array<array<ZN<M>,1>,ts>();
        comp_g_shares<ts,d,delta,M>(*h_vec_shares, *g_vec_shares);
        delete h_vec_shares;

        //std::cout << "sender g vec shares" << std::endl;
        //std::cout << g_vec_shares[1][0].to_uint64_t() << std::endl;

        z_vec_shares = new array<array<block,1>,ts>();
        prt = sender_comp_z_shares<tr,ts,d,delta,M>(oprfSenders[1], sock, *(this->prng), ordIndexSet, *g_vec_shares, *z_vec_shares);
        MC_AWAIT(prt);
        delete g_vec_shares;

        hashed_z_shares.resize(ts);
        hash_z_shares<ts>(*(this->aes), *z_vec_shares, hashed_z_shares);

        delete z_vec_shares;

        prt = sparse_comp::send<block,sparse_comp::COPROTO_MAX_SEND_SIZE_BYTES>(sock, hashed_z_shares);

        MC_AWAIT(prt);

    MC_END();
}

template<size_t ts, size_t tr, size_t d, uint8_t delta, uint8_t ssp>
Proto sparse_comp::sp_l1::Receiver<ts,tr,d,delta,ssp>::receive(
                                                    Socket& sock, 
                                                    array<point,tr>& ordIndexSet, 
                                                    array<array<uint32_t,d>,tr>& in_values, 
                                                    vector<size_t>& intersec_pos) {
    static_assert(ssp <= MAX_SSP,"ssp must be less or equal to 128");

    constexpr const uint8_t l = IN_COMP_BIT_LEN;
    static_assert(l <= 8,"l must be less or equal to 8");
    constexpr const uint16_t twotol = (uint16_t) pow(2, l);
    constexpr const uint64_t M = d*(delta + 1) + 1;

    //constexpr const uint64_t two_to_ssp = std::pow(2,ssp);

    constexpr const size_t oprf_instances = 2;

    MC_BEGIN(Proto, this, &sock, &ordIndexSet, &in_values, &intersec_pos,
             zn_in_values = (array<array<ZN<twotol>,d>,tr>*) nullptr,
             h_vec_shares = (array<array<ZN<M>,d>,tr>*) nullptr,
             g_vec_shares = (array<array<ZN<M>,1>,tr>*) nullptr,
             z_vec_shares =  (array<array<block,1>,tr>*) nullptr,
             sender_hashed_z_sender_shares = vector<block>(),
             oprfReceivers = std::vector<OprfReceiver*>(oprf_instances),
             prt = Proto());

        MC_AWAIT(OprfReceiver::setup(sock, *(this->prng), oprf_instances, oprfReceivers)); // Setup OPRFs

        zn_in_values = new array<array<ZN<twotol>,d>,tr>();
        h_vec_shares = new array<array<ZN<M>,d>,tr>();

        in_values_to_zn<tr,d,twotol>(in_values, *zn_in_values);

        prt = recvr_compute_h_shares<ts,tr,d,twotol,delta,M>(oprfReceivers[0], sock, *(this->prng), ordIndexSet, *zn_in_values, *h_vec_shares);
        MC_AWAIT(prt);

        delete zn_in_values;

        //std::cout << "receiver h vec shares" << std::endl;
        //std::cout << h_vec_shares->at(0)[0].to_uint64_t() << "," << h_vec_shares->at(0)[1].to_uint64_t() << std::endl;
        
        g_vec_shares = new array<array<ZN<M>,1>,tr>();

        comp_g_shares<tr,d,delta,M>(*h_vec_shares,*g_vec_shares);
        delete h_vec_shares;

        //std::cout << "receiver g vec shares" << std::endl;
        //std::cout << g_vec_shares->at(0)[0].to_uint64_t() << std::endl;

        z_vec_shares = new array<array<block,1>,tr>();
        prt = receiver_comp_z_shares<ts,tr,d,delta,M>(oprfReceivers[1], sock, *(this->prng), ordIndexSet,*g_vec_shares,*z_vec_shares);
        MC_AWAIT(prt);
        delete g_vec_shares;

        sender_hashed_z_sender_shares.resize(ts);

        prt = sparse_comp::receive<block,sparse_comp::COPROTO_MAX_SEND_SIZE_BYTES>(sock, ts, sender_hashed_z_sender_shares);
        MC_AWAIT(prt);

        compute_intersec_hashed_z_shares<ts,tr>(*(this->aes), *z_vec_shares, sender_hashed_z_sender_shares, intersec_pos);

        delete z_vec_shares;
    
    MC_END();
}