#include "CustomizedOPRF.h"
#include "volePSI/RsOprf.h"
#include "cryptoTools/Common/Defines.h"
#include "cryptoTools/Common/block.h"
#include "../Common/HashUtils.h"
#include "../Common/Common.h"
#include <vector>
#include <iostream>


using PRNG = osuCrypto::PRNG;
using oprf_point = sparse_comp::custom_oprf::oprf_point;
using MultiOprfRecvr = sparse_comp::multi_oprf::Receiver;
using MultiOprfSender = sparse_comp::multi_oprf::Sender;


template<typename T>
using vector = std::vector<T>;

Proto sparse_comp::custom_oprf::Sender::setup(coproto::Socket& sock, PRNG& prng, size_t num_instances, std::vector<Sender*>& senders) {
    MC_BEGIN(Proto, &sock, &prng, num_instances, &senders,
             oprfSenders = (std::vector<MultiOprfSender*>*) nullptr,
             i = (size_t) 0);
    
        oprfSenders = new std::vector<MultiOprfSender*>(num_instances);

        MC_AWAIT(MultiOprfSender::setup(sock, prng, num_instances, *oprfSenders));

        for(i=0;i < num_instances;i++) {
            senders[i] = new Sender(oprfSenders->at(i));
        } 

        delete oprfSenders;

    MC_END();
}

Proto sparse_comp::custom_oprf::Receiver::setup(coproto::Socket& sock, PRNG& prng, size_t num_instances, std::vector<Receiver*>& receivers) {
    MC_BEGIN(Proto, &sock, &prng, num_instances, &receivers,
             oprfReceivers = (std::vector<MultiOprfRecvr*>*) nullptr,
             i = (size_t) 0);
    
        oprfReceivers = new std::vector<MultiOprfRecvr*>(num_instances);
        // std::vector<MultiOprfRecvr*> oprfReceivers(num_instances);

        MC_AWAIT(MultiOprfRecvr::setup(sock, prng, num_instances, *oprfReceivers));

        for(i=0;i < num_instances;i++) {
            receivers[i] = new Receiver(oprfReceivers->at(i));
        } 

        delete oprfReceivers;

    MC_END();
}

block encode_point_as_block(const AES& aes,const sparse_comp::point& pot, size_t sot_idx, size_t msg_vec_idx) {

    return sparse_comp::hash_point(aes, pot, sot_idx, msg_vec_idx);

}

/*block encode_point_as_block(const AES& aes,const sparse_comp::point& pot, size_t sot_idx) {

    return sparse_comp::hash_point(aes, pot, sot_idx);

}
*/
sparse_comp::custom_oprf::Sender::Sender(MultiOprfSender* oprfSender) {
    this->oprfSender = oprfSender;
}

sparse_comp::custom_oprf::Sender::~Sender() {
    delete (this->oprfSender);
}

Proto sparse_comp::custom_oprf::Sender::send(coproto::Socket& sock, uint_fast32_t n) {
    MC_BEGIN(Proto, this, &sock, n);

    MC_AWAIT(this->oprfSender->send(sock, n));

    MC_END();
}

void sparse_comp::custom_oprf::Sender::eval(point& point, size_t k, size_t n, VecMatrix<block>& out) {
    vector<block> point_digests(k*n);
    vector<block> rsOprfOut(k*n);

    size_t g = 0;
    for (size_t i=0;i < k;i++) {

        for (size_t j=0;j < n;j++) {
            point_digests[g] = encode_point_as_block(Sender::aes, point, i, j);
            
            g++;
        }

    }

    this->oprfSender->eval(point_digests, rsOprfOut);

    g = 0;
    for (size_t i=0;i < k;i++) {
        
        vector<block>& out_row = out[i];

        for (size_t j=0;j < n;j++) {
            out_row[j] =  rsOprfOut[g];
            g++;
        }

    }

}

void sparse_comp::custom_oprf::Sender::eval(point& point, size_t k, vector<block>& out) {
    vector<block> point_digests(k);

    for (size_t i=0;i < k;i++) {
        point_digests[i] = encode_point_as_block(Sender::aes, point, i, 0);
    }

    this->oprfSender->eval(point_digests, out);

}

sparse_comp::custom_oprf::Receiver::Receiver(MultiOprfRecvr* oprfRecvr) {
    this->oprfRecvr = oprfRecvr;
}

sparse_comp::custom_oprf::Receiver::~Receiver() {
    delete (this->oprfRecvr);
}

Proto sparse_comp::custom_oprf::Receiver::receive(coproto::Socket& sock, uint32_t n, vector<oprf_point>& points, vector<block>& outs) {
    MC_BEGIN(Proto, this, &sock, n, &points, &outs,    
    i = (size_t) 0,
    point_digests = vector<block>(n)
    );

        i = 0;
        for (const oprf_point& pot : points) {
            point_digests[i] = encode_point_as_block(Receiver::aes, pot.point, pot.sot_idx, pot.sot_choice_share);
            i++;
        }

        MC_AWAIT(this->oprfRecvr->receive(sock, point_digests, outs));
     
    MC_END();
}
