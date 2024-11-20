#include "cryptoTools/Crypto/PRNG.h"
#include "cryptoTools/Common/block.h"
#include "CustomOPRF/CustomizedOPRF.h"
#include "./SpBSOT/SpBSOT.h"
#include "./Common/Common.h"
#include "./Common/VecMatrix.h"
#include "./Common/ZN.h"
#include "./Common/HashUtils.h"
#include "volePSI/Paxos.h"
#include "./SpLInf/SpLInf.h"
//#include "./SpL1/SpL1.h"
#include "./MultiOPRF/MultiOPRF.h"
#include "./Common/SockUtils.h"
#include <cstdint>
#include <iostream> 
#include <vector>
#include <array>
#include <random>
#include <limits>
#include <chrono>
#include <algorithm>
#include <set>
#include <limits>
#include <math.h>

//#define TR 1024
//#define TS 256
#define TR 33554432
#define TS 32768
#define DELTA 10
#define D 10
//#define D 2
#define ssp 40

using namespace sparse_comp;
using std::chrono::high_resolution_clock;
using std::chrono::duration_cast;
using std::chrono::duration;
using std::chrono::milliseconds;
using std::min;
using std::max;

using OprfSender = sparse_comp::custom_oprf::Sender;
using OprfReceiver = sparse_comp::custom_oprf::Receiver;

using sparse_comp::VecMatrix;
using sparse_comp::ZN;
using sparse_comp::hash_point;

using PRNG = osuCrypto::PRNG;

using coproto::LocalAsyncSocket;

using macoro::sync_wait;
using macoro::when_all_ready;

template <typename T>
using vector = std::vector<T>;

template <typename T, size_t N>
using array = std::array<T, N>;


template <size_t tr, size_t ts, size_t d, u_int8_t delta>
void gen_random_party_inputs(array<point,tr>& rcvr_sparse_points,
                            array<array<uint32_t,d>,tr>& rcvr_in_values,
                            array<point,ts>& sndr_sparse_points,
                            array<array<uint32_t,d>,ts>& sndr_in_values) {
    AES aes = AES(block(13133210048402866,17132091720387928));
    std::set<block> s;
    std::pair<std::set<block>::iterator,bool> ret;
    std::random_device dev;
    std::mt19937 rng(dev());
    std::uniform_int_distribution<std::mt19937::result_type> dst(0, std::numeric_limits<uint32_t>::max()); // distribution in range [1, 6]
    uint32_t coords[point::MAX_DIM] = {};

    size_t i=0;
    while (s.size() < tr) {
        for (size_t j=0;j < d;j++) {
            coords[j] = dst(rng);
            rcvr_in_values[i][j] = coords[j] >> 24;
        }

        rcvr_sparse_points[i] = point(d,coords);

        ret = s.insert(hash_point(aes, point(d,coords)));

        if (ret.second)
            i++;
    }

    s.clear();
    i=0;
    while (s.size() < ts) {
        for (size_t j=0;j < d;j++) {
            coords[j] = dst(rng);
            sndr_in_values[i][j] = coords[j] >> 24;
        }
        
        sndr_sparse_points[i] = point(d,coords);

        ret = s.insert(hash_point(aes, point(d,coords)));

        if (ret.second)
            i++;
    }

}

static Proto sendVec(coproto::Socket& sock, vector<block>& m) {

     MC_BEGIN(Proto, &sock, &m);
    
        MC_AWAIT(sock.send(m));
    
    MC_END();
    
}

static Proto recvVec(coproto::Socket& sock, vector<block>& m) {
     MC_BEGIN(Proto, &sock, &m);
    
        MC_AWAIT(sock.recvResize(m));

     MC_END();
    
}
/*
int main(int argc, char** argv) {

    auto socks = LocalAsyncSocket::makePair();

    std::vector<uint64_t> v_recv(1000);
    std::vector<uint64_t> v_send(1000);
    for (size_t i = 0; i < v_send.size(); i++) {
        v_send[i] = i;
    }

    auto send_task = sparse_comp::send<uint64_t,17>(socks[0], v_send);
    auto recv_task = sparse_comp::receive<uint64_t,17>(socks[1], v_recv.size(), v_recv);

    coproto::sync_wait(coproto::when_all_ready(std::move(send_task),std::move(recv_task)));


    for (size_t i = 0; i < v_recv.size(); i++) {
        //std::cout << v_recv[i] << " = " << i << std::endl;
    }

    return 0;
}*/

/*
int main(int argc, char** argv) {

    const size_t ms = 2*2*2*2*26214400;

    auto socks = LocalAsyncSocket::makePair();

    auto senderPRNG = PRNG(block(5,6));
    auto receiverPRNG = PRNG(block(37,44));
    auto aes = AES(block(311,127));

    vector<block>* m = new vector<block>(ms); 
    vector<block>* m2 = new vector<block>(); 

    for (size_t i=0;i < ms;i++) {
        m->at(i) = receiverPRNG.get<block>();
    }

    auto sender_proto = sendVec(socks[0], *m);
    auto receiver_proto = recvVec(socks[1], *m2);

    sync_wait(when_all_ready(std::move(sender_proto),std::move(receiver_proto)));

    delete m;
    delete m2;

    return 0;
}*/

/*
int main(int argc, char** argv) {

    auto socks = LocalAsyncSocket::makePair();

    auto senderPRNG = PRNG(block(5,6));
    auto receiverPRNG = PRNG(block(37,44));
    auto aes = AES(block(311,127));

    std::array<point,TS>* senderSparsePoints = new std::array<point,TS>();
    std::array<point,TR>* receiverSparsePoints = new std::array<point,TR>();
    array<array<uint32_t,D>,TS>* sender_in_values = new array<array<uint32_t,D>,TS>();
    array<array<uint32_t,D>,TR>* receiver_in_values = new array<array<uint32_t,D>,TR>();

    gen_random_party_inputs<TR,TS,D,DELTA>(
                                            *receiverSparsePoints,
                                            *receiver_in_values,
                                            *senderSparsePoints,
                                            *sender_in_values);






    vector<size_t> intersec;

    auto t1 = high_resolution_clock::now();

    sparse_comp::sp_l1::Sender<TR,TS,D,DELTA,ssp> spL1Sender(senderPRNG, aes);
    sparse_comp::sp_l1::Receiver<TS,TR,D,DELTA,ssp> spL1Recvr(receiverPRNG, aes);

    auto sender_proto = spL1Sender.send(socks[0], *senderSparsePoints, *sender_in_values);
    auto receiver_proto = spL1Recvr.receive(socks[1], *receiverSparsePoints, *receiver_in_values, intersec);

    sync_wait(when_all_ready(std::move(sender_proto),std::move(receiver_proto)));

    auto t2 = high_resolution_clock::now();

    auto ms_int = duration_cast<milliseconds>(t2 - t1);

    std::cout << TR << ":" << TS << ":" << D << ":" << DELTA << ";";
    std::cout << ms_int.count() << ";";
    std::cout << (((double)(socks[0].bytesSent()+socks[0].bytesReceived()))/1024.0/1024.0); 

    delete senderSparsePoints; 
    delete receiverSparsePoints;
    delete sender_in_values;
    delete receiver_in_values;

    return 0;
}

*/

int main(int argc, char** argv) {

    auto socks = LocalAsyncSocket::makePair();

    auto senderPRNG = PRNG(block(5,6));
    auto receiverPRNG = PRNG(block(37,44));
    auto aes = AES(block(311,127));

    std::array<point,TS>* senderSparsePoints = new std::array<point,TS>();
    std::array<point,TR>* receiverSparsePoints = new std::array<point,TR>();
    array<array<uint32_t,D>,TS>* sender_in_values = new array<array<uint32_t,D>,TS>();
    array<array<uint32_t,D>,TR>* receiver_in_values = new array<array<uint32_t,D>,TR>();

/*
    uint32_t c[point::MAX_DIM] = {};
    uint32_t c1[point::MAX_DIM] = {};
    uint32_t c2[point::MAX_DIM] = {};
    c[0] = 4;
    c[1] = 5;
    c1[0] = 10;
    c1[1] = 5;
    c2[0] = 68;
    c2[1] = 44;

    senderSparsePoints->at(0) = point(2,c1);
    senderSparsePoints->at(1) = point(2,c);
    receiverSparsePoints->at(0) = point(2,c);
    receiverSparsePoints->at(1) = point(2,c1);

    sender_in_values->at(0)[0] = 77;
    sender_in_values->at(0)[1] = 32;
    sender_in_values->at(1)[0] = 50;
    sender_in_values->at(1)[1] = 66;

    receiver_in_values->at(0)[0] = 59;
    receiver_in_values->at(0)[1] = 56;
    receiver_in_values->at(1)[0] = 77;
    receiver_in_values->at(1)[1] = 32;*/
    

    gen_random_party_inputs<TR,TS,D,DELTA>(
                                            *receiverSparsePoints,
                                            *receiver_in_values,
                                            *senderSparsePoints,
                                            *sender_in_values);


    vector<size_t> intersec;

    auto t1 = high_resolution_clock::now();

    sparse_comp::sp_linf::Sender<TR,TS,D,DELTA,ssp> spLinfSender(senderPRNG, aes);
    sparse_comp::sp_linf::Receiver<TS,TR,D,DELTA,ssp> spLinfRecvr(receiverPRNG, aes);

    auto sender_proto = spLinfSender.send(socks[0], *senderSparsePoints, *sender_in_values);
    auto receiver_proto = spLinfRecvr.receive(socks[1], *receiverSparsePoints, *receiver_in_values, intersec);

    sync_wait(when_all_ready(std::move(sender_proto),std::move(receiver_proto)));

    auto t2 = high_resolution_clock::now();

    auto ms_int = duration_cast<milliseconds>(t2 - t1);

    //std::cout << ms_int.count() << "ms" << std::endl; 

    //std::cout << "intersec size: " << intersec.size() << std::endl;

    /*for (size_t i=0;i < intersec.size();i++) {
        std::cout << intersec[i] << std::endl;
    }*/

    std::cout << TR << ":" << TS << ":" << D << ":" << DELTA << ";";
    std::cout << ms_int.count() << ";";
    std::cout << (((double)(socks[0].bytesSent()+socks[0].bytesReceived()))/1024.0/1024.0); 

    delete senderSparsePoints; 
    delete receiverSparsePoints;
    delete sender_in_values;
    delete receiver_in_values;

  return 0;
}

