#include "catch2/catch_test_macros.hpp"
#include "cryptoTools/Crypto/PRNG.h"
#include "cryptoTools/Common/block.h"
#include "coproto/Socket/LocalAsyncSock.h"
#include "cryptoTools/Crypto/AES.h"
#include "../sparseComp/Common/Common.h"
#include "../sparseComp/SpLInf/SpLInf.h"
#include <cstdint>
#include <chrono>
#include <utility>
#include <vector>

using sparse_comp::point;

using coproto::LocalAsyncSocket;

using PRNG = osuCrypto::PRNG;
using osuCrypto::block;
using AES = osuCrypto::AES;

using std::chrono::high_resolution_clock;
using std::chrono::milliseconds;

int add(int a, int b)
{
    return a + b;
}

TEST_CASE("Sparse L_inf : simple test (t_s=2, t_r=2, d=2, delta=2, ssp=40)") {
    constexpr size_t TS = 2;
    constexpr size_t TR = 2;
    constexpr size_t D = 2;
    constexpr size_t DELTA = 2;
    constexpr size_t ssp = 40;

    auto socks = LocalAsyncSocket::makePair();

    PRNG senderPRNG = PRNG(block(5,6));
    PRNG receiverPRNG = PRNG(block(37,44));
    AES aes = AES(block(311,127));

    std::array<point,TS>* senderSparsePoints = new std::array<point,TS>();
    std::array<point,TR>* receiverSparsePoints = new std::array<point,TR>();
    array<array<uint32_t,D>,TS>* sender_in_values = new array<array<uint32_t,D>,TS>();
    array<array<uint32_t,D>,TR>* receiver_in_values = new array<array<uint32_t,D>,TR>();

    uint32_t c[point::MAX_DIM];
    uint32_t c1[point::MAX_DIM];
    uint32_t c2[point::MAX_DIM];
    c[0] = 4;
    c[1] = 5;
    c1[0] = 10;
    c1[1] = 5;
    c2[0] = 68;
    c2[1] = 44;

    senderSparsePoints->at(0) = point(D,c1);
    senderSparsePoints->at(1) = point(D,c);
    receiverSparsePoints->at(0) = point(D,c);
    receiverSparsePoints->at(1) = point(D,c1);

    sender_in_values->at(0)[0] = 77;
    sender_in_values->at(0)[1] = 32;
    sender_in_values->at(1)[0] = 50;
    sender_in_values->at(1)[1] = 66;

    receiver_in_values->at(0)[0] = 59;
    receiver_in_values->at(0)[1] = 56;
    receiver_in_values->at(1)[0] = 77;
    receiver_in_values->at(1)[1] = 32;

    vector<size_t> intersec;

    sparse_comp::sp_linf::Sender<TR,TS,D,DELTA,ssp> spLinfSender(senderPRNG, aes);
    sparse_comp::sp_linf::Receiver<TS,TR,D,DELTA,ssp> spLinfRecvr(receiverPRNG, aes);

    auto sender_proto = spLinfSender.send(socks[0], *senderSparsePoints, *sender_in_values);
    auto receiver_proto = spLinfRecvr.receive(socks[1], *receiverSparsePoints, *receiver_in_values, intersec);

    sync_wait(when_all_ready(std::move(sender_proto),std::move(receiver_proto)));

    delete senderSparsePoints;
    delete receiverSparsePoints;
    delete sender_in_values;
    delete receiver_in_values;

    REQUIRE(intersec.size() == 1);
   
}