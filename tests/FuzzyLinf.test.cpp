#include "catch2/catch_test_macros.hpp"
#include "cryptoTools/Crypto/PRNG.h"
#include "cryptoTools/Common/block.h"
#include "coproto/Socket/LocalAsyncSock.h"
#include "cryptoTools/Crypto/AES.h"
#include "../sparseComp/Common/Common.h"
#include "../sparseComp/Common/HashUtils.h"
#include "../sparseComp/FuzzyLinf/FuzzyLinf.h"
#include <cstdint>
#include <array>
#include <vector>
#include <set>

using sparse_comp::point;

using coproto::LocalAsyncSocket;

using PRNG = osuCrypto::PRNG;
using osuCrypto::block;
using AES = osuCrypto::AES;

using std::set;

using macoro::sync_wait;
using macoro::when_all_ready;

template<size_t tr, size_t ts, size_t d, uint8_t delta>
static void expected_linf_intersect(std::array<point, tr> & rcvr_points,
                                    std::array<point, ts> & sndr_points,
                                    set<size_t>& intersec_idxs) {
    
}

TEST_CASE("Fuzzy L_inf : simple test (t_s=2, t_r=2, d=2, delta=10, ssp=40)","[splinf][simple]")
{
    constexpr size_t TS = 2;
    constexpr size_t TR = 2;
    constexpr size_t D = 2;
    constexpr size_t DELTA = 10;
    constexpr size_t ssp = 40;

    auto socks = LocalAsyncSocket::makePair();

    PRNG senderPRNG = PRNG(block(50, 6));
    PRNG receiverPRNG = PRNG(block(37, 44));
    AES aes = AES(block(311, 127));

    std::array<point, TS> *senderPoints = new std::array<point, TS>();
    std::array<point, TR> *receiverPoints = new std::array<point, TR>();
    std::vector<size_t> intersec;

    uint32_t c[point::MAX_DIM];
    uint32_t c1[point::MAX_DIM];
    uint32_t c2[point::MAX_DIM];
    c[0] = 4;
    c[1] = 25;
    c1[0] = 10;
    c1[1] = 5;
    c2[0] = 68;
    c2[1] = 44;

    senderPoints->at(0) = point(D, c1);
    senderPoints->at(1) = point(D, c);
    receiverPoints->at(0) = point(D, c);
    receiverPoints->at(1) = point(D, c1);

    sparse_comp::fuzzy_linf::Sender<TR, TS, D, DELTA, ssp> fuzzyLinfSender(senderPRNG, aes);
    sparse_comp::fuzzy_linf::Receiver<TS, TR, D, DELTA, ssp> fuzzyLinfRecvr(receiverPRNG, aes);

    auto sender_proto = fuzzyLinfSender.send(socks[0], *senderPoints);
    auto receiver_proto = fuzzyLinfRecvr.receive(socks[1], *receiverPoints, intersec);

    sync_wait(when_all_ready(std::move(sender_proto), std::move(receiver_proto)));

    std::cout << "Intersection set size: " << intersec.size() << std::endl;

    for (const auto& idx : intersec) {
        std::cout << "Intersection index: " << idx << std::endl;
    }

    set<size_t> expected_intersec;

    // COMPUTE EXPECTED INTERSECTION [TODO]

    set<size_t> intersec_set(intersec.begin(), intersec.end());


    delete senderPoints;
    delete receiverPoints;

    //REQUIRE(intersec_set == expected_intersec);

}