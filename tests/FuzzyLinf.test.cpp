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
#include <unordered_map>
#include <utility>


using sparse_comp::point;

using coproto::LocalAsyncSocket;

using PRNG = osuCrypto::PRNG;
using osuCrypto::block;
using AES = osuCrypto::AES;

using std::set;
using std::unordered_map;

using macoro::sync_wait;
using macoro::when_all_ready;

static const int64_t MAX_8_BIT_VAL = 255;

template<size_t tr, size_t ts, size_t d, uint8_t delta>
static void expected_linf_intersect(std::array<point, tr> & rcvr_points,
                                    std::array<point, ts> & sndr_points,
                                    set<size_t>& intersec_idxs) {
    
}

template<size_t d>
inline static bool is_linf_close(point& pt, point& ball_center, uint8_t delta) {
    int64_t int64_delta = (int64_t) delta;
    
    for (size_t i = 0; i < d; i++) {
        // std::cout << "pt[" << i << "]: " << pt[i] << " ball_center[" << i << "]: " << ball_center[i] << std::endl;
        // assert(ball_center[i] >= delta && ball_center[i] <= MAX_8_BIT_VAL - delta);
        int64_t ub = (((int64_t) ball_center.coords[i]) % (MAX_8_BIT_VAL + 1)) + int64_delta;
        int64_t lb = (((int64_t) ball_center.coords[i]) % (MAX_8_BIT_VAL + 1)) - int64_delta;
        int64_t pt_i = ((int64_t) pt.coords[i]) % (MAX_8_BIT_VAL + 1);
        
        if (!(lb <= pt_i && pt_i <= ub)) {
            return false;
        }
    }

    return true;

}

template<size_t tr, size_t ts, size_t d, uint8_t delta>
static void expected_linf_intersect(AES& aes,
                                    array<point, tr> & rcvr_sparse_points,
                                    array<point, ts> & sndr_sparse_points,
                                    set<size_t>& intersec_idxs) {
    unordered_map<block,size_t> rcvr_pts_hashes;

    std::array<point, tr>* rcvr_st_hashes = new std::array<point, tr>();
    std::array<point, ts>* sndr_st_hashes = new std::array<point, ts>();

    sparse_comp::spatial_hash<tr>(rcvr_sparse_points, *rcvr_st_hashes, d, delta);
    sparse_comp::spatial_hash<ts>(sndr_sparse_points, *sndr_st_hashes, d, delta);

    for (size_t i = 0; i < tr; i++) {
        rcvr_pts_hashes.insert(std::make_pair(hash_point(aes, (*rcvr_st_hashes)[i]),i));
    }

    for (size_t i = 0; i < ts; i++) {
        block hash = hash_point(aes, (*sndr_st_hashes)[i]);
        if (rcvr_pts_hashes.contains(hash)) {
            size_t r_idx = rcvr_pts_hashes[hash];

            if (is_linf_close<d>(rcvr_sparse_points[r_idx], sndr_sparse_points[i], delta)) {
                intersec_idxs.insert(r_idx);
            }
        }
    }

    delete rcvr_st_hashes;
    delete sndr_st_hashes;

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

    set<size_t> expected_intersec;

    expected_linf_intersect<TR, TS, D, DELTA>(aes,
                                    *receiverPoints,
                                    *senderPoints,
                                    expected_intersec);

    set<size_t> intersec_set(intersec.begin(), intersec.end());


    delete senderPoints;
    delete receiverPoints;

    REQUIRE(intersec_set == expected_intersec);

}