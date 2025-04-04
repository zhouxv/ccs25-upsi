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
#include <cmath>

using sparse_comp::point;

using coproto::LocalAsyncSocket;

using PRNG = osuCrypto::PRNG;
using osuCrypto::block;
using AES = osuCrypto::AES;

using std::set;
using std::unordered_map;

using macoro::sync_wait;
using macoro::when_all_ready;

using sparse_comp::spatial_hash;

static const int64_t MAX_U8_BIT_VAL = 255;
static const int64_t MAX_U32_BIT_VAL = 4294967295;

template<size_t d>
static point gen_rand_rcvr_point(PRNG& prng, uint32_t delta) {
    assert(d <= point::MAX_DIM);

    uint32_t lb = 2*delta;
    uint32_t ub = MAX_U32_BIT_VAL - 2*delta;

    uint32_t c[d];
    for (size_t i = 0; i < d; i++) {
        c[i] = (prng.get<uint32_t>() % (ub - lb + 1)) + lb;
    }

    return point(d, c);
}

template<size_t d>
static point gen_rand_sndr_point(PRNG& prng) {
    assert(d <= point::MAX_DIM);

    uint32_t c[d];
    for (size_t i = 0; i < d; i++) {
        c[i] = prng.get<uint32_t>();
    }

    return point(d, c);
}

template <size_t tr, size_t d>
static void samp_rcvr_pts(AES& aes, PRNG & prng, array<point, tr> & rcvr_points, uint8_t delta) {
    assert(d <= point::MAX_DIM);

    for (size_t i = 0; i < tr; i++) {
        rcvr_points[i] = gen_rand_rcvr_point<d>(prng, delta);
    }

}

template <typename T>
static void shuffle_vector(PRNG & prng, vector<T> & vec) {

    for (size_t i = 0; i < vec.size(); i++) {
        size_t j = prng.get<size_t>() % vec.size();
        std::swap(vec[i], vec[j]);
    }

}

template <typename T, size_t n>
static void shuffle_array(PRNG & prng, array<T, n>& arr) {

    for (size_t i = 0; i < n; i++) {
        size_t j = prng.get<size_t>() % n;
        std::swap(arr[i], arr[j]);
    }
}

// Picks n random non repeating element in the range [start,end].
static void pick_rand_from_seq(PRNG &prng, 
    size_t start, 
    size_t end, 
    size_t n, 
    vector<size_t> &out) {
    
    assert(start <= end);
    assert(n <= (end - start + 1));

    std::vector<size_t> seq;

    for (size_t i = start; i <= end; i++) {
    seq.push_back(i);
    }

    shuffle_vector(prng, seq);

    out.resize(n);

    for (size_t i = 0; i < n; i++) out[i] = seq[i];

}

/* Picks n random non repeating receiver points and values and puts them 
   in the beginning of the respective sender output arrays*/
template <size_t tr, size_t ts, size_t d, uint8_t delta>
static void pick_rand_rcvr_pts(AES& aes, PRNG& prng,
                                           array<point, tr> & rcvr_points,
                                           array<point, ts> & sndr_points,
                                           size_t n) {
    vector<size_t> rand_idxs(n);
    pick_rand_from_seq(prng, 0, tr - 1, n, rand_idxs);
   
    for (size_t i = 0; i < n; i++) {
        sndr_points[i] = rcvr_points[rand_idxs[i]];
    }
   
}

template <size_t d, uint8_t delta>
static void samp_rand_linf_matching_ball(PRNG& prng,
                                         point& pt) {
    
    for (size_t i = 0; i < d; i++) {
        uint32_t lb = std::max((int64_t) 0, ((int64_t) pt[i]) - ((int64_t) delta));
        uint32_t ub = std::min(MAX_U32_BIT_VAL, ((int64_t) pt[i]) + ((int64_t) delta));

        pt.coords[i] = (prng.get<uint32_t>() % (ub - lb + 1)) + lb;
    }

}


template <size_t d, uint8_t delta>
static void make_ball_almost_matching(PRNG& prng,
                                      point& pt) {
    int64_t int64_delta = (int64_t) delta;

    for (size_t i = 0; i < d; i++) {
        int64_t mod256_pt_i = ((int64_t) pt[i]) % (MAX_U8_BIT_VAL + 1);

        if (mod256_pt_i - delta - 1 < 0) {
            pt.coords[i] = MAX_U8_BIT_VAL;
        } else {
            pt.coords[i] = 0;
        }
        
    }

}

template <size_t tr, size_t ts, size_t d, uint8_t delta>
static void smpl_sndr_rand_pts(AES& aes, 
                                        PRNG& prng,
                                        size_t target_num_matching_pts,
                                        array<point, tr>& rcvr_points,
                                        array<point, ts>& sndr_points) {
    assert(target_num_matching_pts <= ts && target_num_matching_pts <= tr);
    assert(d <= point::MAX_DIM);
                                            
    // Copies target_num_matching_pts from rcvr_points to beginning of sndr_points.
    pick_rand_rcvr_pts<tr, ts, d, delta>(aes, prng, rcvr_points, sndr_points, target_num_matching_pts);

    // Randomizes target_num_matching_pts first points in sndr_points such that they are within delta of the corresponding points in rcvr_points.
    for (size_t i = 0; i < target_num_matching_pts; i++) {    
        samp_rand_linf_matching_ball<d,delta>(prng, sndr_points[i]);
    }


    
    // Samples the rest of the points in sndr_points such that they are not within delta of the corresponding points in rcvr_points.
    for (size_t i = target_num_matching_pts; i < ts; i++) {
        point pt = gen_rand_sndr_point<d>(prng);
        sndr_points[i] = pt;
    }

    shuffle_array<point,ts>(prng, sndr_points);
}

template <size_t tr, size_t ts, size_t d, uint8_t delta>
static void gen_constrained_rand_inputs(block seed,
                                        size_t target_num_matching_pts,
                                        array<point, tr> & rcvr_points,
                                        array<point, ts> & sndr_points) {
    assert(target_num_matching_pts <= ts && target_num_matching_pts <= tr);

    PRNG prng(seed);
    AES aes(prng.get<block>());

    samp_rcvr_pts<tr,d>(aes, prng, rcvr_points, delta);

    smpl_sndr_rand_pts<tr, ts, d, delta>(aes,
                                         prng, 
                                         target_num_matching_pts,
                                         rcvr_points, 
                                         sndr_points);
    
}

template<size_t d>
inline static bool is_linf_close(point& pt, point& ball_center, int64_t delta) {
    
    for (size_t i = 0; i < d; i++) {
        int64_t ball_center_i_mod256 = ((int64_t) ball_center[i]) % (MAX_U8_BIT_VAL + 1);
        int64_t pt_i_mod256 = ((int64_t) pt[i]) % (MAX_U8_BIT_VAL + 1);
        int64_t dist_mod256_0 = ball_center_i_mod256 - pt_i_mod256;
        int64_t dist_mod256_1 = pt_i_mod256 - ball_center_i_mod256;

        if (dist_mod256_0 < 0) dist_mod256_0 += (MAX_U8_BIT_VAL + 1);
        if (dist_mod256_1 < 0) dist_mod256_1 += (MAX_U8_BIT_VAL + 1);

        if (dist_mod256_0 > delta && dist_mod256_1 > delta) {
            return false;
        }
    }

    return true;

}

template<size_t tr, size_t ts, size_t d, uint8_t delta>
static void expected_linf_intersect(AES& aes,
                                    array<point, tr> & rcvr_points,
                                    array<point, ts> & sndr_points,
                                    set<size_t>& intersec_idxs) {
    constexpr const size_t twotod = ((size_t) pow(2, d));
    constexpr const size_t recvr_cell_count = ((size_t) pow(2, d)) * tr;
                                        
    unordered_map<block,size_t> rcvr_pts_hashes;

    std::array<point, recvr_cell_count>* rcvr_stcell_hashes = new std::array<point, recvr_cell_count>();
    std::array<point, ts>* sndr_st_hashes = new std::array<point, ts>();

    sparse_comp::spatial_cell_hash<tr, d>(rcvr_points, *rcvr_stcell_hashes, delta);
    sparse_comp::spatial_hash<ts>(sndr_points, *sndr_st_hashes, d, delta);

    for (size_t i=0;i < tr;i++) {
        for (size_t j=0;j < twotod;j++) {
            block hash = hash_point(aes, (*rcvr_stcell_hashes)[twotod*i+j]);
            rcvr_pts_hashes.insert(std::make_pair(hash,i));
        }
    }

    for (size_t i = 0; i < ts; i++) {
        block hash = hash_point(aes, (*sndr_st_hashes)[i]);
        if (rcvr_pts_hashes.contains(hash)) {
            size_t r_idx = rcvr_pts_hashes[hash];

            if (is_linf_close<d>(rcvr_points[r_idx], sndr_points[i], delta)) {
                intersec_idxs.insert(r_idx);
            }
        }
    }

    delete rcvr_stcell_hashes;
    delete sndr_st_hashes;

}

template<size_t tr, size_t ts, size_t d, uint8_t delta>
static void safer_expected_linf_intersect(AES& aes,
                                    array<point, tr> & rcvr_points,
                                    array<point, ts> & sndr_points,
                                    set<size_t>& intersec_idxs) {
    int64_t int64_delta = (int64_t) delta;

    for (size_t i = 0; i < ts; i++) {
        for (size_t j = 0; j < tr; j++) {
            bool close = true;            
            
            for (size_t k = 0; k < d; k++) {
                int64_t dist = std::abs((int64_t) sndr_points[i][k] - (int64_t) rcvr_points[j][k]);

                if (dist > int64_delta) {
                    close = false;
                    break;
                }
            }

            if (close) {
                intersec_idxs.insert(j);
            }
        }
    }
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
    c[0] = 10000;
    c[1] = 25000;
    c1[0] = 18;
    c1[1] = 25;
    c2[0] = 25;
    c2[1] = 25;

    senderPoints->at(0) = point(D, c1);
    senderPoints->at(1) = point(D, c);
    receiverPoints->at(0) = point(D, c);
    receiverPoints->at(1) = point(D, c2);

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

TEST_CASE("Fuzzy L_inf : random test (t_s=2, t_r=2, d=2, delta=10, ssp=40)","[splinf][random]")
{
    constexpr size_t TS = 2;
    constexpr size_t TR = 2;
    constexpr size_t D = 2;
    constexpr size_t DELTA = 10;
    constexpr size_t ssp = 40;
    size_t target_matching_points = 1;

    auto socks = LocalAsyncSocket::makePair();

    block seed = block(9536629126107612350ULL,2721317164341290561ULL);
    PRNG senderPRNG = PRNG(block(742130310438916676ULL, 11803924226990735076ULL));
    PRNG receiverPRNG = PRNG(block(2457938039974938056ULL, 17910068785450354990ULL));
    AES aes = AES(block(14034463513942181890ULL, 16276202269246990858ULL));

    std::array<point, TS> *senderPoints = new std::array<point, TS>();
    std::array<point, TR> *receiverPoints = new std::array<point, TR>();
    std::vector<size_t> intersec;

    gen_constrained_rand_inputs<TR, TS, D, DELTA>(seed,
                                                target_matching_points,
                                                *receiverPoints,
                                                *senderPoints);

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
    REQUIRE(intersec_set.size() == target_matching_points);

}

TEST_CASE("Fuzzy L_inf : random test (t_s=2, t_r=2, d=6, delta=10, ssp=40)","[splinf][random]")
{
    constexpr size_t TS = 2;
    constexpr size_t TR = 2;
    constexpr size_t D = 6;
    constexpr size_t DELTA = 10;
    constexpr size_t ssp = 40;
    size_t target_matching_points = 1;

    auto socks = LocalAsyncSocket::makePair();

    block seed = block(9536629126107612350ULL,2721317164341290561ULL);
    PRNG senderPRNG = PRNG(block(742130310438916676ULL, 11803924226990735076ULL));
    PRNG receiverPRNG = PRNG(block(2457938039974938056ULL, 17910068785450354990ULL));
    AES aes = AES(block(14034463513942181890ULL, 16276202269246990858ULL));

    std::array<point, TS> *senderPoints = new std::array<point, TS>();
    std::array<point, TR> *receiverPoints = new std::array<point, TR>();
    std::vector<size_t> intersec;

    gen_constrained_rand_inputs<TR, TS, D, DELTA>(seed,
                                                target_matching_points,
                                                *receiverPoints,
                                                *senderPoints);

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
    REQUIRE(intersec_set.size() == target_matching_points);

}

TEST_CASE("Fuzzy L_inf : random test (t_s=10, t_r=10, d=2, delta=10, ssp=40)","[splinf][random]")
{
    constexpr size_t TS = 10;
    constexpr size_t TR = 10;
    constexpr size_t D = 2;
    constexpr size_t DELTA = 10;
    constexpr size_t ssp = 40;
    size_t target_matching_points = 3;

    auto socks = LocalAsyncSocket::makePair();

    block seed = block(9536629026107651350ULL,2724119864341290560ULL);
    PRNG senderPRNG = PRNG(block(742130310438916676ULL, 11803924226990735076ULL));
    PRNG receiverPRNG = PRNG(block(2457938039974938056ULL, 17910068785450354990ULL));
    AES aes = AES(block(14034463513942181890ULL, 16276202269246990858ULL));

    std::array<point, TS> *senderPoints = new std::array<point, TS>();
    std::array<point, TR> *receiverPoints = new std::array<point, TR>();
    std::vector<size_t> intersec;

    gen_constrained_rand_inputs<TR, TS, D, DELTA>(seed,
                                                  target_matching_points,
                                                  *receiverPoints,
                                                  *senderPoints);

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
    REQUIRE(intersec_set.size() == target_matching_points);

}

TEST_CASE("Fuzzy L_inf : random test (t_s=10, t_r=10, d=6, delta=10, ssp=40)","[splinf][random]")
{
    constexpr size_t TS = 10;
    constexpr size_t TR = 10;
    constexpr size_t D = 6;
    constexpr size_t DELTA = 10;
    constexpr size_t ssp = 40;
    size_t target_matching_points = 3;

    auto socks = LocalAsyncSocket::makePair();

    block seed = block(9536629026107651350ULL,2724119864341290560ULL);
    PRNG senderPRNG = PRNG(block(742130310438916676ULL, 11803924226990735076ULL));
    PRNG receiverPRNG = PRNG(block(2457938039974938056ULL, 17910068785450354990ULL));
    AES aes = AES(block(14034463513942181890ULL, 16276202269246990858ULL));

    std::array<point, TS> *senderPoints = new std::array<point, TS>();
    std::array<point, TR> *receiverPoints = new std::array<point, TR>();
    std::vector<size_t> intersec;

    gen_constrained_rand_inputs<TR, TS, D, DELTA>(seed,
                                                  target_matching_points,
                                                  *receiverPoints,
                                                  *senderPoints);

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
    REQUIRE(intersec_set.size() == target_matching_points);

}

TEST_CASE("Fuzzy L_inf : safer random test (t_s=10, t_r=10, d=2, delta=10, ssp=40)","[splinf][random]")
{
    constexpr size_t TS = 10;
    constexpr size_t TR = 10;
    constexpr size_t D = 2;
    constexpr size_t DELTA = 10;
    constexpr size_t ssp = 40;
    size_t target_matching_points = 3;

    auto socks = LocalAsyncSocket::makePair();

    block seed = block(9532629066107651350ULL,2724141865341230560ULL);
    PRNG senderPRNG = PRNG(block(742130310438916676ULL, 11803924226990735076ULL));
    PRNG receiverPRNG = PRNG(block(2457938039974938056ULL, 17910068785450354990ULL));
    AES aes = AES(block(14034463513942181890ULL, 16276202269246990858ULL));

    std::array<point, TS> *senderPoints = new std::array<point, TS>();
    std::array<point, TR> *receiverPoints = new std::array<point, TR>();
    std::vector<size_t> intersec;

    gen_constrained_rand_inputs<TR, TS, D, DELTA>(seed,
                                                  target_matching_points,
                                                  *receiverPoints,
                                                  *senderPoints);

    sparse_comp::fuzzy_linf::Sender<TR, TS, D, DELTA, ssp> fuzzyLinfSender(senderPRNG, aes);
    sparse_comp::fuzzy_linf::Receiver<TS, TR, D, DELTA, ssp> fuzzyLinfRecvr(receiverPRNG, aes);

    auto sender_proto = fuzzyLinfSender.send(socks[0], *senderPoints);
    auto receiver_proto = fuzzyLinfRecvr.receive(socks[1], *receiverPoints, intersec);

    sync_wait(when_all_ready(std::move(sender_proto), std::move(receiver_proto)));

    set<size_t> expected_intersec;

    safer_expected_linf_intersect<TR, TS, D, DELTA>(aes,
                                    *receiverPoints,
                                    *senderPoints,
                                    expected_intersec);

    set<size_t> intersec_set(intersec.begin(), intersec.end());

    delete senderPoints;
    delete receiverPoints;

    REQUIRE(intersec_set == expected_intersec);
    REQUIRE(intersec_set.size() == target_matching_points);

}

TEST_CASE("Fuzzy L_inf : safer random test (t_s=10, t_r=10, d=6, delta=10, ssp=40)","[splinf][random]")
{
    constexpr size_t TS = 10;
    constexpr size_t TR = 10;
    constexpr size_t D = 6;
    constexpr size_t DELTA = 10;
    constexpr size_t ssp = 40;
    size_t target_matching_points = 3;

    auto socks = LocalAsyncSocket::makePair();

    block seed = block(9532629066107651350ULL,2724141865341230560ULL);
    PRNG senderPRNG = PRNG(block(742130310438916676ULL, 11803924226990735076ULL));
    PRNG receiverPRNG = PRNG(block(2457938039974938056ULL, 17910068785450354990ULL));
    AES aes = AES(block(14034463513942181890ULL, 16276202269246990858ULL));

    std::array<point, TS> *senderPoints = new std::array<point, TS>();
    std::array<point, TR> *receiverPoints = new std::array<point, TR>();
    std::vector<size_t> intersec;

    gen_constrained_rand_inputs<TR, TS, D, DELTA>(seed,
                                                  target_matching_points,
                                                  *receiverPoints,
                                                  *senderPoints);

    sparse_comp::fuzzy_linf::Sender<TR, TS, D, DELTA, ssp> fuzzyLinfSender(senderPRNG, aes);
    sparse_comp::fuzzy_linf::Receiver<TS, TR, D, DELTA, ssp> fuzzyLinfRecvr(receiverPRNG, aes);

    auto sender_proto = fuzzyLinfSender.send(socks[0], *senderPoints);
    auto receiver_proto = fuzzyLinfRecvr.receive(socks[1], *receiverPoints, intersec);

    sync_wait(when_all_ready(std::move(sender_proto), std::move(receiver_proto)));

    set<size_t> expected_intersec;

    safer_expected_linf_intersect<TR, TS, D, DELTA>(aes,
                                    *receiverPoints,
                                    *senderPoints,
                                    expected_intersec);

    set<size_t> intersec_set(intersec.begin(), intersec.end());

    delete senderPoints;
    delete receiverPoints;

    REQUIRE(intersec_set == expected_intersec);
    REQUIRE(intersec_set.size() == target_matching_points);

}

TEST_CASE("Fuzzy L_inf : safer random test (t_s=256, t_r=256, d=2, delta=10, ssp=40)","[splinf][random]")
{
    constexpr size_t TS = 256;
    constexpr size_t TR = 256;
    constexpr size_t D = 2;
    constexpr size_t DELTA = 10;
    constexpr size_t ssp = 40;
    size_t target_matching_points = 3;

    auto socks = LocalAsyncSocket::makePair();

    block seed = block(4337051091801073793ULL,11168696276776517519ULL);
    PRNG senderPRNG = PRNG(block(742130310438916676ULL, 11803924226990735076ULL));
    PRNG receiverPRNG = PRNG(block(2457938039974938056ULL, 17910068785450354990ULL));
    AES aes = AES(block(14034463513942181890ULL, 16276202269246990858ULL));

    std::array<point, TS> *senderPoints = new std::array<point, TS>();
    std::array<point, TR> *receiverPoints = new std::array<point, TR>();
    std::vector<size_t> intersec;

    gen_constrained_rand_inputs<TR, TS, D, DELTA>(seed,
                                                  target_matching_points,
                                                  *receiverPoints,
                                                  *senderPoints);

    sparse_comp::fuzzy_linf::Sender<TR, TS, D, DELTA, ssp> fuzzyLinfSender(senderPRNG, aes);
    sparse_comp::fuzzy_linf::Receiver<TS, TR, D, DELTA, ssp> fuzzyLinfRecvr(receiverPRNG, aes);

    auto sender_proto = fuzzyLinfSender.send(socks[0], *senderPoints);
    auto receiver_proto = fuzzyLinfRecvr.receive(socks[1], *receiverPoints, intersec);

    sync_wait(when_all_ready(std::move(sender_proto), std::move(receiver_proto)));

    set<size_t> expected_intersec;

    safer_expected_linf_intersect<TR, TS, D, DELTA>(aes,
                                    *receiverPoints,
                                    *senderPoints,
                                    expected_intersec);

    set<size_t> intersec_set(intersec.begin(), intersec.end());

    delete senderPoints;
    delete receiverPoints;

    REQUIRE(intersec_set == expected_intersec);
    REQUIRE(intersec_set.size() == target_matching_points);

}