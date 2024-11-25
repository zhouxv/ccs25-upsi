#include "catch2/catch_test_macros.hpp"
#include "cryptoTools/Crypto/PRNG.h"
#include "cryptoTools/Common/block.h"
#include "coproto/Socket/LocalAsyncSock.h"
#include "cryptoTools/Crypto/AES.h"
#include "../sparseComp/Common/Common.h"
#include "../sparseComp/Common/HashUtils.h"
#include "../sparseComp/SpLInf/SpLInf.h"
#include <cstdint>
#include <chrono>
#include <utility>
#include <vector>
#include <algorithm>
#include <set>
#include <unordered_map>

using sparse_comp::point;

using coproto::LocalAsyncSocket;

using PRNG = osuCrypto::PRNG;
using osuCrypto::block;
using AES = osuCrypto::AES;

using std::chrono::high_resolution_clock;
using std::chrono::milliseconds;
using sparse_comp::hash_point;

using std::array;
using std::set;
using std::unordered_map;

using macoro::sync_wait;
using macoro::when_all_ready;

static const int64_t MAX_8_BIT_VAL = 255;


template <typename T>
static void shuffle_vector(PRNG & prng, vector<T> & vec) {

    for (size_t i = 0; i < vec.size(); i++) {
        size_t j = prng.get<size_t>() % vec.size();
        std::swap(vec[i], vec[j]);
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
static void pick_rand_rcvr_pts_and_vals(AES& aes, PRNG& prng,
                                        array<point, tr> & rcvr_sparse_points,
                                        array<array<uint32_t, d>, tr> & rcvr_in_values,
                                        array<point, ts> & sndr_sparse_points,
                                        array<array<uint32_t, d>, ts> & sndr_in_values,
                                        size_t n) {
    vector<size_t> rand_idxs(n);
    pick_rand_from_seq(prng, 0, tr - 1, n, rand_idxs);

    for (size_t i = 0; i < n; i++) {
        sndr_sparse_points[i] = rcvr_sparse_points[rand_idxs[i]];
        sndr_in_values[i] = rcvr_in_values[rand_idxs[i]];
    }

}

template <size_t d, uint8_t delta>
static array<uint32_t,d> samp_rand_linf_matching_ball(PRNG& prng,
                                                      array<uint32_t,d>& pt) {
    array<uint32_t,d> ball;
    
    for (size_t i = 0; i < d; i++) {
        uint32_t lb = std::max((int64_t) delta, ((int64_t) pt[i]) - (int64_t) delta);
        uint32_t ub = std::min(MAX_8_BIT_VAL-delta, ((int64_t) pt[i]) + (int64_t) delta);

        ball[i] = (prng.get<uint32_t>() % (ub - lb + 1)) + lb;
    }

    return ball;
}

template <typename T, size_t n>
static void shuffle_array(PRNG & prng, array<T, n>& arr) {

    for (size_t i = 0; i < n; i++) {
        size_t j = prng.get<size_t>() % n;
        std::swap(arr[i], arr[j]);
    }
}

template <typename T1, typename T2, size_t n>
static void shuffle_array_pair(PRNG & prng, array<T1, n>& arr1, array<T2, n>& arr2) {

    for (size_t i = 0; i < n; i++) {
        size_t j = prng.get<size_t>() % n;
        std::swap(arr1[i], arr1[j]);
        std::swap(arr2[i], arr2[j]);
    }

}

template<size_t d>
static point gen_rand_point(PRNG& prng) {
    assert(d <= point::MAX_DIM);

    uint32_t c[d];
    for (size_t i = 0; i < d; i++) {
        c[i] = prng.get<uint32_t>();
    }

    return point(d, c);
}

template <size_t t, size_t d>
static void gen_rand_array_mtx(PRNG & prng, array<array<uint32_t, d>, t> & arr)
{

    for (size_t i = 0; i < t; i++)
    {
        for (size_t j = 0; j < d; j++)
        {
            arr[i][j] = prng.get<uint32_t>();
        }
    }
}

template <size_t n, size_t d>
static void gen_rand_point_array(PRNG& prng, array<point, n>& pt_arr) {

    for (size_t i = 0; i < n; i++) {
        pt_arr[i] = gen_rand_point<d>(prng);
    }
}

template <size_t tr, size_t d>
static void samp_rcvr_in_vals(PRNG & prng, array<array<uint32_t, d>, tr> & rcvr_in_vals)
{

    for (size_t i = 0; i < tr; i++)
    {
        for (size_t j = 0; j < d; j++)
        {
            rcvr_in_vals[i][j] = prng.get<uint8_t>();
        }
    }
}

template <size_t d, uint8_t delta>
static array<uint32_t,d> smpl_sndr_rand_val(PRNG& prng) {
    array<uint32_t,d> pt;
    uint32_t lb = delta;
    uint32_t ub = MAX_8_BIT_VAL - delta;
    
    for (size_t i = 0; i < d; i++) {
        pt[i] = (prng.get<uint32_t>() % (ub - lb + 1)) + lb;
        //std::cout << "pt[" << i << "]: " << pt[i] << std::endl;
    }

    return pt;
}

template <size_t tr, size_t d>
static void samp_rcvr_sparse_pts(AES& aes, PRNG & prng, array<point, tr> & rcvr_sparse_points) {

    set<block> existing_pts;

    for (size_t i = 0; i < tr; i++) {
        point pt = gen_rand_point<d>(prng);
        block hash = hash_point(aes, pt);

        while (existing_pts.contains(hash)) {
            pt = gen_rand_point<d>(prng);
            hash = hash_point(aes, pt);
        }

        rcvr_sparse_points[i] = pt;
        existing_pts.insert(hash);
    }

}

template <size_t tr, size_t ts, size_t d, uint8_t delta>
static void smpl_sndr_rand_pts_and_vals(AES& aes, 
                                        PRNG& prng,
                                        size_t min_num_matching_bins,
                                        size_t min_num_matching_pts,
                                        array<point, tr> & rcvr_sparse_points,
                                        array<array<uint32_t, d>, tr> & rcvr_in_values,
                                        array<point, ts>& sndr_sparse_points,
                                        array<array<uint32_t, d>, ts> & sndr_in_values) {
    // Copies min_num_matching_random bins from rcvr_sparse_points to beginning of sndr_sparse_points.
    pick_rand_rcvr_pts_and_vals<tr, ts, d, delta>(aes, prng, rcvr_sparse_points, rcvr_in_values, sndr_sparse_points, sndr_in_values, min_num_matching_bins);

    // Randomizes min_num_matching_pts first points in sndr_in_values such that they are within delta of the corresponding points in rcvr_in_values.
    for (size_t i = 0; i < min_num_matching_pts; i++) {
        sndr_in_values[i] = samp_rand_linf_matching_ball<d,delta>(prng, sndr_in_values[i]);
    }

    // Randomizes the remaining points in sndr_in_values copied picked randomly from rcvr_in_values.
    for (size_t i = min_num_matching_pts; i < min_num_matching_bins; i++) {
        sndr_in_values[i] = smpl_sndr_rand_val<d,delta>(prng);
    }

    set<block> existing_pts;

    // Stores the hash of the first min_num_matching_bins points in sndr_sparse_points.
    for (size_t i=0;i < min_num_matching_bins;i++) {
        existing_pts.insert(hash_point(aes, sndr_sparse_points[i]));
    }

    // Genrates remaning sender points randomly and stores them in sndr_sparse_points such that no repeating points happen.
    size_t i = min_num_matching_bins;
    while (i < ts) {
        point pt = gen_rand_point<d>(prng);
        block hash = hash_point(aes, pt);

        if (!existing_pts.contains(hash)) {
            sndr_sparse_points[i] = pt;
            existing_pts.insert(hash);
            i++;
        }
    }

    // Generates random values for the remaining values in sndr_in_values such that they respect the delta value.
    for (size_t i = min_num_matching_bins; i < ts; i++) {
        sndr_in_values[i] = smpl_sndr_rand_val<d,delta>(prng);
    }

    for (size_t i = 0; i < ts; i++) {
        for (size_t j = 0; j < d; j++) {
            if (!(sndr_in_values.at(i)[j] >= delta && sndr_in_values.at(i)[j] <= MAX_8_BIT_VAL - delta)) {
                std::cout << "sndr_in_values[" << i << "][" << j << "]: " << sndr_in_values.at(i)[j] << std::endl;
            }
        }
    }

    shuffle_array_pair<point,array<uint32_t,d>,ts>(prng, sndr_sparse_points, sndr_in_values);

}

template <size_t tr, size_t ts, size_t d, uint8_t delta>
static void gen_constrained_rand_inputs(block seed,
                                        size_t min_num_matching_bins,
                                        size_t min_num_matching_pts,
                                        array<point, tr> & rcvr_sparse_points,
                                        array<array<uint32_t, d>, tr> & rcvr_in_values,
                                        array<point, ts> & sndr_sparse_points,
                                        array<array<uint32_t, d>, ts> & sndr_in_values) {
    size_t min_t = std::min(tr, ts);
    min_num_matching_bins = std::min(min_num_matching_bins, min_t);
    min_num_matching_pts = std::min(min_num_matching_pts, min_t);

    assert(min_num_matching_pts <= min_num_matching_bins);

    PRNG prng(seed);
    AES aes(prng.get<block>());

    samp_rcvr_sparse_pts<tr,d>(aes, prng, rcvr_sparse_points);
    samp_rcvr_in_vals<tr,d>(prng, rcvr_in_values);

    smpl_sndr_rand_pts_and_vals<tr, ts, d, delta>(aes,
                                                  prng, 
                                                  min_num_matching_bins, 
                                                  min_num_matching_pts, 
                                                  rcvr_sparse_points, 
                                                  rcvr_in_values, 
                                                  sndr_sparse_points, 
                                                  sndr_in_values);
    
}

template<size_t d>
static bool is_linf_close(array<uint32_t, d>& pt, array<uint32_t, d>& ball_center, uint8_t delta) {
    for (size_t i = 0; i < d; i++) {
        //std::cout << "pt[" << i << "]: " << pt[i] << " ball_center[" << i << "]: " << ball_center[i] << std::endl;
        assert(ball_center[i] >= delta && ball_center[i] <= MAX_8_BIT_VAL - delta);
        if (!(ball_center[i] - delta <= pt[i] && pt[i] <= ball_center[i] + delta)) {
            return false;
        }
    }

    return true;

}

template<size_t tr, size_t ts, size_t d, uint8_t delta>
static void expected_linf_intersect(AES& aes,
                                    array<point, tr> & rcvr_sparse_points,
                                    array<array<uint32_t, d>, tr> & rcvr_in_values,
                                    array<point, ts> & sndr_sparse_points,
                                    array<array<uint32_t, d>, ts> & sndr_in_values,
                                    set<size_t>& intersec_idxs) {
    unordered_map<block,size_t> rcvr_pts_hashes;

    for (size_t i = 0; i < tr; i++) {
        rcvr_pts_hashes.insert(std::make_pair(hash_point(aes, rcvr_sparse_points[i]),i));
    }

    for (size_t i = 0; i < ts; i++) {
        block hash = hash_point(aes, sndr_sparse_points[i]);
        if (rcvr_pts_hashes.contains(hash)) {
            size_t r_idx = rcvr_pts_hashes[hash];

            if (is_linf_close<d>(rcvr_in_values[r_idx], sndr_in_values[i], delta)) {
                intersec_idxs.insert(r_idx);
            }
        }
    }

}

TEST_CASE("Sparse L_inf : simple test (t_s=2, t_r=2, d=2, delta=10, ssp=40)")
{
    constexpr size_t TS = 2;
    constexpr size_t TR = 2;
    constexpr size_t D = 2;
    constexpr size_t DELTA = 10;
    constexpr size_t ssp = 40;

    auto socks = LocalAsyncSocket::makePair();

    PRNG senderPRNG = PRNG(block(5, 6));
    PRNG receiverPRNG = PRNG(block(37, 44));
    AES aes = AES(block(311, 127));

    std::array<point, TS> *senderSparsePoints = new std::array<point, TS>();
    std::array<point, TR> *receiverSparsePoints = new std::array<point, TR>();
    array<array<uint32_t, D>, TS> *sender_in_values = new array<array<uint32_t, D>, TS>();
    array<array<uint32_t, D>, TR> *receiver_in_values = new array<array<uint32_t, D>, TR>();

    uint32_t c[point::MAX_DIM];
    uint32_t c1[point::MAX_DIM];
    uint32_t c2[point::MAX_DIM];
    c[0] = 4;
    c[1] = 5;
    c1[0] = 10;
    c1[1] = 5;
    c2[0] = 68;
    c2[1] = 44;

    senderSparsePoints->at(0) = point(D, c1);
    senderSparsePoints->at(1) = point(D, c);
    receiverSparsePoints->at(0) = point(D, c);
    receiverSparsePoints->at(1) = point(D, c1);

    sender_in_values->at(0)[0] = 77;
    sender_in_values->at(0)[1] = 32;
    sender_in_values->at(1)[0] = 50;
    sender_in_values->at(1)[1] = 66;

    receiver_in_values->at(0)[0] = 59;
    receiver_in_values->at(0)[1] = 56;
    receiver_in_values->at(1)[0] = 77;
    receiver_in_values->at(1)[1] = 32;

    vector<size_t> intersec;

    sparse_comp::sp_linf::Sender<TR, TS, D, DELTA, ssp> spLinfSender(senderPRNG, aes);
    sparse_comp::sp_linf::Receiver<TS, TR, D, DELTA, ssp> spLinfRecvr(receiverPRNG, aes);

    auto sender_proto = spLinfSender.send(socks[0], *senderSparsePoints, *sender_in_values);
    auto receiver_proto = spLinfRecvr.receive(socks[1], *receiverSparsePoints, *receiver_in_values, intersec);

    sync_wait(when_all_ready(std::move(sender_proto), std::move(receiver_proto)));

    set<size_t> expected_intersec;

    expected_linf_intersect<TR, TS, D, DELTA>(aes, 
                                              *receiverSparsePoints, 
                                              *receiver_in_values, 
                                              *senderSparsePoints, 
                                              *sender_in_values, 
                                              expected_intersec);

    set<size_t> intersec_set(intersec.begin(), intersec.end());

    delete senderSparsePoints;
    delete receiverSparsePoints;
    delete sender_in_values;
    delete receiver_in_values;

    REQUIRE(intersec_set == expected_intersec);
}

TEST_CASE("Sparse L_inf : random constrained test (t_s=2, t_r=2, d=2, delta=2, ssp=40)") {
    constexpr size_t TS = 2;
    constexpr size_t TR = 2;
    constexpr size_t D = 2;
    constexpr size_t DELTA = 2;
    constexpr size_t ssp = 40;
    size_t min_num_matching_bins = 1;
    size_t min_num_matching_pts = 1;

    auto socks = LocalAsyncSocket::makePair();
    block seed = block(9536629026107651350ULL,2724119864341290560ULL);
    PRNG senderPRNG = PRNG(block(742130310438916676ULL, 11803924226990735076ULL));
    PRNG receiverPRNG = PRNG(block(2457938039974938056ULL, 17910068785450354990ULL));
    AES aes = AES(block(14034463513942181890ULL, 16276202269246990858ULL));

    std::array<point, TS> *senderSparsePoints = new std::array<point, TS>();
    std::array<point, TR> *receiverSparsePoints = new std::array<point, TR>();
    array<array<uint32_t, D>, TS> *sender_in_values = new array<array<uint32_t, D>, TS>();
    array<array<uint32_t, D>, TR> *receiver_in_values = new array<array<uint32_t, D>, TR>();
    vector<size_t> intersec;

    gen_constrained_rand_inputs<TR, TS, D, DELTA>(seed, 
                                                  min_num_matching_bins,
                                                  min_num_matching_pts,
                                                  *receiverSparsePoints,
                                                  *receiver_in_values, 
                                                  *senderSparsePoints, 
                                                  *sender_in_values);

    sparse_comp::sp_linf::Sender<TR, TS, D, DELTA, ssp> spLinfSender(senderPRNG, aes);
    sparse_comp::sp_linf::Receiver<TS, TR, D, DELTA, ssp> spLinfRecvr(receiverPRNG, aes);

    auto sender_proto = spLinfSender.send(socks[0], *senderSparsePoints, *sender_in_values);
    auto receiver_proto = spLinfRecvr.receive(socks[1], *receiverSparsePoints, *receiver_in_values, intersec);

    sync_wait(when_all_ready(std::move(sender_proto), std::move(receiver_proto)));

    set<size_t> expected_intersec;

    expected_linf_intersect<TR, TS, D, DELTA>(aes, 
                                              *receiverSparsePoints, 
                                              *receiver_in_values, 
                                              *senderSparsePoints, 
                                              *sender_in_values, 
                                              expected_intersec);
    REQUIRE(expected_intersec.size() >= min_num_matching_pts);


    set<size_t> intersec_set(intersec.begin(), intersec.end());

    delete senderSparsePoints;
    delete receiverSparsePoints;
    delete sender_in_values;
    delete receiver_in_values;

    REQUIRE(intersec_set == expected_intersec);
}

TEST_CASE("Sparse L_inf : random constrained test (t_s=10, t_r=10, d=2, delta=10, ssp=40)") {
    constexpr size_t TS = 10;
    constexpr size_t TR = 10;
    constexpr size_t D = 2;
    constexpr size_t DELTA = 10;
    constexpr size_t ssp = 40;
    size_t min_num_matching_bins = 4;
    size_t min_num_matching_pts = 3;

    auto socks = LocalAsyncSocket::makePair();
    block seed = block(15118332354506195568ULL,16271682051025781703ULL);
    PRNG senderPRNG = PRNG(block(6173834459723358941ULL, 13868447840116429460ULL));
    PRNG receiverPRNG = PRNG(block(15049811300244079322ULL, 15605034700940240023ULL));
    AES aes = AES(block(14034463513942181890ULL, 16276202269246990858ULL));

    std::array<point, TS> *senderSparsePoints = new std::array<point, TS>();
    std::array<point, TR> *receiverSparsePoints = new std::array<point, TR>();
    array<array<uint32_t, D>, TS> *sender_in_values = new array<array<uint32_t, D>, TS>();
    array<array<uint32_t, D>, TR> *receiver_in_values = new array<array<uint32_t, D>, TR>();
    vector<size_t> intersec;

    gen_constrained_rand_inputs<TR, TS, D, DELTA>(seed, 
                                                  min_num_matching_bins,
                                                  min_num_matching_pts,
                                                  *receiverSparsePoints,
                                                  *receiver_in_values, 
                                                  *senderSparsePoints, 
                                                  *sender_in_values);

    sparse_comp::sp_linf::Sender<TR, TS, D, DELTA, ssp> spLinfSender(senderPRNG, aes);
    sparse_comp::sp_linf::Receiver<TS, TR, D, DELTA, ssp> spLinfRecvr(receiverPRNG, aes);

    auto sender_proto = spLinfSender.send(socks[0], *senderSparsePoints, *sender_in_values);
    auto receiver_proto = spLinfRecvr.receive(socks[1], *receiverSparsePoints, *receiver_in_values, intersec);

    sync_wait(when_all_ready(std::move(sender_proto), std::move(receiver_proto)));

    set<size_t> expected_intersec;

    expected_linf_intersect<TR, TS, D, DELTA>(aes, 
                                              *receiverSparsePoints, 
                                              *receiver_in_values, 
                                              *senderSparsePoints, 
                                              *sender_in_values, 
                                              expected_intersec);
    REQUIRE(expected_intersec.size() >= min_num_matching_pts);


    set<size_t> intersec_set(intersec.begin(), intersec.end());

    delete senderSparsePoints;
    delete receiverSparsePoints;
    delete sender_in_values;
    delete receiver_in_values;

    REQUIRE(intersec_set == expected_intersec);
}


TEST_CASE("Sparse L_inf : random constrained test (t_s=256, t_r=1024, d=2, delta=10, ssp=40)") {
    constexpr size_t TS = 256;
    constexpr size_t TR = 1024;
    constexpr size_t D = 2;
    constexpr size_t DELTA = 10;
    constexpr size_t ssp = 40;
    size_t min_num_matching_bins = 53;
    size_t min_num_matching_pts = 17;

    auto socks = LocalAsyncSocket::makePair();
    block seed = block(9536629026107651350ULL,2724119864341290560ULL);
    PRNG senderPRNG = PRNG(block(15914074867899273501ULL, 6004108516319388444ULL));
    PRNG receiverPRNG = PRNG(block(6427781726132732903ULL, 8471345356057289138ULL));
    AES aes = AES(block(14034463513942181890ULL, 16276202269246990858ULL));

    std::array<point, TS> *senderSparsePoints = new std::array<point, TS>();
    std::array<point, TR> *receiverSparsePoints = new std::array<point, TR>();
    array<array<uint32_t, D>, TS> *sender_in_values = new array<array<uint32_t, D>, TS>();
    array<array<uint32_t, D>, TR> *receiver_in_values = new array<array<uint32_t, D>, TR>();
    vector<size_t> intersec;

    gen_constrained_rand_inputs<TR, TS, D, DELTA>(seed, 
                                                  min_num_matching_bins,
                                                  min_num_matching_pts,
                                                  *receiverSparsePoints,
                                                  *receiver_in_values, 
                                                  *senderSparsePoints, 
                                                  *sender_in_values);

    sparse_comp::sp_linf::Sender<TR, TS, D, DELTA, ssp> spLinfSender(senderPRNG, aes);
    sparse_comp::sp_linf::Receiver<TS, TR, D, DELTA, ssp> spLinfRecvr(receiverPRNG, aes);

    auto sender_proto = spLinfSender.send(socks[0], *senderSparsePoints, *sender_in_values);
    auto receiver_proto = spLinfRecvr.receive(socks[1], *receiverSparsePoints, *receiver_in_values, intersec);

    sync_wait(when_all_ready(std::move(sender_proto), std::move(receiver_proto)));

    set<size_t> expected_intersec;
    expected_linf_intersect<TR, TS, D, DELTA>(aes, 
                                              *receiverSparsePoints, 
                                              *receiver_in_values, 
                                              *senderSparsePoints, 
                                              *sender_in_values, 
                                              expected_intersec);
    REQUIRE(expected_intersec.size() >= min_num_matching_pts);

    set<size_t> intersec_set(intersec.begin(), intersec.end());

    delete senderSparsePoints;
    delete receiverSparsePoints;
    delete sender_in_values;
    delete receiver_in_values;

    REQUIRE(intersec_set == expected_intersec);
}