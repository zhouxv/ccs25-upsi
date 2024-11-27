#include "catch2/catch_test_macros.hpp"
#include "catch2/benchmark/catch_benchmark.hpp"
#include "cryptoTools/Crypto/PRNG.h"
#include "cryptoTools/Common/block.h"
#include "coproto/Socket/LocalAsyncSock.h"
#include "cryptoTools/Crypto/AES.h"
#include "../sparseComp/Common/Common.h"
#include "../sparseComp/Common/HashUtils.h"
#include "../sparseComp/SpL1/SpL1.h"
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

template <typename T>
static void shuffle_vector(PRNG & prng, vector<T> & vec) {

    for (size_t i = 0; i < vec.size(); i++) {
        size_t j = prng.get<size_t>() % vec.size();
        std::swap(vec[i], vec[j]);
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

template<size_t d>
static point gen_rand_point(PRNG& prng) {
    assert(d <= point::MAX_DIM);

    uint32_t c[d];
    for (size_t i = 0; i < d; i++) {
        c[i] = prng.get<uint32_t>();
    }

    return point(d, c);
}

static const int64_t MAX_8_BIT_VAL = 255;
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
static array<uint32_t,d> samp_rand_linf_matching_ball(PRNG& prng,
                                                      array<uint32_t,d>& pt) {
    array<uint32_t,d> ball;
    array<uint32_t,d> offsets;
    std::fill(std::begin(offsets), std::end(offsets), 0);

    for (size_t i=0;i < delta;i++) {
        size_t j = prng.get<size_t>() % d;
        offsets[j]++;
    }

    for (size_t i = 0; i < d; i++) {
        bool orient = prng.get<bool>();
        orient = (orient && pt[i] + offsets[i] <= MAX_8_BIT_VAL) ? orient : false;
        orient = (!orient && (MAX_8_BIT_VAL + pt[i] - offsets[i]) >= MAX_8_BIT_VAL) ? orient : true;

        ball[i] = orient ? pt[i] + offsets[i] : pt[i] - offsets[i];
    }

    return ball;
}

template <size_t d, uint8_t delta>
static array<uint32_t,d> smpl_sndr_rand_val(PRNG& prng) {
    array<uint32_t,d> pt;
    
    for (size_t i = 0; i < d; i++) {
        pt[i] = prng.get<uint8_t>();
    }

    return pt;
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
static bool is_l1_close(array<uint32_t, d>& pt,
                          array<uint32_t, d>& ball_center, 
                          uint8_t delta) {
    uint32_t acc_delta = 0;
    
    for (size_t i = 0; i < d; i++) {
        const uint32_t l = std::min(pt[i] , ball_center[i]);
        const uint32_t h = std::max(pt[i] , ball_center[i]);

        acc_delta += h - l;
    }

    return acc_delta <= delta;

}

template<size_t tr, size_t ts, size_t d, uint8_t delta>
static void expected_l1_intersect(AES& aes,
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

            if (is_l1_close<d>(rcvr_in_values[r_idx], sndr_in_values[i], delta)) {
                intersec_idxs.insert(r_idx);
            }
        }
    }

}

// START OF TESTS FOR N=M=2^8

TEST_CASE("spl1 (t_s=256 t_r=1024 d=2 delta=10 ssp=40)","[spl1][n=m=2^8]") {

    BENCHMARK_ADVANCED("t_s=256 t_r=1024 d=2 delta=10 ssp=40")(Catch::Benchmark::Chronometer meter) {
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

        sparse_comp::sp_l1::Sender<TR,TS,D,DELTA,ssp> spL1Sender(senderPRNG, aes);
        sparse_comp::sp_l1::Receiver<TS,TR,D,DELTA,ssp> spL1Recvr(receiverPRNG, aes);

        auto sender_proto = spL1Sender.send(socks[0], *senderSparsePoints, *sender_in_values);
        auto receiver_proto = spL1Recvr.receive(socks[1], *receiverSparsePoints, *receiver_in_values, intersec);

        meter.measure([&sender_proto,&receiver_proto]() { sync_wait(when_all_ready(std::move(sender_proto), std::move(receiver_proto))); });

        set<size_t> expected_intersec;
        expected_l1_intersect<TR, TS, D, DELTA>(aes, 
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

        const double nMBsExchanged = ((double)(socks[0].bytesSent()+socks[0].bytesReceived()))/1024.0/1024.0; 

        SUCCEED("Number of MBs exchanged: " << nMBsExchanged);
    };
}

TEST_CASE("spl1 (t_s=256 t_r=16384 d=6 delta=10 ssp=40)","[spl1][n=m=2^8]") {


    BENCHMARK_ADVANCED("t_s=256, t_r=16384, d=6, delta=10, ssp=40")(Catch::Benchmark::Chronometer meter) {
        constexpr size_t TS = 256;
        constexpr size_t TR = 16384;
        constexpr size_t D = 6;
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

        sparse_comp::sp_l1::Sender<TR,TS,D,DELTA,ssp> spL1Sender(senderPRNG, aes);
        sparse_comp::sp_l1::Receiver<TS,TR,D,DELTA,ssp> spL1Recvr(receiverPRNG, aes);

        auto sender_proto = spL1Sender.send(socks[0], *senderSparsePoints, *sender_in_values);
        auto receiver_proto = spL1Recvr.receive(socks[1], *receiverSparsePoints, *receiver_in_values, intersec);

        meter.measure([&sender_proto,&receiver_proto]() { sync_wait(when_all_ready(std::move(sender_proto), std::move(receiver_proto))); });

        set<size_t> expected_intersec;
        expected_l1_intersect<TR, TS, D, DELTA>(aes, 
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

        const double nMBsExchanged = ((double)(socks[0].bytesSent()+socks[0].bytesReceived()))/1024.0/1024.0; 

        SUCCEED("Number of MBs exchanged: " << nMBsExchanged);
    };
}

TEST_CASE("spl1 (t_s=256 t_r=262144 d=10 delta=10 ssp=40)","[spl1][n=m=2^8]") {

    BENCHMARK_ADVANCED("t_s=256, t_r=262144, d=10, delta=10, ssp=40")(Catch::Benchmark::Chronometer meter) {
        constexpr size_t TS = 256;
        constexpr size_t TR = 262144;
        constexpr size_t D = 10;
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

        sparse_comp::sp_l1::Sender<TR,TS,D,DELTA,ssp> spL1Sender(senderPRNG, aes);
        sparse_comp::sp_l1::Receiver<TS,TR,D,DELTA,ssp> spL1Recvr(receiverPRNG, aes);

        auto sender_proto = spL1Sender.send(socks[0], *senderSparsePoints, *sender_in_values);
        auto receiver_proto = spL1Recvr.receive(socks[1], *receiverSparsePoints, *receiver_in_values, intersec);

        meter.measure([&sender_proto,&receiver_proto]() { sync_wait(when_all_ready(std::move(sender_proto), std::move(receiver_proto))); });

        set<size_t> expected_intersec;
        expected_l1_intersect<TR, TS, D, DELTA>(aes, 
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

        const double nMBsExchanged = ((double)(socks[0].bytesSent()+socks[0].bytesReceived()))/1024.0/1024.0; 

        SUCCEED("Number of MBs exchanged: " << nMBsExchanged);
    };
}

TEST_CASE("spl1 (t_s=256 t_r=1024 d=2 delta=30 ssp=40)","[spl1][n=m=2^8]") {

    BENCHMARK_ADVANCED("t_s=256, t_r=1024, d=2, delta=30, ssp=40")(Catch::Benchmark::Chronometer meter) {
        constexpr size_t TS = 256;
        constexpr size_t TR = 1024;
        constexpr size_t D = 2;
        constexpr size_t DELTA = 30;
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

        sparse_comp::sp_l1::Sender<TR,TS,D,DELTA,ssp> spL1Sender(senderPRNG, aes);
        sparse_comp::sp_l1::Receiver<TS,TR,D,DELTA,ssp> spL1Recvr(receiverPRNG, aes);

        auto sender_proto = spL1Sender.send(socks[0], *senderSparsePoints, *sender_in_values);
        auto receiver_proto = spL1Recvr.receive(socks[1], *receiverSparsePoints, *receiver_in_values, intersec);

        meter.measure([&sender_proto,&receiver_proto]() { sync_wait(when_all_ready(std::move(sender_proto), std::move(receiver_proto))); });

        set<size_t> expected_intersec;
        expected_l1_intersect<TR, TS, D, DELTA>(aes, 
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

        const double nMBsExchanged = ((double)(socks[0].bytesSent()+socks[0].bytesReceived()))/1024.0/1024.0; 

        SUCCEED("Number of MBs exchanged: " << nMBsExchanged);
    };
}

TEST_CASE("spl1 (t_s=256 t_r=16384 d=6 delta=30 ssp=40)","[spl1][n=m=2^8]") {


    BENCHMARK_ADVANCED("t_s=256, t_r=16384, d=6, delta=30, ssp=40")(Catch::Benchmark::Chronometer meter) {
        constexpr size_t TS = 256;
        constexpr size_t TR = 16384;
        constexpr size_t D = 6;
        constexpr size_t DELTA = 30;
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

        sparse_comp::sp_l1::Sender<TR,TS,D,DELTA,ssp> spL1Sender(senderPRNG, aes);
        sparse_comp::sp_l1::Receiver<TS,TR,D,DELTA,ssp> spL1Recvr(receiverPRNG, aes);

        auto sender_proto = spL1Sender.send(socks[0], *senderSparsePoints, *sender_in_values);
        auto receiver_proto = spL1Recvr.receive(socks[1], *receiverSparsePoints, *receiver_in_values, intersec);

        meter.measure([&sender_proto,&receiver_proto]() { sync_wait(when_all_ready(std::move(sender_proto), std::move(receiver_proto))); });

        set<size_t> expected_intersec;
        expected_l1_intersect<TR, TS, D, DELTA>(aes, 
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

        const double nMBsExchanged = ((double)(socks[0].bytesSent()+socks[0].bytesReceived()))/1024.0/1024.0; 

        SUCCEED("Number of MBs exchanged: " << nMBsExchanged);
    };

}

TEST_CASE("spl1 (t_s=256 t_r=262144 d=10 delta=30 ssp=40)","[spl1][n=m=2^8]") {

    BENCHMARK_ADVANCED("t_s=256, t_r=262144, d=10, delta=30, ssp=40")(Catch::Benchmark::Chronometer meter) {
        constexpr size_t TS = 256;
        constexpr size_t TR = 262144;
        constexpr size_t D = 10;
        constexpr size_t DELTA = 30;
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

        sparse_comp::sp_l1::Sender<TR,TS,D,DELTA,ssp> spL1Sender(senderPRNG, aes);
        sparse_comp::sp_l1::Receiver<TS,TR,D,DELTA,ssp> spL1Recvr(receiverPRNG, aes);

        auto sender_proto = spL1Sender.send(socks[0], *senderSparsePoints, *sender_in_values);
        auto receiver_proto = spL1Recvr.receive(socks[1], *receiverSparsePoints, *receiver_in_values, intersec);

        meter.measure([&sender_proto,&receiver_proto]() { sync_wait(when_all_ready(std::move(sender_proto), std::move(receiver_proto))); });

        set<size_t> expected_intersec;
        expected_l1_intersect<TR, TS, D, DELTA>(aes, 
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

        const double nMBsExchanged = ((double)(socks[0].bytesSent()+socks[0].bytesReceived()))/1024.0/1024.0; 

        SUCCEED("Number of MBs exchanged: " << nMBsExchanged);
    };

}

// END OF TESTS FOR N=M=2^8

// START OF TESTS FOR N=M=2^12

TEST_CASE("spl1 (t_s=4096 t_r=16384 d=2 delta=10 ssp=40)","[spl1][n=m=2^12]") {

    BENCHMARK_ADVANCED("t_s=4096 t_r=16384 d=2 delta=10 ssp=40")(Catch::Benchmark::Chronometer meter) {
        constexpr size_t TS = 4096;
        constexpr size_t TR = 16384;
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

        sparse_comp::sp_l1::Sender<TR,TS,D,DELTA,ssp> spL1Sender(senderPRNG, aes);
        sparse_comp::sp_l1::Receiver<TS,TR,D,DELTA,ssp> spL1Recvr(receiverPRNG, aes);

        auto sender_proto = spL1Sender.send(socks[0], *senderSparsePoints, *sender_in_values);
        auto receiver_proto = spL1Recvr.receive(socks[1], *receiverSparsePoints, *receiver_in_values, intersec);

        meter.measure([&sender_proto,&receiver_proto]() { sync_wait(when_all_ready(std::move(sender_proto), std::move(receiver_proto))); });

        set<size_t> expected_intersec;
        expected_l1_intersect<TR, TS, D, DELTA>(aes, 
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

        const double nMBsExchanged = ((double)(socks[0].bytesSent()+socks[0].bytesReceived()))/1024.0/1024.0; 

        SUCCEED("Number of MBs exchanged: " << nMBsExchanged);
    };
}

TEST_CASE("spl1 (t_s=4096 t_r=262144 d=6 delta=10 ssp=40)","[spl1][n=m=2^12]") {


    BENCHMARK_ADVANCED("t_s=4096, t_r=262144, d=6, delta=10, ssp=40")(Catch::Benchmark::Chronometer meter) {
        constexpr size_t TS = 4096;
        constexpr size_t TR = 262144;
        constexpr size_t D = 6;
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

        sparse_comp::sp_l1::Sender<TR,TS,D,DELTA,ssp> spL1Sender(senderPRNG, aes);
        sparse_comp::sp_l1::Receiver<TS,TR,D,DELTA,ssp> spL1Recvr(receiverPRNG, aes);

        auto sender_proto = spL1Sender.send(socks[0], *senderSparsePoints, *sender_in_values);
        auto receiver_proto = spL1Recvr.receive(socks[1], *receiverSparsePoints, *receiver_in_values, intersec);

        meter.measure([&sender_proto,&receiver_proto]() { sync_wait(when_all_ready(std::move(sender_proto), std::move(receiver_proto))); });

        set<size_t> expected_intersec;
        expected_l1_intersect<TR, TS, D, DELTA>(aes, 
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

        const double nMBsExchanged = ((double)(socks[0].bytesSent()+socks[0].bytesReceived()))/1024.0/1024.0; 

        SUCCEED("Number of MBs exchanged: " << nMBsExchanged);
    };
}

TEST_CASE("spl1 (t_s=4096 t_r=4194304 d=10 delta=10 ssp=40)","[spl1][n=m=2^12]") {

    BENCHMARK_ADVANCED("t_s=4096, t_r=4194304, d=10, delta=10, ssp=40")(Catch::Benchmark::Chronometer meter) {
        constexpr size_t TS = 4096;
        constexpr size_t TR = 4194304;
        constexpr size_t D = 10;
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

        sparse_comp::sp_l1::Sender<TR,TS,D,DELTA,ssp> spL1Sender(senderPRNG, aes);
        sparse_comp::sp_l1::Receiver<TS,TR,D,DELTA,ssp> spL1Recvr(receiverPRNG, aes);

        auto sender_proto = spL1Sender.send(socks[0], *senderSparsePoints, *sender_in_values);
        auto receiver_proto = spL1Recvr.receive(socks[1], *receiverSparsePoints, *receiver_in_values, intersec);

        meter.measure([&sender_proto,&receiver_proto]() { sync_wait(when_all_ready(std::move(sender_proto), std::move(receiver_proto))); });

        set<size_t> expected_intersec;
        expected_l1_intersect<TR, TS, D, DELTA>(aes, 
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

        const double nMBsExchanged = ((double)(socks[0].bytesSent()+socks[0].bytesReceived()))/1024.0/1024.0; 

        SUCCEED("Number of MBs exchanged: " << nMBsExchanged);
    };
}

TEST_CASE("spl1 (t_s=4096 t_r=16384 d=2 delta=30 ssp=40)","[spl1][n=m=2^12]") {

    BENCHMARK_ADVANCED("t_s=4096, t_r=16384, d=2, delta=30, ssp=40")(Catch::Benchmark::Chronometer meter) {
        constexpr size_t TS = 4096;
        constexpr size_t TR = 16384;
        constexpr size_t D = 2;
        constexpr size_t DELTA = 30;
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

        sparse_comp::sp_l1::Sender<TR,TS,D,DELTA,ssp> spL1Sender(senderPRNG, aes);
        sparse_comp::sp_l1::Receiver<TS,TR,D,DELTA,ssp> spL1Recvr(receiverPRNG, aes);

        auto sender_proto = spL1Sender.send(socks[0], *senderSparsePoints, *sender_in_values);
        auto receiver_proto = spL1Recvr.receive(socks[1], *receiverSparsePoints, *receiver_in_values, intersec);

        meter.measure([&sender_proto,&receiver_proto]() { sync_wait(when_all_ready(std::move(sender_proto), std::move(receiver_proto))); });

        set<size_t> expected_intersec;
        expected_l1_intersect<TR, TS, D, DELTA>(aes, 
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

        const double nMBsExchanged = ((double)(socks[0].bytesSent()+socks[0].bytesReceived()))/1024.0/1024.0; 

        SUCCEED("Number of MBs exchanged: " << nMBsExchanged);
    };
}

TEST_CASE("spl1 (t_s=4096 t_r=262144 d=6 delta=30 ssp=40)","[spl1][n=m=2^12]") {


    BENCHMARK_ADVANCED("t_s=4096, t_r=262144, d=6, delta=30, ssp=40")(Catch::Benchmark::Chronometer meter) {
        constexpr size_t TS = 4096;
        constexpr size_t TR = 262144;
        constexpr size_t D = 6;
        constexpr size_t DELTA = 30;
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

        sparse_comp::sp_l1::Sender<TR,TS,D,DELTA,ssp> spL1Sender(senderPRNG, aes);
        sparse_comp::sp_l1::Receiver<TS,TR,D,DELTA,ssp> spL1Recvr(receiverPRNG, aes);

        auto sender_proto = spL1Sender.send(socks[0], *senderSparsePoints, *sender_in_values);
        auto receiver_proto = spL1Recvr.receive(socks[1], *receiverSparsePoints, *receiver_in_values, intersec);

        meter.measure([&sender_proto,&receiver_proto]() { sync_wait(when_all_ready(std::move(sender_proto), std::move(receiver_proto))); });

        set<size_t> expected_intersec;
        expected_l1_intersect<TR, TS, D, DELTA>(aes, 
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

        const double nMBsExchanged = ((double)(socks[0].bytesSent()+socks[0].bytesReceived()))/1024.0/1024.0; 

        SUCCEED("Number of MBs exchanged: " << nMBsExchanged);
    };

}

TEST_CASE("spl1 (t_s=4096 t_r=4194304 d=10 delta=30 ssp=40)","[spl1][n=m=2^12]") {

    BENCHMARK_ADVANCED("t_s=4096, t_r=4194304, d=10, delta=30, ssp=40")(Catch::Benchmark::Chronometer meter) {
        constexpr size_t TS = 4096;
        constexpr size_t TR = 4194304;
        constexpr size_t D = 10;
        constexpr size_t DELTA = 30;
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

        sparse_comp::sp_l1::Sender<TR,TS,D,DELTA,ssp> spL1Sender(senderPRNG, aes);
        sparse_comp::sp_l1::Receiver<TS,TR,D,DELTA,ssp> spL1Recvr(receiverPRNG, aes);

        auto sender_proto = spL1Sender.send(socks[0], *senderSparsePoints, *sender_in_values);
        auto receiver_proto = spL1Recvr.receive(socks[1], *receiverSparsePoints, *receiver_in_values, intersec);

        meter.measure([&sender_proto,&receiver_proto]() { sync_wait(when_all_ready(std::move(sender_proto), std::move(receiver_proto))); });

        set<size_t> expected_intersec;
        expected_l1_intersect<TR, TS, D, DELTA>(aes, 
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

        const double nMBsExchanged = ((double)(socks[0].bytesSent()+socks[0].bytesReceived()))/1024.0/1024.0; 

        SUCCEED("Number of MBs exchanged: " << nMBsExchanged);
    };

}

// END OF TESTS FOR N=M=2^12

// START OF TESTS FOR N=M=2^16

TEST_CASE("spl1 (t_s=65536 t_r=262144 d=2 delta=10 ssp=40)","[spl1][n=m=2^16]") {

    BENCHMARK_ADVANCED("t_s=65536 t_r=262144 d=2 delta=10 ssp=40")(Catch::Benchmark::Chronometer meter) {
        constexpr size_t TS = 65536;
        constexpr size_t TR = 262144;
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

        sparse_comp::sp_l1::Sender<TR,TS,D,DELTA,ssp> spL1Sender(senderPRNG, aes);
        sparse_comp::sp_l1::Receiver<TS,TR,D,DELTA,ssp> spL1Recvr(receiverPRNG, aes);

        auto sender_proto = spL1Sender.send(socks[0], *senderSparsePoints, *sender_in_values);
        auto receiver_proto = spL1Recvr.receive(socks[1], *receiverSparsePoints, *receiver_in_values, intersec);

        meter.measure([&sender_proto,&receiver_proto]() { sync_wait(when_all_ready(std::move(sender_proto), std::move(receiver_proto))); });

        set<size_t> expected_intersec;
        expected_l1_intersect<TR, TS, D, DELTA>(aes, 
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

        const double nMBsExchanged = ((double)(socks[0].bytesSent()+socks[0].bytesReceived()))/1024.0/1024.0; 

        SUCCEED("Number of MBs exchanged: " << nMBsExchanged);
    };
}

TEST_CASE("spl1 (t_s=65536 t_r=4194304 d=6 delta=10 ssp=40)","[spl1][n=m=2^16]") {


    BENCHMARK_ADVANCED("t_s=65536, t_r=4194304, d=6, delta=10, ssp=40")(Catch::Benchmark::Chronometer meter) {
        constexpr size_t TS = 65536;
        constexpr size_t TR = 4194304;
        constexpr size_t D = 6;
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

        sparse_comp::sp_l1::Sender<TR,TS,D,DELTA,ssp> spL1Sender(senderPRNG, aes);
        sparse_comp::sp_l1::Receiver<TS,TR,D,DELTA,ssp> spL1Recvr(receiverPRNG, aes);

        auto sender_proto = spL1Sender.send(socks[0], *senderSparsePoints, *sender_in_values);
        auto receiver_proto = spL1Recvr.receive(socks[1], *receiverSparsePoints, *receiver_in_values, intersec);

        meter.measure([&sender_proto,&receiver_proto]() { sync_wait(when_all_ready(std::move(sender_proto), std::move(receiver_proto))); });

        set<size_t> expected_intersec;
        expected_l1_intersect<TR, TS, D, DELTA>(aes, 
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

        const double nMBsExchanged = ((double)(socks[0].bytesSent()+socks[0].bytesReceived()))/1024.0/1024.0; 

        SUCCEED("Number of MBs exchanged: " << nMBsExchanged);
    };
}

TEST_CASE("spl1 (t_s=65536 t_r=67108864 d=10 delta=10 ssp=40)","[spl1][n=m=2^16]") {

    BENCHMARK_ADVANCED("t_s=65536, t_r=67108864, d=10, delta=10, ssp=40")(Catch::Benchmark::Chronometer meter) {
        constexpr size_t TS = 65536;
        constexpr size_t TR = 67108864;
        constexpr size_t D = 10;
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

        sparse_comp::sp_l1::Sender<TR,TS,D,DELTA,ssp> spL1Sender(senderPRNG, aes);
        sparse_comp::sp_l1::Receiver<TS,TR,D,DELTA,ssp> spL1Recvr(receiverPRNG, aes);

        auto sender_proto = spL1Sender.send(socks[0], *senderSparsePoints, *sender_in_values);
        auto receiver_proto = spL1Recvr.receive(socks[1], *receiverSparsePoints, *receiver_in_values, intersec);

        meter.measure([&sender_proto,&receiver_proto]() { sync_wait(when_all_ready(std::move(sender_proto), std::move(receiver_proto))); });

        set<size_t> expected_intersec;
        expected_l1_intersect<TR, TS, D, DELTA>(aes, 
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

        const double nMBsExchanged = ((double)(socks[0].bytesSent()+socks[0].bytesReceived()))/1024.0/1024.0; 

        SUCCEED("Number of MBs exchanged: " << nMBsExchanged);
    };
}

TEST_CASE("spl1 (t_s=65536 t_r=262144 d=2 delta=30 ssp=40)","[spl1][n=m=2^16]") {

    BENCHMARK_ADVANCED("t_s=65536, t_r=262144, d=2, delta=30, ssp=40")(Catch::Benchmark::Chronometer meter) {
        constexpr size_t TS = 65536;
        constexpr size_t TR = 262144;
        constexpr size_t D = 2;
        constexpr size_t DELTA = 30;
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

        sparse_comp::sp_l1::Sender<TR,TS,D,DELTA,ssp> spL1Sender(senderPRNG, aes);
        sparse_comp::sp_l1::Receiver<TS,TR,D,DELTA,ssp> spL1Recvr(receiverPRNG, aes);

        auto sender_proto = spL1Sender.send(socks[0], *senderSparsePoints, *sender_in_values);
        auto receiver_proto = spL1Recvr.receive(socks[1], *receiverSparsePoints, *receiver_in_values, intersec);

        meter.measure([&sender_proto,&receiver_proto]() { sync_wait(when_all_ready(std::move(sender_proto), std::move(receiver_proto))); });

        set<size_t> expected_intersec;
        expected_l1_intersect<TR, TS, D, DELTA>(aes, 
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

        const double nMBsExchanged = ((double)(socks[0].bytesSent()+socks[0].bytesReceived()))/1024.0/1024.0; 

        SUCCEED("Number of MBs exchanged: " << nMBsExchanged);
    };
}

TEST_CASE("spl1 (t_s=65536 t_r=4194304 d=6 delta=30 ssp=40)","[spl1][n=m=2^16]") {


    BENCHMARK_ADVANCED("t_s=65536, t_r=4194304, d=6, delta=30, ssp=40")(Catch::Benchmark::Chronometer meter) {
        constexpr size_t TS = 65536;
        constexpr size_t TR = 4194304;
        constexpr size_t D = 6;
        constexpr size_t DELTA = 30;
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

        sparse_comp::sp_l1::Sender<TR,TS,D,DELTA,ssp> spL1Sender(senderPRNG, aes);
        sparse_comp::sp_l1::Receiver<TS,TR,D,DELTA,ssp> spL1Recvr(receiverPRNG, aes);

        auto sender_proto = spL1Sender.send(socks[0], *senderSparsePoints, *sender_in_values);
        auto receiver_proto = spL1Recvr.receive(socks[1], *receiverSparsePoints, *receiver_in_values, intersec);

        meter.measure([&sender_proto,&receiver_proto]() { sync_wait(when_all_ready(std::move(sender_proto), std::move(receiver_proto))); });

        set<size_t> expected_intersec;
        expected_l1_intersect<TR, TS, D, DELTA>(aes, 
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

        const double nMBsExchanged = ((double)(socks[0].bytesSent()+socks[0].bytesReceived()))/1024.0/1024.0; 

        SUCCEED("Number of MBs exchanged: " << nMBsExchanged);
    };

}

TEST_CASE("spl1 (t_s=65536 t_r=67108864 d=10 delta=30 ssp=40)","[spl1][n=m=2^16]") {

    BENCHMARK_ADVANCED("t_s=65536, t_r=67108864, d=10, delta=30, ssp=40")(Catch::Benchmark::Chronometer meter) {
        constexpr size_t TS = 65536;
        constexpr size_t TR = 67108864;
        constexpr size_t D = 10;
        constexpr size_t DELTA = 30;
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

        sparse_comp::sp_l1::Sender<TR,TS,D,DELTA,ssp> spL1Sender(senderPRNG, aes);
        sparse_comp::sp_l1::Receiver<TS,TR,D,DELTA,ssp> spL1Recvr(receiverPRNG, aes);

        auto sender_proto = spL1Sender.send(socks[0], *senderSparsePoints, *sender_in_values);
        auto receiver_proto = spL1Recvr.receive(socks[1], *receiverSparsePoints, *receiver_in_values, intersec);

        meter.measure([&sender_proto,&receiver_proto]() { sync_wait(when_all_ready(std::move(sender_proto), std::move(receiver_proto))); });

        set<size_t> expected_intersec;
        expected_l1_intersect<TR, TS, D, DELTA>(aes, 
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

        const double nMBsExchanged = ((double)(socks[0].bytesSent()+socks[0].bytesReceived()))/1024.0/1024.0; 

        SUCCEED("Number of MBs exchanged: " << nMBsExchanged);
    };

}


// END OF TESTS FOR N=M=2^16