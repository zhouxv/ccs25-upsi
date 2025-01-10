#include "catch2/catch_test_macros.hpp"
#include "catch2/benchmark/catch_benchmark.hpp"
#include "cryptoTools/Crypto/PRNG.h"
#include "cryptoTools/Common/block.h"
#include "cryptoTools/Crypto/AES.h"

using PRNG = osuCrypto::PRNG;
using osuCrypto::block;
using AES = osuCrypto::AES;
using AESDec = osuCrypto::AESDec;

using std::vector;

static void rand_block_fill(vector<block>& vec) {
    PRNG prng = PRNG(block(9536629726117351353ULL,2724349864741298565ULL));

    for (size_t i = 0; i < vec.size(); i++) {
        vec[i] = prng.get<block>();
    }

}

static void enc(vector<block>& ks, vector<block>& pts) {
    auto cs = new vector<block>(pts.size());

    AES aes;

    for (size_t i = 0; i < pts.size(); i++) {
        aes.setKey(ks[i]);
        aes.ecbEncBlock(pts[i],cs->at(i));
    }

    delete cs;
}

static void dec(vector<block>& ks, vector<block>& cs) {
    auto pts = new vector<block>(cs.size());

    AESDec aes;

    for (size_t i = 0; i < cs.size(); i++) {
        aes.setKey(ks[i]);
        aes.ecbDecBlock(cs[i],pts->at(i));
    }

    delete pts;
}

static void bench_enc_dec(Catch::Benchmark::Chronometer meter, size_t n_enc, size_t n_dec) {
    auto enc_ks = new vector<block>(n_enc), dec_ks = new vector<block>(n_dec), pts = new vector<block>(n_enc), cs = new vector<block>(n_dec);

    rand_block_fill(*enc_ks);
    rand_block_fill(*dec_ks);
    rand_block_fill(*pts);
    rand_block_fill(*cs);

    meter.measure([&] {
        enc(*enc_ks,*pts);
        dec(*dec_ks,*cs);
    });

}

TEST_CASE("enc and dec (nA=nB=2^8, d=2)", "[enc_dec][nA=nB=2^8][d=2]") {
    BENCHMARK_ADVANCED("nA=nB=2^8, d=2")(Catch::Benchmark::Chronometer meter) {
        size_t n_enc = 256;
        size_t n_dec = 256*4;
        
        bench_enc_dec(meter, n_enc, n_dec);
    };
}

TEST_CASE("enc and dec (nA=nB=2^8, d=6)", "[enc_dec][nA=nB=2^8][d=6]") {
    BENCHMARK_ADVANCED("nA=nB=2^8, d=6")(Catch::Benchmark::Chronometer meter) {
        size_t n_enc = 256;
        size_t n_dec = 256*64;
        
        bench_enc_dec(meter, n_enc, n_dec);
    };
}

TEST_CASE("enc and dec (nA=nB=2^8, d=10)", "[enc_dec][nA=nB=2^8][d=10]") {
    BENCHMARK_ADVANCED("nA=nB=2^8, d=10")(Catch::Benchmark::Chronometer meter) {
        size_t n_enc = 256;
        size_t n_dec = 256*1024;
        
        bench_enc_dec(meter, n_enc, n_dec);
    };
}

TEST_CASE("enc and dec (nA=nB=2^12, d=2)", "[enc_dec][nA=nB=2^12][d=2]") {
    BENCHMARK_ADVANCED("nA=nB=2^12, d=2")(Catch::Benchmark::Chronometer meter) {
        size_t n_enc = 4096;
        size_t n_dec = 4096*4;
        
        bench_enc_dec(meter, n_enc, n_dec);
    };
}

TEST_CASE("enc and dec (nA=nB=2^12, d=6)", "[enc_dec][nA=nB=2^12][d=6]") {
    BENCHMARK_ADVANCED("nA=nB=2^12, d=6")(Catch::Benchmark::Chronometer meter) {
        size_t n_enc = 4096;
        size_t n_dec = 4096*64;
        
        bench_enc_dec(meter, n_enc, n_dec);
    };
}

TEST_CASE("enc and dec (nA=nB=2^12, d=10)", "[enc_dec][nA=nB=2^12][d=10]") {
    BENCHMARK_ADVANCED("nA=nB=2^12, d=10")(Catch::Benchmark::Chronometer meter) {
        size_t n_enc = 4096;
        size_t n_dec = 4096*1024;
        
        bench_enc_dec(meter, n_enc, n_dec);
    };
}

TEST_CASE("enc and dec (nA=nB=2^16, d=2)", "[enc_dec][nA=nB=2^16][d=2]") {
    BENCHMARK_ADVANCED("nA=nB=2^16, d=2")(Catch::Benchmark::Chronometer meter) {
        size_t n_enc = 65536;
        size_t n_dec = 65536*4;
        
        bench_enc_dec(meter, n_enc, n_dec);
    };
}

TEST_CASE("enc and dec (nA=nB=2^16, d=6)", "[enc_dec][nA=nB=2^16][d=6]") {
    BENCHMARK_ADVANCED("nA=nB=2^16, d=6")(Catch::Benchmark::Chronometer meter) {
        size_t n_enc = 65536;
        size_t n_dec = 65536*64;
        
        bench_enc_dec(meter, n_enc, n_dec);
    };
}

TEST_CASE("enc and dec (nA=nB=2^16, d=10)", "[enc_dec][nA=nB=2^16][d=10]") {
    BENCHMARK_ADVANCED("nA=nB=2^16, d=10")(Catch::Benchmark::Chronometer meter) {
        size_t n_enc = 65536;
        size_t n_dec = 65536*1024;
        
        bench_enc_dec(meter, n_enc, n_dec);
    };
}