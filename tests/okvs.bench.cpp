#include "catch2/catch_test_macros.hpp"
#include "catch2/benchmark/catch_benchmark.hpp"
#include "cryptoTools/Common/block.h"
#include "volePSI/Paxos.h"
#include "../sparseComp/Common/BaxosUtils.h"

using osuCrypto::block;
using osuCrypto::PRNG;
using Baxos = volePSI::Baxos;
using PaxosParam = volePSI::PaxosParam;
using std::vector;

static void gen_keys_vals(vector<block>& keys, vector<block>& vals, size_t n) {
    PRNG prng = PRNG(block(9536629726117351353ULL,2724349864741298565ULL));

    keys.resize(n);
    vals.resize(n);

    for (size_t i = 0; i < n; i++) {
        keys[i] = prng.get<block>();
        vals[i] = prng.get<block>();
    }
}

static void gen_decode_keys(vector<block>& keys, size_t n) {
    PRNG prng = PRNG(block(9537729726117351353ULL,2724319864747298360ULL));

    keys.resize(n);

    for (size_t i = 0; i < n; i++) {
        keys[i] = prng.get<block>();
    }
}


static size_t baxos_enc_dec(vector<block>& keys, vector<block>& vals, vector<block>& decoded_keys, size_t ssp) {
    std::vector<block>* okvs = new std::vector<block>();
    size_t n_encd_items = keys.size();

    Baxos senderPaxos;
    senderPaxos.init(n_encd_items, sparse_comp::baxosBinSize(n_encd_items), 3, ssp, PaxosParam::GF128, oc::ZeroBlock);
    
    okvs->resize(senderPaxos.size());

    senderPaxos.solve<block>(keys, vals, *okvs, nullptr, 1);

    Baxos receiverPaxos;
    receiverPaxos.init(n_encd_items, sparse_comp::baxosBinSize(n_encd_items), 3, ssp, PaxosParam::GF128, oc::ZeroBlock);
    
    vector<block>* decoded_vals = new vector<block>(decoded_keys.size());

    receiverPaxos.decode<block>(decoded_keys, *decoded_vals, *okvs);

    size_t okvs_size = okvs->size();

    delete okvs;
    delete decoded_vals;

    return okvs_size;
}

static size_t bench_baxos_enc_dec(Catch::Benchmark::Chronometer meter, size_t n_encd_items, size_t n_decoded_items, size_t ssp) {
    auto keys = new vector<block>(n_encd_items), vals = new vector<block>(n_encd_items), decoded_keys = new vector<block>(n_decoded_items);
    
    gen_keys_vals(*keys, *vals, n_encd_items);
    gen_decode_keys(*decoded_keys, n_decoded_items);

    size_t paxos_size;

    meter.measure([keys, vals, decoded_keys, &paxos_size, ssp] {
        paxos_size = baxos_enc_dec(*keys, *vals, *decoded_keys, ssp);
    });

    delete keys;
    delete vals;
    delete decoded_keys;

    return paxos_size;
}

TEST_CASE("baxos enc/dec (nA=nB=2^8, d=2, ssp=40)", "[baxos][nA=nB=2^8][d=2]") {
    BENCHMARK_ADVANCED("nA=nB=2^8, d=2, ssp = 40")(Catch::Benchmark::Chronometer meter) {
        size_t n_encd_items = 256;
        size_t n_decoded_items = 256*4;
        
        size_t paxos_size = bench_baxos_enc_dec(meter, n_encd_items, n_decoded_items, 40);

        double nKBsPaxos = ((double) paxos_size) * 16.0 / 1024.0;

        SUCCEED("Paxos size (KBs): " << nKBsPaxos);
    };
}

TEST_CASE("baxos enc/dec (nA=nB=2^8, d=6, ssp=40)", "[baxos][nA=nB=2^8][d=6]") {
    BENCHMARK_ADVANCED("nA=nB=2^8, d=6, ssp=40")(Catch::Benchmark::Chronometer meter) {
        size_t n_encd_items = 256;
        size_t n_decoded_items = 256*64;
        
        size_t paxos_size = bench_baxos_enc_dec(meter, n_encd_items, n_decoded_items, 40);

        double nKBsPaxos = ((double) paxos_size) * 16.0 / 1024.0;

        SUCCEED("Paxos size (KBs): " << nKBsPaxos);
    };
}

TEST_CASE("baxos enc/dec (nA=nB=2^8, d=10, ssp=40)", "[baxos][nA=nB=2^8][d=10]") {

    BENCHMARK_ADVANCED("nA=nB=2^8, d=10, ssp=40")(Catch::Benchmark::Chronometer meter) {
        size_t n_encd_items = 256;
        size_t n_decoded_items = 256*1024;
        
        size_t paxos_size = bench_baxos_enc_dec(meter, n_encd_items, n_decoded_items, 40);

        double nKBsPaxos = ((double) paxos_size) * 16.0 / 1024.0;

        SUCCEED("Paxos size (KBs): " << nKBsPaxos);
    };
}

TEST_CASE("baxos enc/dec (nA=nB=2^12, d=2, ssp=40)", "[baxos][nA=nB=2^12][d=2]") {
    BENCHMARK_ADVANCED("nA=nB=2^12, d=2, ssp = 40")(Catch::Benchmark::Chronometer meter) {
        size_t n_encd_items = 4096;
        size_t n_decoded_items = 4096*4;
        
        size_t paxos_size = bench_baxos_enc_dec(meter, n_encd_items, n_decoded_items, 40);

        double nKBsPaxos = ((double) paxos_size) * 16.0 / 1024.0;

        SUCCEED("Paxos size (KBs): " << nKBsPaxos);
    };
}

TEST_CASE("baxos enc/dec (nA=nB=2^12, d=6, ssp=40)", "[baxos][nA=nB=2^12][d=6]") {
    BENCHMARK_ADVANCED("nA=nB=2^12, d=6, ssp=40")(Catch::Benchmark::Chronometer meter) {
        size_t n_encd_items = 4096;
        size_t n_decoded_items = 4096*64;
        
        size_t paxos_size = bench_baxos_enc_dec(meter, n_encd_items, n_decoded_items, 40);

        double nKBsPaxos = ((double) paxos_size) * 16.0 / 1024.0;

        SUCCEED("Paxos size (KBs): " << nKBsPaxos);
    };
}

TEST_CASE("baxos enc/dec (nA=nB=2^12, d=10, ssp=40)", "[baxos][nA=nB=2^12][d=10]") {

    BENCHMARK_ADVANCED("nA=nB=2^12, d=10, ssp=40")(Catch::Benchmark::Chronometer meter) {
        size_t n_encd_items = 4096;
        size_t n_decoded_items = 4096*1024;
        
        size_t paxos_size = bench_baxos_enc_dec(meter, n_encd_items, n_decoded_items, 40);

        double nKBsPaxos = ((double) paxos_size) * 16.0 / 1024.0;

        SUCCEED("Paxos size (KBs): " << nKBsPaxos);
    };
}

TEST_CASE("baxos enc/dec (nA=nB=2^16, d=2, ssp=40)", "[baxos][nA=nB=2^16][d=2]") {
    BENCHMARK_ADVANCED("nA=nB=2^16, d=2, ssp = 40")(Catch::Benchmark::Chronometer meter) {
        size_t n_encd_items = 65536;
        size_t n_decoded_items = 65536*4;
        
        size_t paxos_size = bench_baxos_enc_dec(meter, n_encd_items, n_decoded_items, 40);

        double nKBsPaxos = ((double) paxos_size) * 16.0 / 1024.0;

        SUCCEED("Paxos size (KBs): " << nKBsPaxos);
    };
}

TEST_CASE("baxos enc/dec (nA=nB=2^16, d=6, ssp=40)", "[baxos][nA=nB=2^16][d=6]") {
    BENCHMARK_ADVANCED("nA=nB=2^16, d=6, ssp=40")(Catch::Benchmark::Chronometer meter) {
        size_t n_encd_items = 65536;
        size_t n_decoded_items = 65536*64;
        
        size_t paxos_size = bench_baxos_enc_dec(meter, n_encd_items, n_decoded_items, 40);

        double nKBsPaxos = ((double) paxos_size) * 16.0 / 1024.0;

        SUCCEED("Paxos size (KBs): " << nKBsPaxos);
    };
}

TEST_CASE("baxos enc/dec (nA=nB=2^16, d=10, ssp=40)", "[baxos][nA=nB=2^16][d=10]") {

    BENCHMARK_ADVANCED("nA=nB=2^16, d=10, ssp=40")(Catch::Benchmark::Chronometer meter) {
        size_t n_encd_items = 65536;
        size_t n_decoded_items = 65536*1024;
        
        size_t paxos_size = bench_baxos_enc_dec(meter, n_encd_items, n_decoded_items, 40);

        double nKBsPaxos = ((double) paxos_size) * 16.0 / 1024.0;

        SUCCEED("Paxos size (KBs): " << nKBsPaxos);
    };
}

