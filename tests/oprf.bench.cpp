#include "catch2/catch_test_macros.hpp"
#include "catch2/benchmark/catch_benchmark.hpp"
#include "../sparseComp/MultiOPRF/MultiOPRF.h"
#include "coproto/Socket/Socket.h"
#include "cryptoTools/Common/block.h"
#include <vector>
#include <cmath>
#include <iostream>

using std::pow;
using std::vector;
using coproto::LocalAsyncSocket;

using macoro::sync_wait;
using macoro::when_all_ready;

using namespace sparse_comp::multi_oprf;

double bench_oprf(Catch::Benchmark::Chronometer meter, 
                size_t n_oprf_instances, 
                vector<size_t> n_qs, 
                vector<size_t> n_es) {

    assert(n_qs.size() == n_oprf_instances && n_es.size() == n_oprf_instances);

    auto socks = LocalAsyncSocket::makePair();
    PRNG senderPRNG = PRNG(block(742130310438916676ULL, 11803924226990735076ULL));
    PRNG receiverPRNG = PRNG(block(2457938039974938056ULL, 17910068785450354990ULL));

    vector<vector<block>*> qpts(n_oprf_instances);
    vector<vector<block>*> qvals(n_oprf_instances);
    vector<vector<block>*> epts(n_oprf_instances);
    vector<vector<block>*> evals(n_oprf_instances);

    for (size_t i = 0; i < n_oprf_instances; i++) {
        qpts[i] = new vector<block>(n_qs[i]);
        qvals[i] = new vector<block>(n_qs[i]);
        epts[i] = new vector<block>(n_es[i]);
        evals[i] = new vector<block>(n_es[i]);

        for (size_t j = 0; j < n_qs[i]; j++) {
            qpts[i]->at(j) = senderPRNG.get<block>();
        }

        for (size_t j = 0; j < n_es[i]; j++) {
            epts[i]->at(j) = senderPRNG.get<block>();
        }
    }

    vector<Sender*> senders(n_oprf_instances);
    vector<Receiver*> receivers(n_oprf_instances);

    meter.measure([n_oprf_instances, &n_es, &n_qs, &socks, &senders, &receivers, &senderPRNG, &receiverPRNG, &qpts, &qvals, &epts, &evals] {
        Proto sender_setup = Sender::setup(socks[0], senderPRNG, n_oprf_instances, senders);
        Proto receiver_setup = Receiver::setup(socks[1], receiverPRNG, n_oprf_instances, receivers);

        sync_wait(when_all_ready(sender_setup, receiver_setup));

        for (size_t i = 0; i < n_oprf_instances; i++) {
            auto p1 = senders[i]->send(socks[0], n_qs[i]);
            auto p2 = receivers[i]->receive(socks[1], *qpts[i], *qvals[i]);
        
            sync_wait(when_all_ready(p1, p2));

            senders[i]->eval(*epts[i], *evals[i]);
        }
    });

    const double nMBsExchanged = ((double)(socks[0].bytesSent()+socks[0].bytesReceived()))/1024.0/1024.0; 

    for (size_t i = 0; i < n_oprf_instances; i++) {
        delete qpts[i];
        delete qvals[i];
        delete epts[i];
        delete evals[i];
        delete senders[i];
        delete receivers[i];
    }

    return nMBsExchanged;
}

TEST_CASE("oprf (n=1, q=1, e=1)", "[oprf][n=1][q=1][e=1]") {
    BENCHMARK_ADVANCED("n=1, q=1, e=1")(Catch::Benchmark::Chronometer meter) {
        size_t n_oprf_instances = 1;
        vector<size_t> n_qs = {100};
        vector<size_t> n_es = {100};

        double nMBsExchanged = bench_oprf(meter, n_oprf_instances, n_qs, n_es);

        SUCCEED("Number of MBs exchanged: " << nMBsExchanged);
    };
}

TEST_CASE("nA=nB=2^8, d=2", "[oprf][nA=nB=2^8][d=2]") {
    double nMBsExchanged = -1;
    
    BENCHMARK_ADVANCED("nA=nB=2^8, d=2")(Catch::Benchmark::Chronometer meter) {
        size_t set_size = 256;
        size_t d = 2;

        size_t two_to_d = pow(2,d);
        size_t n_oprf_instances = 2;
        vector<size_t> n_qs = {set_size*d, set_size};
        vector<size_t> n_es = {set_size*two_to_d*d, set_size*two_to_d};

        nMBsExchanged = bench_oprf(meter, n_oprf_instances, n_qs, n_es);

    };

    std::cout << "Number of MBs exchanged: " << nMBsExchanged << std::endl;
}

TEST_CASE("nA=nB=2^8, d=6", "[oprf][nA=nB=2^8][d=6]") {
    double nMBsExchanged = -1;

    BENCHMARK_ADVANCED("nA=nB=2^8, d=6")(Catch::Benchmark::Chronometer meter) {
        size_t set_size = 256;
        size_t d = 6;
        
        size_t two_to_d = pow(2,d);
        size_t n_oprf_instances = 2;
        vector<size_t> n_qs = {set_size*d, set_size};
        vector<size_t> n_es = {set_size*two_to_d*d, set_size*two_to_d};

        nMBsExchanged = bench_oprf(meter, n_oprf_instances, n_qs, n_es);

        // SUCCEED("Number of MBs exchanged: " << nMBsExchanged);
    };

    std::cout << "Number of MBs exchanged: " << nMBsExchanged << std::endl;
}

TEST_CASE("nA=nB=2^8, d=10", "[oprf][nA=nB=2^8][d=10]") {
    double nMBsExchanged = -1;

    BENCHMARK_ADVANCED("nA=nB=2^8, d=10")(Catch::Benchmark::Chronometer meter) {
        size_t set_size = 256;
        size_t d = 10;
        
        size_t two_to_d = pow(2,d);
        size_t n_oprf_instances = 2;
        vector<size_t> n_qs = {set_size*d, set_size};
        vector<size_t> n_es = {set_size*two_to_d*d, set_size*two_to_d};

        nMBsExchanged = bench_oprf(meter, n_oprf_instances, n_qs, n_es);

        //SUCCEED("Number of MBs exchanged: " << nMBsExchanged);
    };

    std::cout << "Number of MBs exchanged: " << nMBsExchanged << std::endl;
}

TEST_CASE("nA=nB=2^12, d=2", "[oprf][nA=nB=2^12][d=2]") {
    double nMBsExchanged = -1;
    
    BENCHMARK_ADVANCED("nA=nB=2^12, d=2")(Catch::Benchmark::Chronometer meter) {
        size_t set_size = 4096;
        size_t d = 2;

        size_t two_to_d = pow(2,d);
        size_t n_oprf_instances = 2;
        vector<size_t> n_qs = {set_size*d, set_size};
        vector<size_t> n_es = {set_size*two_to_d*d, set_size*two_to_d};

        nMBsExchanged = bench_oprf(meter, n_oprf_instances, n_qs, n_es);

    };

    std::cout << "Number of MBs exchanged: " << nMBsExchanged << std::endl;
}

TEST_CASE("nA=nB=2^12, d=6", "[oprf][nA=nB=2^12][d=6]") {
    double nMBsExchanged = -1;

    BENCHMARK_ADVANCED("nA=nB=2^12, d=6")(Catch::Benchmark::Chronometer meter) {
        size_t set_size = 4096;
        size_t d = 6;
        
        size_t two_to_d = pow(2,d);
        size_t n_oprf_instances = 2;
        vector<size_t> n_qs = {set_size*d, set_size};
        vector<size_t> n_es = {set_size*two_to_d*d, set_size*two_to_d};

        nMBsExchanged = bench_oprf(meter, n_oprf_instances, n_qs, n_es);

        // SUCCEED("Number of MBs exchanged: " << nMBsExchanged);
    };

    std::cout << "Number of MBs exchanged: " << nMBsExchanged << std::endl;
}

TEST_CASE("nA=nB=2^12, d=10", "[oprf][nA=nB=2^12][d=10]") {
    double nMBsExchanged = -1;

    BENCHMARK_ADVANCED("nA=nB=2^12, d=10")(Catch::Benchmark::Chronometer meter) {
        size_t set_size = 4096;
        size_t d = 10;
        
        size_t two_to_d = pow(2,d);
        size_t n_oprf_instances = 2;
        vector<size_t> n_qs = {set_size*d, set_size};
        vector<size_t> n_es = {set_size*two_to_d*d, set_size*two_to_d};

        nMBsExchanged = bench_oprf(meter, n_oprf_instances, n_qs, n_es);

        //SUCCEED("Number of MBs exchanged: " << nMBsExchanged);
    };

    std::cout << "Number of MBs exchanged: " << nMBsExchanged << std::endl;
}

TEST_CASE("nA=nB=2^16, d=2", "[oprf][nA=nB=2^16][d=2]") {
    double nMBsExchanged = -1;
    
    BENCHMARK_ADVANCED("nA=nB=2^16, d=2")(Catch::Benchmark::Chronometer meter) {
        size_t set_size = 65536;
        size_t d = 2;

        size_t two_to_d = pow(2,d);
        size_t n_oprf_instances = 2;
        vector<size_t> n_qs = {set_size*d, set_size};
        vector<size_t> n_es = {set_size*two_to_d*d, set_size*two_to_d};

        nMBsExchanged = bench_oprf(meter, n_oprf_instances, n_qs, n_es);

    };

    std::cout << "Number of MBs exchanged: " << nMBsExchanged << std::endl;
}

TEST_CASE("nA=nB=2^16, d=6", "[oprf][nA=nB=2^16][d=6]") {
    double nMBsExchanged = -1;

    BENCHMARK_ADVANCED("nA=nB=2^16, d=6")(Catch::Benchmark::Chronometer meter) {
        size_t set_size = 65536;
        size_t d = 6;
        
        size_t two_to_d = pow(2,d);
        size_t n_oprf_instances = 2;
        vector<size_t> n_qs = {set_size*d, set_size};
        vector<size_t> n_es = {set_size*two_to_d*d, set_size*two_to_d};

        nMBsExchanged = bench_oprf(meter, n_oprf_instances, n_qs, n_es);

        // SUCCEED("Number of MBs exchanged: " << nMBsExchanged);
    };

    std::cout << "Number of MBs exchanged: " << nMBsExchanged << std::endl;
}

TEST_CASE("nA=nB=2^16, d=10", "[oprf][nA=nB=2^16][d=10]") {
    double nMBsExchanged = -1;

    BENCHMARK_ADVANCED("nA=nB=2^16, d=10")(Catch::Benchmark::Chronometer meter) {
        size_t set_size = 65536;
        size_t d = 10;
        
        size_t two_to_d = pow(2,d);
        size_t n_oprf_instances = 2;
        vector<size_t> n_qs = {set_size*d, set_size};
        vector<size_t> n_es = {set_size*two_to_d*d, set_size*two_to_d};

        nMBsExchanged = bench_oprf(meter, n_oprf_instances, n_qs, n_es);

        //SUCCEED("Number of MBs exchanged: " << nMBsExchanged);
    };

    std::cout << "Number of MBs exchanged: " << nMBsExchanged << std::endl;
}
