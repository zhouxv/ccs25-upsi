#include "catch2/catch_test_macros.hpp"
#include "../sparseComp/Common/SockUtils.h"
#include "cryptoTools/Crypto/PRNG.h"
#include "cryptoTools/Common/block.h"
#include "coproto/Socket/Socket.h"
#include "coproto/Socket/LocalAsyncSock.h"
#include "macoro/task.h"
#include <utility>

using osuCrypto::block;
using coproto::LocalAsyncSocket;

TEST_CASE("sequential uint64_t transmission") {
    constexpr int64_t max_send_size_bytes = 17;

    const size_t seq_len = 1000;

    auto socks = LocalAsyncSocket::makePair();

    std::vector<uint64_t> v_recv(seq_len);
    std::vector<uint64_t> v_send(seq_len);
    for (size_t i = 0; i < v_send.size(); i++) {
        v_send[i] = i;
    }

    auto send_task = sparse_comp::send<uint64_t,max_send_size_bytes>(socks[0], v_send);
    auto recv_task = sparse_comp::receive<uint64_t,max_send_size_bytes>(socks[1], v_recv.size(), v_recv);

    coproto::sync_wait(coproto::when_all_ready(std::move(send_task),std::move(recv_task)));

    for (size_t i = 0; i < v_recv.size(); i++) {
        REQUIRE(v_recv[i] == i);
    }
}

TEST_CASE("random block transmission with n=1000, max_send_size_bytes = 17") {
    constexpr int64_t max_send_size_bytes = 33;
    const size_t seq_len = 1000;

    auto socks = LocalAsyncSocket::makePair();
    auto prng1 = osuCrypto::PRNG(block(13133210048402866,17132091720387928));
    auto prng2 = osuCrypto::PRNG(block(13133210048402866,17132091720387928));

    std::vector<block> v_recv(seq_len);
    std::vector<block> v_send(seq_len);
    for (size_t i = 0; i < v_send.size(); i++) {
        v_send[i] = prng1.get<block>();
    }

    auto send_task = sparse_comp::send<block,max_send_size_bytes>(socks[0], v_send);
    auto recv_task = sparse_comp::receive<block,max_send_size_bytes>(socks[1], v_recv.size(), v_recv);

    coproto::sync_wait(coproto::when_all_ready(std::move(send_task),std::move(recv_task)));

    for (size_t i = 0; i < v_recv.size(); i++) {
        REQUIRE(v_recv[i] == prng2.get<block>());
    }
}

TEST_CASE("random block transmission with n = 1000, max_send_size_bytes = 33") {
    constexpr int64_t max_send_size_bytes = 33;
    const size_t seq_len = 1000;

    auto socks = LocalAsyncSocket::makePair();
    auto prng1 = osuCrypto::PRNG(block(13133210048402866,17132091720387928));
    auto prng2 = osuCrypto::PRNG(block(13133210048402866,17132091720387928));

    std::vector<block> v_recv(seq_len);
    std::vector<block> v_send(seq_len);
    for (size_t i = 0; i < v_send.size(); i++) {
        v_send[i] = prng1.get<block>();
    }

    auto send_task = sparse_comp::send<block,max_send_size_bytes>(socks[0], v_send);
    auto recv_task = sparse_comp::receive<block,max_send_size_bytes>(socks[1], v_recv.size(), v_recv);

    coproto::sync_wait(coproto::when_all_ready(std::move(send_task),std::move(recv_task)));

    for (size_t i = 0; i < v_recv.size(); i++) {
        REQUIRE(v_recv[i] == prng2.get<block>());
    }
}

TEST_CASE("random block transmission with n = 1000, max_send_size_bytes = 13000") {
    constexpr int64_t max_send_size_bytes = 13000;
    const size_t seq_len = 1000;

    auto socks = LocalAsyncSocket::makePair();
    auto prng1 = osuCrypto::PRNG(block(13133210048402866,17132091720387928));
    auto prng2 = osuCrypto::PRNG(block(13133210048402866,17132091720387928));

    std::vector<block> v_recv(seq_len);
    std::vector<block> v_send(seq_len);
    for (size_t i = 0; i < v_send.size(); i++) {
        v_send[i] = prng1.get<block>();
    }

    auto send_task = sparse_comp::send<block,max_send_size_bytes>(socks[0], v_send);
    auto recv_task = sparse_comp::receive<block,max_send_size_bytes>(socks[1], v_recv.size(), v_recv);

    coproto::sync_wait(coproto::when_all_ready(std::move(send_task),std::move(recv_task)));

    for (size_t i = 0; i < v_recv.size(); i++) {
        REQUIRE(v_recv[i] == prng2.get<block>());
    }
}

TEST_CASE("random block transmission with n = 2797, max_send_size_bytes = 2147483648L") {

    constexpr int64_t max_send_size_bytes = 2147483648L;
    const size_t seq_len = 2797;

    auto socks = LocalAsyncSocket::makePair();
    auto prng1 = osuCrypto::PRNG(block(13133210048402866,17132091720387928));
    auto prng2 = osuCrypto::PRNG(block(13133210048402866,17132091720387928));

    std::vector<block> v_recv(3);
    std::vector<block> v_send(seq_len);
    for (size_t i = 0; i < v_send.size(); i++) {
        v_send[i] = prng1.get<block>();
    }

    auto send_task = sparse_comp::send<block,max_send_size_bytes>(socks[0], v_send);
    auto recv_task = sparse_comp::receive<block,max_send_size_bytes>(socks[1], seq_len, v_recv);

    coproto::sync_wait(coproto::when_all_ready(std::move(send_task),std::move(recv_task)));

    for (size_t i = 0; i < v_recv.size(); i++) {
        REQUIRE(v_recv[i] == prng2.get<block>());
    }

}