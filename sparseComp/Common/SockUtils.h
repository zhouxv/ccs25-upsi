#pragma once

#include "coproto/Socket/Socket.h"
#include <vector>
#include <cmath>
#include <span>

//static int64_t MAX_SEND_SIZE_BYTES(2147483648L); // 2 GBs

namespace sparse_comp {

    constexpr int64_t COPROTO_MAX_SEND_SIZE_BYTES(2147483648L); // 2 GBs

    template <typename T, int64_t max_send_size_bytes>
    coproto::task<void> send(coproto::Socket& sock, std::vector<T>& v) {
        MC_BEGIN(coproto::task<void>, &sock, &v,
            total_send_size = size_t(0),
            i = size_t(0),
            max_item_per_round = size_t(0),
            n_items_togo = size_t(0),
            n_items_sent = size_t(0),
            v_span = std::span<T>(v),
            send_size = size_t(0),
            send_span = std::span<T>{});

    
            max_item_per_round = max_send_size_bytes / sizeof(T);

            n_items_togo = v.size();
            n_items_sent = 0;

            /*std::cout << "(s) max_item_per_round: " << max_item_per_round << std::endl;
            std::cout << "(s) n_items_togo: " << n_items_togo << std::endl;
            std::cout << "(s) sizeof(T): " << sizeof(T) << std::endl;*/

            while (n_items_togo > 0) {
                send_size = std::min(max_item_per_round, n_items_togo);
                send_span = v_span.subspan(n_items_sent, send_size);

                MC_AWAIT(sock.send(send_span));

                n_items_togo -= send_size;
                n_items_sent += send_size;
            }

        MC_END();
    }

    template <typename T, int64_t max_send_size_bytes>
    coproto::task<void> receive(coproto::Socket& sock, size_t recv_vec_size, std::vector<T>& v) {
        MC_BEGIN(coproto::task<void>, &sock, &v, recv_vec_size,
    v_span = std::span<T>{},
    recv_span = std::span<T>{},
    num_rounds = size_t(0),
    max_items_per_round = size_t(0),
    n_items_togo = size_t(0),
    n_items_sent = size_t(0),
    n_items_recv = size_t(0),
    i = size_t(0));
    
        v.resize(recv_vec_size);
        v_span = std::span<T>(v);

        max_items_per_round = max_send_size_bytes / sizeof(T);
        num_rounds = std::ceil(double(recv_vec_size) / double(max_items_per_round));

        n_items_togo = v.size();
        n_items_sent = 0;

        /*std::cout << "(r) recv_vec_size: " << recv_vec_size << std::endl;
        std::cout << "(r) max_items_per_round: " << max_items_per_round << std::endl;
        std::cout << "(r) n_items_togo: " << n_items_togo << std::endl;
        std::cout << "(r) sizeof(T): " << sizeof(T) << std::endl;
        std::cout << "(r) num_rounds: " << num_rounds << std::endl;*/

        for (i = 0; i < num_rounds; i++) {
            n_items_recv = std::min(max_items_per_round, n_items_togo);
            recv_span = v_span.subspan(n_items_sent, n_items_recv);

            //std::cout << "before recv" << std::endl;
            MC_AWAIT(sock.recv(recv_span));
            //std::cout << "after recv" << std::endl;
            

            n_items_togo -= n_items_recv;
            n_items_sent += n_items_recv;
        }

    MC_END();
    }

}