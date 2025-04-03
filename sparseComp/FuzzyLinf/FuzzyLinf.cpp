#include "./FuzzyLinf.h"
#include "../SpLInf/SpLInf.h"
#include "../Common/HashUtils.h"
#include "../Common/Common.h"
#include <array>
#include <cstdint>
#include <iostream>

template<size_t tr, size_t t, size_t d, uint8_t delta, uint8_t ssp>
using SpLinfSender = sparse_comp::sp_linf::Sender<tr,t,d,delta,ssp>;

template<size_t ts, size_t t, size_t d, uint8_t delta, uint8_t ssp>
using SpLinfReceiver = sparse_comp::sp_linf::Receiver<ts,t,d,delta,ssp>;

template<size_t t, size_t d>
static void points_to_in_values(std::array<point,t>& points, std::array<std::array<uint32_t,d>,t>& in_values) {

    for (size_t i = 0; i < t; i++) {
        for (size_t j = 0; j < d; j++) {
            in_values[i][j] = points[i].coords[j];
        }
    }

}

template<size_t tr, size_t t, size_t d, uint8_t delta, uint8_t ssp>
Proto sparse_comp::fuzzy_linf::Sender<tr,t,d,delta,ssp>::send(
                                                     Socket& sock, 
                                                     array<point,t>& points) {
    MC_BEGIN(Proto, this, &sock, &points, 
             spLinfSender = (SpLinfSender<tr,t,d,delta,ssp>*) nullptr,
             cells = (array<point,t>*) nullptr,
             in_values = (array<array<uint32_t,d>,t>*) nullptr,
             prt = Proto());

        spLinfSender = new SpLinfSender<tr,t,d,delta,ssp>(*(this->prng), *(this->aes));
        cells = new array<point,t>();
        in_values = new array<array<uint32_t,d>,t>();

        // Maps points to cells using spatial hashing
        sparse_comp::spatial_hash<t>(points, *cells, d, delta);

        // Maps points to in_values
        points_to_in_values<t,d>(points, *in_values);

        prt = spLinfSender->send(sock, *cells, *in_values);

        MC_AWAIT(prt);

        delete spLinfSender;
        delete cells;
        delete in_values;
    
    MC_END();
}

template<size_t ts, size_t t, size_t d, uint8_t delta, uint8_t ssp>
Proto sparse_comp::fuzzy_linf::Receiver<ts,t,d,delta,ssp>::receive(
                                                     Socket& sock, 
                                                     array<point,t>& points,
                                                     vector<size_t>& intersec_pos) {
    MC_BEGIN(Proto, this, &sock, &points, &intersec_pos,
             spLinfReceiver = (SpLinfReceiver<ts,t,d,delta,ssp>*) nullptr,
             cells = (array<point,t>*) nullptr,
             in_values = (array<array<uint32_t,d>,t>*) nullptr,
             prt = Proto());

        spLinfReceiver = new SpLinfReceiver<ts,t,d,delta,ssp>(*(this->prng), *(this->aes));
        cells = new array<point,t>();
        in_values = new array<array<uint32_t,d>,t>();

        // Maps points to cells using spatial hashing
        sparse_comp::spatial_hash<t>(points, *cells, d, delta);

        // Maps points to in_values
        points_to_in_values<t,d>(points, *in_values);

        prt = spLinfReceiver->receive(sock, *cells, *in_values, intersec_pos);

        MC_AWAIT(prt);

        delete spLinfReceiver;
        delete cells;
        delete in_values;

    MC_END();
}