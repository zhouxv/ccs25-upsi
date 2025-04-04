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
static void sndr_points_to_in_values(std::array<point,t>& points, std::array<std::array<uint32_t,d>,t>& in_values) {

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
    constexpr const size_t twotod = (size_t) pow(2, d);
    constexpr const size_t rcvr_cell_count = twotod * tr;
    
    MC_BEGIN(Proto, this, &sock, &points, 
             spLinfSender = (SpLinfSender<rcvr_cell_count,t,d,delta,ssp>*) nullptr,
             point_hashs = (array<point,t>*) nullptr,
             in_values = (array<array<uint32_t,d>,t>*) nullptr,
             prt = Proto());

        spLinfSender = new SpLinfSender<rcvr_cell_count, t, d, delta, ssp>(*(this->prng), *(this->aes));
        point_hashs = new array<point, t>();
        in_values = new array<array<uint32_t,d>,t>();

        // Maps points to cells using spatial hashing
        sparse_comp::spatial_hash<t>(points, *point_hashs, d, delta);

        // Maps points to in_values
        sndr_points_to_in_values<t,d>(points, *in_values);

        prt = spLinfSender->send(sock, *point_hashs, *in_values);

        MC_AWAIT(prt);

        delete spLinfSender;
        delete point_hashs;
        delete in_values;
    
    MC_END();
}



template<size_t t, size_t d, size_t cell_count>
static void sndr_points_to_in_values(std::array<point,t>& points, std::array<std::array<uint32_t,d>,t>& in_values) {

    for (size_t i = 0; i < t; i++) {
        for (size_t j = 0; j < d; j++) {
            in_values[i][j] = points[i].coords[j];
        }
    }

}

template<size_t t, size_t d, size_t cell_count>
static void rcvr_points_to_in_values(std::array<point, t>& point_center, std::array<std::array<uint32_t,d>, cell_count>& in_values) {
        constexpr const size_t twotod = (size_t) pow(2, d);

        static_assert(cell_count == twotod*t);
        static_assert(d <= point::MAX_DIM);

        for (size_t i = 0; i < t; i++) {
            for (size_t j = 0; j < twotod; j++) {
                for (size_t k = 0; k < d; k++) {
                    in_values[twotod*i+j][k] = point_center[i].coords[k];
                }
            }
        }

}

template<size_t ts, size_t t, size_t d, uint8_t delta, uint8_t ssp>
Proto sparse_comp::fuzzy_linf::Receiver<ts,t,d,delta,ssp>::receive(
                                                     Socket& sock, 
                                                     array<point,t>& points,
                                                     vector<size_t>& intersec_pos) {
    constexpr const size_t twotod = (size_t) pow(2, d);
    constexpr const size_t cell_count = twotod * t;
    
    MC_BEGIN(Proto, this, &sock, &points, &intersec_pos,
             spLinfReceiver = (SpLinfReceiver<ts,cell_count,d,delta,ssp>*) nullptr,
             cells = (array<point, cell_count>*) nullptr,
             in_values = (array<array<uint32_t,d>,cell_count>*) nullptr,
             prt = Proto());

        spLinfReceiver = new SpLinfReceiver<ts,cell_count,d,delta,ssp>(*(this->prng), *(this->aes));
        cells = new array<point,cell_count>();
        in_values = new array<array<uint32_t,d>,cell_count>();

        // Maps points to adjcent cells using spatial hashing
        sparse_comp::spatial_cell_hash<t,d,cell_count>(points, *cells, delta);

        // Maps points to in_values
        rcvr_points_to_in_values<t,d,cell_count>(points, *in_values);

        prt = spLinfReceiver->receive(sock, *cells, *in_values, intersec_pos);

        MC_AWAIT(prt);

        for (size_t i = 0; i < intersec_pos.size(); i++) {
            intersec_pos[i] = intersec_pos[i] / twotod;
        }

        delete spLinfReceiver;
        delete cells;
        delete in_values;

    MC_END();
}