#include "./SpL2.h"
#include "../Common/ZN.h"
#include "../SpBSOT/SpBSOT.h"
#include "../BlockSpBSOT/BlockSpBSOT.h"
#include "../CustomOPRF/CustomizedOPRF.h"
#include <array>
#include <cmath>

#define MAX_SSP 128
#define IN_COMP_BIT_LEN 8

template<typename T, size_t N>
using array = std::array<T, N>;

template<size_t tr, size_t ts, size_t k, size_t n>
using BlockSpBSOTSender = sparse_comp::block_sp_bsot::Sender<tr,ts,k,n>; 

template<size_t ts, size_t tr, size_t k, size_t n>
using BlockSpBSOTReceiver = sparse_comp::block_sp_bsot::Receiver<ts,tr,k,n>; 

template<size_t tr, size_t t, size_t k, size_t n, uint64_t M>
using SpBSOTSender = sparse_comp::sp_bsot::Sender<tr,t,k,n,M>;

template<size_t ts, size_t t, size_t k, size_t n, uint64_t M>
using SpBSOTReceiver = sparse_comp::sp_bsot::Receiver<ts,t,k,n,M>;

using OprfSender = sparse_comp::custom_oprf::Sender;
using OprfReceiver = sparse_comp::custom_oprf::Receiver;

template<size_t t, size_t d, uint64_t N>
static void gen_zero_shares(array<array<ZN<N>,d>,t>& zero_shares) {

    for (size_t i=0;i < t;i++) {
        for (size_t j=0;j < d;j++) {
            zero_shares[i][j] = ZN<N>(0);
        }
    }

}

template<size_t t>
void hash_z_shares(AES& aes, array<array<block,1>,t>& z_vec_shares, array<block,t>& hashed_z_shares) {

    for (size_t i=0;i < t;i++) {
        
        hashed_z_shares[i] = aes.hashBlock(z_vec_shares[i][0]);

    }

}

template<uint64_t M>
inline static uint64_t compMinX(uint64_t a0, ZN<M> b0, uint64_t delta) {
    const uint64_t delta_squared = delta*delta;

    uint64_t x_squared = delta_squared - (a0 - b0.to_uint64_t())*(a0 - b0.to_uint64_t());

     std::ceil(std::sqrt((double) x_squared) + 1.0);

}

template<size_t tr, size_t ts, size_t d, uint16_t twotol, uint8_t delta, uint64_t M>
static Proto sender_compute_h_shares(OprfSender* oprfSender,
                                    Socket& sock, 
                                    PRNG& prng,
                                    array<point,ts>& ordIndexSet,
                                    array<array<ZN<twotol>,d>,ts>& in_vals,
                                    array<array<ZN<M>,d>,ts>& h_shares) {
    
    static_assert(M == d*(delta + 1) + 1,"the following identity must be fulfilled: M = d*(delta + 1) + 1");

    MC_BEGIN(Proto, oprfSender, &sock, &prng, &ordIndexSet,&in_vals,&h_shares,
             zero_shares = (array<array<ZN<twotol>,d>,ts>*) nullptr,
             msg_vecs = (array<array<array<ZN<M>,twotol>,d>,ts>*) nullptr,
             bsotSender = (SpBSOTSender<tr,ts,d,twotol,M>*) nullptr,
             diff = int32_t(0),
             abs_diff = int32_t(0));
        zero_shares = new array<array<ZN<twotol>,d>,ts>();
        msg_vecs = new array<array<array<ZN<M>,twotol>,d>,ts>();
        bsotSender = new SpBSOTSender<tr,ts,d,twotol,M>(prng, oprfSender);

        for (size_t i=0;i < ts;i++) {
            array<array<ZN<M>,twotol>,d>& msg_vec_batch = msg_vecs->at(i);

            for (size_t j=0;j < d;j++) {
                array<ZN<M>,twotol>& msg_vec = msg_vec_batch[j];
                
                for (int32_t h=0;h < twotol;h++) {

                    diff = h-((int32_t) in_vals[i][j].to_uint64_t());
                    abs_diff = (int32_t) abs(diff);
                    msg_vec[h] = (abs_diff <= delta) ? ZN<M>(abs_diff) : ZN<M>(delta+1);
                    
                }

            }

        }

        gen_zero_shares<ts,d,twotol>(*zero_shares);

        MC_AWAIT(bsotSender->send(sock, ordIndexSet, *msg_vecs, *zero_shares, h_shares));

        delete zero_shares;
        delete msg_vecs;
        delete bsotSender;

    MC_END();

}


template<size_t ts, size_t tr, size_t d, uint16_t twotol, uint8_t delta, uint64_t M>
static Proto recvr_compute_h_shares(OprfReceiver* oprfReceiver,
                                   Socket& sock, 
                                   PRNG& prng,
                                   array<point,tr>& ordIndexSet,
                                   array<array<ZN<twotol>,d>,tr>& in_vals,
                                   array<array<ZN<M>,d>,tr>& h_shares) {
    
    static_assert(M == d*(delta + 1) + 1,"the following identity must be fulfilled: M = d*(delta + 1) + 1");
    
    MC_BEGIN(Proto, oprfReceiver, &sock, &prng, &ordIndexSet, &in_vals, &h_shares,
             bsotReceiver = (SpBSOTReceiver<ts, tr,d,twotol,M>*) nullptr);
        
        bsotReceiver = new SpBSOTReceiver<ts, tr,d,twotol,M>(prng, oprfReceiver);
       
        MC_AWAIT(bsotReceiver->receive(sock, ordIndexSet, in_vals, h_shares));

        delete bsotReceiver;

    MC_END();

}