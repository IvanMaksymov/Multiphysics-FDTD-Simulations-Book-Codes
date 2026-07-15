#include "ap_fixed.h"

#define NX 32
#define NY 32

typedef ap_fixed<16,6> real_t;

void fdtd_tmz_step(
    real_t Ez[NX][NY],
    real_t Hx[NX][NY],
    real_t Hy[NX][NY],
    const real_t ce,
    const real_t chx,
    const real_t chy)
{
#pragma HLS INTERFACE bram storage_type=ram_t2p port=Ez
#pragma HLS INTERFACE bram storage_type=ram_t2p port=Hx
#pragma HLS INTERFACE bram storage_type=ram_t2p port=Hy

#pragma HLS INTERFACE s_axilite port=ce
#pragma HLS INTERFACE s_axilite port=chx
#pragma HLS INTERFACE s_axilite port=chy
#pragma HLS INTERFACE s_axilite port=return

// Hx
for(int i=0;i<NX;i++)
  for(int j=0;j<NY-1;j++){
#pragma HLS PIPELINE II=1
    Hx[i][j] -= chx * (Ez[i][j+1] - Ez[i][j]);
  }

// Hy
for(int i=0;i<NX-1;i++)
  for(int j=0;j<NY;j++){
#pragma HLS PIPELINE II=1
    Hy[i][j] += chy * (Ez[i+1][j] - Ez[i][j]);
  }

// Ez
for(int i=1;i<NX;i++)
  for(int j=1;j<NY;j++){
#pragma HLS PIPELINE II=1
    Ez[i][j] += ce * ((Hy[i][j] - Hy[i-1][j]) -
                      (Hx[i][j] - Hx[i][j-1]));
  }

// PEC
for(int i=0;i<NX;i++){
#pragma HLS PIPELINE II=1
  Ez[i][0]=0;
  Ez[i][NY-1]=0;
}
for(int j=0;j<NY;j++){
#pragma HLS PIPELINE II=1
  Ez[0][j]=0;
  Ez[NX-1][j]=0;
}
}

