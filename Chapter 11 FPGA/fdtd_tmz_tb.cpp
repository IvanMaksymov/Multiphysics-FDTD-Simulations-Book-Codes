//============================================================
// Minimal TMz FDTD FPGA C simulation testbench
//
// Host-side simulation model
//
// One kernel call = one timestep
//
//============================================================


#include <stdio.h>
#include <math.h>

#include "ap_fixed.h"



#define NX 32
#define NY 32



typedef ap_fixed<16,6> real_t;



//------------------------------------------------------------
// Kernel declaration
//------------------------------------------------------------

void fdtd_tmz_step
(
    real_t Ez[NX][NY],
    real_t Hx[NX][NY],
    real_t Hy[NX][NY],

    const real_t ce,
    const real_t chx,
    const real_t chy
);



//============================================================
// Main
//============================================================

int main()
{


    //--------------------------------------------------------
    // Field arrays
    //--------------------------------------------------------

    static real_t Ez[NX][NY];
    static real_t Hx[NX][NY];
    static real_t Hy[NX][NY];



    //--------------------------------------------------------
    // Initialise fields
    //--------------------------------------------------------

    for(int i=0;i<NX;i++)
    {
        for(int j=0;j<NY;j++)
        {
            Ez[i][j]=0;
            Hx[i][j]=0;
            Hy[i][j]=0;
        }
    }



    //--------------------------------------------------------
    // FDTD coefficients
    //--------------------------------------------------------

    const real_t ce  = 0.5;
    const real_t chx = 0.5;
    const real_t chy = 0.5;



    //--------------------------------------------------------
    // Time stepping
    //--------------------------------------------------------

    const int NSTEPS = 100;



    for(int n=0;n<NSTEPS;n++)
    {


        //----------------------------------------------------
        // Soft Gaussian source
        //----------------------------------------------------

        float pulse =
            expf(
                -0.5f *
                powf((n-20)/6.0f,2.0f)
            );


        Ez[NX/2][NY/2] = pulse;



        //----------------------------------------------------
        // FPGA kernel call
        //----------------------------------------------------

        fdtd_tmz_step
        (
            Ez,
            Hx,
            Hy,

            ce,
            chx,
            chy
        );



        if(n%10==0)
        {

            printf(
                "step %d : Ez(%d,%d) = %f\n",
                n,
                NX/2,
                NY/2,
                (float)Ez[NX/2][NY/2]
            );

        }

    }



    //--------------------------------------------------------
    // Save Ez field
    //--------------------------------------------------------

    FILE *fp=fopen("Ez.dat","w");


    if(fp==NULL)
    {
        printf("ERROR: cannot open Ez.dat\n");
        return 1;
    }



    for(int i=0;i<NX;i++)
    {
        for(int j=0;j<NY;j++)
        {

            fprintf(
                fp,
                "%d %d %.6f\n",
                i,
                j,
                (float)Ez[i][j]
            );

        }

        fprintf(fp,"\n");

    }


    fclose(fp);



    printf("\nSimulation finished\n");


    return 0;

}
