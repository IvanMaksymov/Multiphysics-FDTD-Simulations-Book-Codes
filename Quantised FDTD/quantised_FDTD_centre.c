#include <stdio.h>
#include <stdint.h>
#include <math.h>

#define NX      80
#define NY      80
#define STEPS   500
#define CFL     0.1f

/*---------------------------------------------
  Quantiser threshold
---------------------------------------------*/
#define THRESHOLD         4

/*---------------------------------------------
  Ternary EM fields (-1,0,+1)
    Ez: electric field (out of plane)
    Hx, Hy: magnetic field components
---------------------------------------------*/
static int8_t Ez[NX][NY];
static int8_t Hx[NX][NY];
static int8_t Hy[NX][NY];

/*---------------------------------------------
  Internal accumulators (integer)
---------------------------------------------*/
static int32_t aEz[NX][NY];
static int32_t aHx[NX][NY];
static int32_t aHy[NX][NY];

/*---------------------------------------------
  Floating-point reference EM solver
---------------------------------------------*/
static float fp_Ez[NX][NY];
static float fp_Hx[NX][NY];
static float fp_Hy[NX][NY];

static float fp_aEz[NX][NY];
static float fp_aHx[NX][NY];
static float fp_aHy[NX][NY];

/*---------------------------------------------
  Residual-preserving quantiser
---------------------------------------------*/
static inline int8_t quantise(int32_t *acc, int8_t old)
{
    if (*acc >= THRESHOLD)
    {
        *acc -= THRESHOLD;
        return 1;
    }

    if (*acc <= -THRESHOLD)
    {
        *acc += THRESHOLD;
        return -1;
    }

    return old;
}

/*---------------------------------------------
  Initialisation
---------------------------------------------*/
void init_fields(void)
{
    int i,j;

    for(i=0;i<NX;i++)
    {
        for(j=0;j<NY;j++)
        {
            Ez[i][j] = 0;
            Hx[i][j] = 0;
            Hy[i][j] = 0;

            aEz[i][j] = 0;
            aHx[i][j] = 0;
            aHy[i][j] = 0;

            fp_Ez[i][j] = 0.0f;
            fp_Hx[i][j] = 0.0f;
            fp_Hy[i][j] = 0.0f;

            fp_aEz[i][j] = 0.0f;
            fp_aHx[i][j] = 0.0f;
            fp_aHy[i][j] = 0.0f;
        }
    }

}

/*---------------------------------------------
  Harmonic Ez source
---------------------------------------------*/
#define SOURCE_AMPLITUDE  4.0f
#define SOURCE_FREQUENCY  0.01f

int source(int n)
{
    float s = SOURCE_AMPLITUDE *
              sinf(2.0f * 3.14159265f *
                   SOURCE_FREQUENCY * n);

    return (int)s;
}

float fp_source(int n)
{
    return SOURCE_AMPLITUDE *
           sinf(2.0f * 3.14159265f *
                SOURCE_FREQUENCY * n);
}

/*---------------------------------------------
  Dump Ez field (ternary)
---------------------------------------------*/
void dump_snapshot(int step)
{
    int i,j;
    char filename[64];

    sprintf(filename,"snapshot_%04d.txt",step);

    FILE *fp=fopen(filename,"w");

    if(fp==NULL)
    {
        perror("file");
        return;
    }

    fprintf(fp,"STEP %d\n\n",step);

    for(i=0;i<NX;i++)
    {
        for(j=0;j<NY;j++)
            fprintf(fp,"%6d ",Ez[i][j]);

        fprintf(fp,"\n");
    }

    fclose(fp);
}

/*---------------------------------------------
  Dump Ez field (floating point)
---------------------------------------------*/
void dump_fp_snapshot(int step)
{
    int i,j;
    char filename[64];

    sprintf(filename,"fp_snapshot_%04d.txt",step);

    FILE *fp=fopen(filename,"w");

    if(fp==NULL)
    {
        perror("file");
        return;
    }

    fprintf(fp,"STEP %d\n\n",step);

    for(i=0;i<NX;i++)
    {
        for(j=0;j<NY;j++)
            fprintf(fp,"%12.6f ",fp_Ez[i][j]);

        fprintf(fp,"\n");
    }

    fclose(fp);
}

/*---------------------------------------------
  Main
---------------------------------------------*/
int main(void)
{
    int n;
    int i,j;

    const int sx = NX/2;
    const int sy = NY/2;

    init_fields();

    for(n=0;n<STEPS;n++)
    {
        /*-------------------------------------
          Magnetic field update (Hx, Hy)
          TMz Yee-like:
            Hx^{n+1/2} += -dEz/dy
            Hy^{n+1/2} +=  dEz/dx
        -------------------------------------*/
        for(i=0;i<NX-1;i++)
        {
            for(j=0;j<NY-1;j++)
            {
                /* integer accumulators */
                aHx[i][j] += Ez[i][j] - Ez[i][j+1];   /* -dEz/dy */
                aHy[i][j] += Ez[i+1][j] - Ez[i][j];   /*  dEz/dx */

                Hx[i][j] = quantise(&aHx[i][j], Hx[i][j]);
                Hy[i][j] = quantise(&aHy[i][j], Hy[i][j]);
            }
        }

        /*-------------------------------------
          Electric field update (Ez)
          Ez^{n+1} += dHy/dx - dHx/dy
        -------------------------------------*/
        for(i=1;i<NX-1;i++)
        {
            for(j=1;j<NY-1;j++)
            {
                aEz[i][j] +=
                    (Hy[i][j]   - Hy[i-1][j])   /* dHy/dx */
                  - (Hx[i][j]   - Hx[i][j-1]);  /* dHx/dy */

                Ez[i][j] = quantise(&aEz[i][j], Ez[i][j]);
            }
        }

        /*-------------------------------------
          Source (Ez at sx,sy)
        -------------------------------------*/
        aEz[sx][sy] += source(n);
        Ez[sx][sy] = quantise(&aEz[sx][sy], Ez[sx][sy]);

        /*-------------------------------------
          Boundary conditions (Ez=0)
        -------------------------------------*/
        for(i=0;i<NX;i++)
        {
            Ez[i][0]    = 0;
            Ez[i][NY-1] = 0;
        }
        for(j=0;j<NY;j++)
        {
            Ez[0][j]    = 0;
            Ez[NX-1][j] = 0;
        }

        /*-------------------------------------
          Floating-point reference EM FDTD
        -------------------------------------*/
        /* H update */
        for(i=0;i<NX-1;i++)
        {
            for(j=0;j<NY-1;j++)
            {
                fp_aHx[i][j] += CFL * (fp_Ez[i][j] - fp_Ez[i][j+1]);
                fp_aHy[i][j] += CFL * (fp_Ez[i+1][j] - fp_Ez[i][j]);

                fp_Hx[i][j] = fp_aHx[i][j];
                fp_Hy[i][j] = fp_aHy[i][j];
            }
        }

        /* Ez update */
        for(i=1;i<NX-1;i++)
        {
            for(j=1;j<NY-1;j++)
            {
                fp_aEz[i][j] += CFL *
                    ((fp_Hy[i][j]   - fp_Hy[i-1][j])   /* dHy/dx */
                   - (fp_Hx[i][j]   - fp_Hx[i][j-1])); /* dHx/dy */

                fp_Ez[i][j] = fp_aEz[i][j];
            }
        }

        /* identical source */
        fp_aEz[sx][sy] += fp_source(n);
        fp_Ez[sx][sy] = fp_aEz[sx][sy];

        /* boundaries Ez=0 */
        for(i=0;i<NX;i++)
        {
            fp_Ez[i][0]    = 0.0f;
            fp_Ez[i][NY-1] = 0.0f;
        }
        for(j=0;j<NY;j++)
        {
            fp_Ez[0][j]    = 0.0f;
            fp_Ez[NX-1][j] = 0.0f;
        }

        long energy=0;

        for(i=0;i<NX;i++)
        {
            for(j=0;j<NY;j++)
            {
                energy += Ez[i][j]*Ez[i][j];
                energy += Hx[i][j]*Hx[i][j];
                energy += Hy[i][j]*Hy[i][j];
            }
        }

        printf("%d energy=%ld\n",n,energy);

        if((n % 20) == 0)
        {  
            dump_snapshot(n);
            dump_fp_snapshot(n);
        }
    }

    return 0;
}

