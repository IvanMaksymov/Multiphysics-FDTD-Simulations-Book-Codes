#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>

#define NX 500
#define NZ 500
#define MAX_ITER 10000

// Memory-safe 2D mapping (Row-Major to match C)
#define TXX(x,z)    tau_xx[(x)*NZ + (z)]
#define TZZ(x,z)    tau_zz[(x)*NZ + (z)]
#define TXZ(x,z)    tau_xz[(x)*(NZ-1) + (z)]
#define UX(x,z)     ux[(x)*NZ + (z)]
#define UZ(x,z)     uz[(x)*(NZ-1) + (z)]
#define TXX_MOD(x,z) tau_xx_mod[(x)*NZ + (z)]
#define UZ_MOD(x,z)  uz_mod[(x)*(NZ-1) + (z)]
#define UX_MOD(x,z)  ux_mod[(x)*NZ + (z)]

int main() {
    int xind, zind, n;
    float pi = 3.14159265f, freq = 1.7e6f, cmax = 1600.0f;
    float dx = 5e-6f, dz = 5e-6f;
    float dt = 0.95f * (dx / (sqrtf(2.0f) * cmax));
    
    // Allocate all arrays on heap to avoid stack overflow
    float *l = calloc(NX*NZ, sizeof(float));
    float *m = calloc(NX*NZ, sizeof(float));
    float *rho = calloc(NX*NZ, sizeof(float));
    float *b = calloc(NX*NZ, sizeof(float));
    float *tau_xx = calloc(NX*NZ, sizeof(float));
    float *tau_zz = calloc(NX*NZ, sizeof(float));
    float *tau_xx_n = calloc(NX*NZ, sizeof(float));
    float *tau_zz_n = calloc(NX*NZ, sizeof(float));
    float *tau_xz = calloc((NX-1)*(NZ-1), sizeof(float));
    float *ux = calloc((NX-1)*NZ, sizeof(float));
    float *uz = calloc(NX*(NZ-1), sizeof(float));
    
    float complex *ux_mod = calloc((NX-1)*NZ, sizeof(float complex));
    float complex *uz_mod = calloc(NX*(NZ-1), sizeof(float complex));

    int center_x = NX/2;
    int center_z = NZ/2;

    // 1. Initialize Background (Water)
    for(int i=0; i<NX*NZ; i++) {
        rho[i] = 1000.0f;
        l[i] = 1000.0f * powf(1540.0f, 2.0f);
    }

    // 2. Define Waveguide Structure (Circle)
    for(xind=0; xind<NX; xind++) {
        for(zind=0; zind<NZ; zind++) {
            float dist = sqrtf(powf(center_x-xind, 2) + powf(center_z-zind, 2));
            if(dist < (525e-6f/dx/2.0f)) {
                l[xind*NZ+zind] = 1100.0f * (powf(1400.0f, 2.0f) - 2.0f*powf(440.0f, 2.0f));
                m[xind*NZ+zind] = 1100.0f * powf(440.0f, 2.0f);
                rho[xind*NZ+zind] = 1100.0f;
            }
            if(dist < (250e-6f/dx/2.0f)) {
                l[xind*NZ+zind] = 1000.0f * powf(1540.0f, 2.0f);
                m[xind*NZ+zind] = 0.0f;
                rho[xind*NZ+zind] = 1000.0f;
            }
            b[xind*NZ+zind] = 1.0f / rho[xind*NZ+zind];
        }
    }

    printf("Starting FDTD loop. dt = %e\n", dt);

    // 3. Main FDTD Loop
    for(n=0; n <= MAX_ITER; n++) {
        
        // Update Velocities (UX, UZ)
        for(xind=1; xind<NX-1; xind++) {
            for(zind=1; zind<NZ-1; zind++) {
                UX(xind,zind) += dt * b[xind*NZ+zind] * ((TXX(xind,zind) - TXX(xind-1,zind))/dx + (TXZ(xind,zind) - TXZ(xind,zind-1))/dz);
                UZ(xind,zind) += dt * b[xind*NZ+zind] * ((TXZ(xind,zind) - TXZ(xind-1,zind))/dx + (TZZ(xind,zind) - TZZ(xind,zind-1))/dz);
            }
        }

        // Update Stresses (TXX, TZZ)
        for(xind=0; xind<NX-2; xind++) {
            for(zind=0; zind<NZ-2; zind++) {
                float dux = (UX(xind+1,zind) - UX(xind,zind))/dx;
                float duz = (UZ(xind,zind+1) - UZ(xind,zind))/dz;
                TXX(xind,zind) += dt * (l[xind*NZ+zind] + 2.0f*m[xind*NZ+zind]) * dux + dt * l[xind*NZ+zind] * duz;
                TZZ(xind,zind) += dt * (l[xind*NZ+zind] + 2.0f*m[xind*NZ+zind]) * duz + dt * l[xind*NZ+zind] * dux;
            }
        }

        // Update Shear Stress (TXZ)
        for(xind=0; xind<NX-1; xind++) {
            for(zind=0; zind<NZ-1; zind++) {
                TXZ(xind,zind) += dt * m[xind*NZ+zind] * ((UX(xind,zind+1) - UX(xind,zind))/dz + (UZ(xind+1,zind) - UZ(xind,zind))/dz);
            }
        }

        // Source Injection (Moved to center-top)
        TXX(center_x, 5) += sinf(2.0f * pi * freq * n * dt);
        TZZ(center_x, 5) += sinf(2.0f * pi * freq * n * dt);

        // Mur1 ABC (Note: not all edges of the computational domain are implemented)
        float c_abs = 1540.0f;
        float factor = (c_abs*dt - dx) / (c_abs*dt + dx);
        for(zind=0; zind<NZ; zind++) {
            TXX(0, zind) = tau_xx_n[1*NZ+zind] + factor*(TXX(1,zind) - tau_xx_n[0*NZ+zind]);
            TXX(NX-1, zind) = tau_xx_n[(NX-2)*NZ+zind] + factor*(TXX(NX-2,zind) - tau_xx_n[(NX-1)*NZ+zind]);
        }

        // Store history for ABC
        for(int i=0; i<NX*NZ; i++) { tau_xx_n[i] = tau_xx[i]; tau_zz_n[i] = tau_zz[i]; }

        // Accumulate Complex Fields for FFT/Phasor
        float complex phase = cexpf(I * n * dt * 2.0f * pi * freq);
        for(xind=0; xind<NX-1; xind++) {
            for(zind=0; zind<NZ; zind++) UX_MOD(xind,zind) += UX(xind,zind) * phase;
        }
        for(xind=0; xind<NX; xind++) {
            for(zind=0; zind<NZ-1; zind++) UZ_MOD(xind,zind) += UZ(xind,zind) * phase;
        }

        if (n % 100 == 0) {
            printf("Step %d | Time: %.2f us\n", n, (n*dt)/1e-6);
            FILE *f = fopen("abs_u.out", "w");
            for(zind=0; zind<NZ-1; zind++) {
                for(xind=0; xind<NX-1; xind++) {
                    float val = cabsf(UX_MOD(xind,zind)) + cabsf(UZ_MOD(xind,zind));
                    fprintf(f, "%e ", val);
                }
                fprintf(f, "\n");
            }
            fclose(f);
        }
    }

    return 0;
}
