#include <stdio.h>
#include <math.h>

/* --- Grid and simulation parameters --- */
#define NX 2000           // number of spatial points
#define NSTEPS 2000       // number of time steps
#define OUTPUT_STEP 20    // output interval

#define C0 299792458.0f
#define EPS0 8.854e-12f
#define MU0 (1.0f/(EPS0*C0*C0))

int main() {

    /* --- User-defined wavelength and resolution --- */
    float lambda = 0.2f;   // wavelength in meters
    int N = 20;            // points per wavelength

    /* --- Derived spatial and temporal steps --- */
    float dx = lambda / N;
    float dt = 0.5f * dx / C0;   // CFL condition (S=0.5)

    /* --- Pulse temporal length --- */
    int pulse_steps = 2 * N;     // half-sine pulse spans half wavelength
                                  // full sine would span 2*half_sine

    /* --- Fields --- */
    float Ez[NX];
    float Hy[NX];

    for(int i = 0; i < NX; i++) {
        Ez[i] = 0.0f;
        Hy[i] = 0.0f;
    }

    int src = NX/4;   // source location

    FILE *f = fopen("field.dat","w");

    for(int t = 0; t < NSTEPS; t++) {

        /* --- update H field --- */
        for(int i = 0; i < NX-1; i++)
            Hy[i] += (dt/(MU0*dx)) * (Ez[i+1] - Ez[i]);

        /* --- update E field --- */
        for(int i = 1; i < NX; i++)
            Ez[i] += (dt/(EPS0*dx)) * (Hy[i] - Hy[i-1]);

        /* --- PEC boundaries --- */
        Ez[0] = 0.0f;
        Ez[NX-1] = 0.0f;

        /* --- half-sine pulse source --- */
        if(t < pulse_steps/2) {
            float pulse = sinf(M_PI * t / (pulse_steps/2));
            Ez[src] += pulse;
        }

        /* --- output --- */
        if(t % OUTPUT_STEP == 0) {
            for(int i = 0; i < NX; i++)
                fprintf(f, "%e ", Ez[i]);
            fprintf(f, "\n");
        }
    }

    fclose(f);

    printf("Simulation finished. dx = %e m, dt = %e s, pulse_steps = %d\n", dx, dt, pulse_steps);

    return 0;
}
