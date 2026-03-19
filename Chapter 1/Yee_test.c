#include <stdio.h>
#include <math.h>

#define NX 2000              /* grid size */
#define NSTEPS 2000          /* time steps */
#define OUTPUT_STEP 20       /* output interval */

#define C0 299792458.0f
#define EPS0 8.854e-12f
#define MU0 (1.0f/(EPS0*C0*C0))

int main() {

    float dx = 1.0e-3f;
    float dt = dx/(2.0f*C0);

    float Ez[NX] = {0.0f};
    float Hy[NX] = {0.0f};

    int i, t;
    int src = NX/4;   /* source position */

    FILE *f = fopen("field.dat","w");

    for(t = 0; t < NSTEPS; t++) {

        /* --- update H field --- */
        for(i = 0; i < NX-1; i++) {
            Hy[i] += (dt/(MU0*dx)) * (Ez[i+1] - Ez[i]);
        }

        /* --- update E field --- */
        for(i = 1; i < NX; i++) {
            Ez[i] += (dt/(EPS0*dx)) * (Hy[i] - Hy[i-1]);
        }

        /* --- PEC boundaries --- */
        Ez[0] = 0.0f;
        Ez[NX-1] = 0.0f;

        /* --- half-sine pulse source --- */
        if(t < 100) {
            float pulse = sinf(3.1415926f * t / 100.0f);
            Ez[src] += pulse;
        }

        /* --- output --- */
        if(t % OUTPUT_STEP == 0) {
            for(i = 0; i < NX; i++) {
                fprintf(f, "%e ", Ez[i]);
            }
            fprintf(f, "\n");
        }
    }

    fclose(f);

    return 0;
}
