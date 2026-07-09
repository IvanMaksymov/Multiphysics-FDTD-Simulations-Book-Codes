// em_fdtd_cpu.c
//
// 2D TEz Electromagnetic FDTD on CPU
//
// Fields:
//   Ez(x,y)
//   Hx(x,y)
//   Hy(x,y)
//
// Yee grid, same equations as CUDA version
//
// Compile:
//   gcc em_fdtd_cpu.c -O3 -lm -o em_fdtd_cpu
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// ----------------------------------------------------
// Grid size
// ----------------------------------------------------

#define NX 1024
#define NY 1024

#define NSTEPS 500

// ----------------------------------------------------
// Physical constants
// ----------------------------------------------------

const float eps0 = 8.854187817e-12f;
const float mu0  = 1.2566370614e-6f;
const float c0   = 299792458.0f;

// ----------------------------------------------------
// Spatial steps
// ----------------------------------------------------

const float dx = 1.0e-3f;
const float dy = 1.0e-3f;
float dt;

// ----------------------------------------------------
// Index helper
// ----------------------------------------------------

static inline int IDX(int x, int y)
{
    return y * NX + x;
}

// ----------------------------------------------------
// Update Hx
// ----------------------------------------------------

void updateHx(float *Hx, const float *Ez)
{
    for (int y = 0; y < NY - 1; y++)
    {
        for (int x = 0; x < NX; x++)
        {
            int i = IDX(x, y);

            float dEzdy =
                (Ez[IDX(x, y+1)] - Ez[i]) / dy;

            Hx[i] -= (dt / mu0) * dEzdy;
        }
    }
}

// ----------------------------------------------------
// Update Hy
// ----------------------------------------------------

void updateHy(float *Hy, const float *Ez)
{
    for (int y = 0; y < NY; y++)
    {
        for (int x = 0; x < NX - 1; x++)
        {
            int i = IDX(x, y);

            float dEzdx =
                (Ez[IDX(x+1, y)] - Ez[i]) / dx;

            Hy[i] += (dt / mu0) * dEzdx;
        }
    }
}

// ----------------------------------------------------
// Soft source (Gaussian pulse)
// ----------------------------------------------------

void addSource(float *Ez, int step)
{
    int cx = NX / 2;
    int cy = NY / 2;

    float t0 = 40.0f;
    float spread = 12.0f;

    float arg = (step - t0) / spread;
    float pulse = expf(-arg * arg);

    Ez[IDX(cx, cy)] += pulse;
}

// ----------------------------------------------------
// PEC boundaries (Ez = 0)
// ----------------------------------------------------

void applyPEC(float *Ez)
{
    // Top and bottom
    for (int x = 0; x < NX; x++)
    {
        Ez[IDX(x, 0)]    = 0.0f;
        Ez[IDX(x, NY-1)] = 0.0f;
    }

    // Left and right
    for (int y = 0; y < NY; y++)
    {
        Ez[IDX(0,     y)] = 0.0f;
        Ez[IDX(NX-1,  y)] = 0.0f;
    }
}

// ----------------------------------------------------
// Update Ez
// ----------------------------------------------------

void updateEz(float *Ez, const float *Hx, const float *Hy)
{
    for (int y = 1; y < NY - 1; y++)
    {
        for (int x = 1; x < NX - 1; x++)
        {
            int i = IDX(x, y);

            float dHydx =
                (Hy[i] - Hy[IDX(x-1, y)]) / dx;

            float dHxdy =
                (Hx[i] - Hx[IDX(x, y-1)]) / dy;

            Ez[i] += (dt / eps0) * (dHydx - dHxdy);
        }
    }
}

// ----------------------------------------------------
// MAIN
// ----------------------------------------------------

int main()
{
    size_t bytes = NX * NY * sizeof(float);

    float *Ez = (float*)malloc(bytes);
    float *Hx = (float*)malloc(bytes);
    float *Hy = (float*)malloc(bytes);

    if (!Ez || !Hx || !Hy)
    {
        printf("Memory allocation failed.\n");
        return 1;
    }

    // Initialize fields
    for (int i = 0; i < NX * NY; i++)
    {
        Ez[i] = 0.0f;
        Hx[i] = 0.0f;
        Hy[i] = 0.0f;
    }

	dt = 0.99f / (c0 * sqrtf(1.0f/(dx*dx) + 1.0f/(dy*dy)));

    printf("\nCPU 2D TEz FDTD\n");
    printf("---------------------------\n");
    printf("Grid : %d x %d\n", NX, NY);
    printf("Steps: %d\n", NSTEPS);
    printf("dt   : %e\n\n", dt);

    // Time stepping
    for (int n = 0; n < NSTEPS; n++)
    {
        updateHx(Hx, Ez);
        updateHy(Hy, Ez);
        updateEz(Ez, Hx, Hy);

        addSource(Ez, n);
        applyPEC(Ez);

        if ((n % 500) == 0)
            printf("Step %d\n", n);
    }

    // Save Ez field
    FILE *fp = fopen("Ez_cpu.dat", "w");
    if (!fp)
    {
        printf("Cannot open Ez_cpu.dat\n");
        return 1;
    }

    for (int y = 0; y < NY; y++)
    {
        for (int x = 0; x < NX; x++)
        {
            fprintf(fp, "%e ", Ez[IDX(x, y)]);
        }
        fprintf(fp, "\n");
    }

    fclose(fp);
    printf("Saved Ez.dat\n");

    free(Ez);
    free(Hx);
    free(Hy);

    printf("Finished successfully.\n");
    return 0;
}

