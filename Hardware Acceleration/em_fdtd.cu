// em_fdtd.cu
//
// 2D TEz Electromagnetic FDTD on GPU (Yee grid)
//
// Fields:
//   Ez(x,y)  - electric field (out of plane)
//   Hx(x,y)  - magnetic field (x-direction)
//   Hy(x,y)  - magnetic field (y-direction)
//
// Scheme:
//   Standard Yee FDTD for TEz in 2D
//
// Compile:
//   nvcc em_fdtd.cu -O3 -use_fast_math -o em_fdtd
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cuda_runtime.h>

// ----------------------------------------------------
// Grid and time-stepping parameters
// ----------------------------------------------------

#define NX 1024          // number of cells in x
#define NY 1024          // number of cells in y

#define NSTEPS 500      // number of time steps

// ----------------------------------------------------
// Physical constants (host-side)
// ----------------------------------------------------

const float eps0 = 8.854187817e-12f;   // vacuum permittivity
const float mu0  = 1.2566370614e-6f;   // vacuum permeability
const float c0   = 299792458.0f;       // speed of light

// ----------------------------------------------------
// Spatial steps (host-side)
// ----------------------------------------------------

const float dx = 1.0e-3f;              // cell size in x
const float dy = 1.0e-3f;              // cell size in y

// CFL time step (host-side)
const float dt =
    0.99f /
    (c0 * sqrtf(1.0f/(dx*dx) + 1.0f/(dy*dy)));

// ----------------------------------------------------
// Device-side constants (visible in kernels)
// ----------------------------------------------------

// These are copies of the host constants, stored in GPU constant memory.
__constant__ float dt_d;
__constant__ float dx_d;
__constant__ float dy_d;
__constant__ float eps0_d;
__constant__ float mu0_d;

// ----------------------------------------------------
// CUDA block size
// ----------------------------------------------------

#define BX 16
#define BY 16

// ----------------------------------------------------
// Index helper (device)
// ----------------------------------------------------
//
// Flatten (x,y) into a single index for 1D arrays of size NX*NY.
//

__device__ __forceinline__
int IDX(int x, int y)
{
    return y * NX + x;
}

// ----------------------------------------------------
// Update Hx kernel
//
// Hx(i,j) = Hx(i,j) - (dt/mu0) * dEz/dy
//
// Uses forward difference in y for Ez.
// ----------------------------------------------------

__global__
void updateHx(
    float *Hx,
    const float *Ez
)
{
    int x = blockIdx.x * blockDim.x + threadIdx.x;
    int y = blockIdx.y * blockDim.y + threadIdx.y;

    // Hx is defined for all y except the last (forward difference)
    if (x >= NX)   return;
    if (y >= NY-1) return;

    int i = IDX(x, y);

    float dEzdy =
        (Ez[IDX(x, y+1)] - Ez[i]) / dy_d;

    Hx[i] -= (dt_d / mu0_d) * dEzdy;
}

// ----------------------------------------------------
// Update Hy kernel
//
// Hy(i,j) = Hy(i,j) + (dt/mu0) * dEz/dx
//
// Uses forward difference in x for Ez.
// ----------------------------------------------------

__global__
void updateHy(
    float *Hy,
    const float *Ez
)
{
    int x = blockIdx.x * blockDim.x + threadIdx.x;
    int y = blockIdx.y * blockDim.y + threadIdx.y;

    // Hy is defined for all x except the last (forward difference)
    if (x >= NX-1) return;
    if (y >= NY)   return;

    int i = IDX(x, y);

    float dEzdx =
        (Ez[IDX(x+1, y)] - Ez[i]) / dx_d;

    Hy[i] += (dt_d / mu0_d) * dEzdx;
}

// ----------------------------------------------------
// Soft source kernel
//
// Adds a Gaussian pulse at the center of the grid.
// ----------------------------------------------------

__global__
void sourceKernel(
    float *Ez,
    int step
)
{
    int cx = NX / 2;
    int cy = NY / 2;

    float t0     = 40.0f;   // pulse center in time
    float spread = 12.0f;   // pulse width

    float arg = (step - t0) / spread;
    float pulse = expf(-arg * arg);

    Ez[IDX(cx, cy)] += pulse;
}

// ----------------------------------------------------
// PEC boundary kernel
//
// Enforces Ez = 0 at all boundaries (perfect electric conductor).
// ----------------------------------------------------

__global__
void applyPEC(float *Ez)
{
    int x = blockIdx.x * blockDim.x + threadIdx.x;

    // Top and bottom rows
    if (x < NX)
    {
        Ez[IDX(x, 0)]    = 0.0f;
        Ez[IDX(x, NY-1)] = 0.0f;
    }

    // Left and right columns
    if (x < NY)
    {
        Ez[IDX(0,     x)] = 0.0f;
        Ez[IDX(NX-1,  x)] = 0.0f;
    }
}

// ----------------------------------------------------
// Error checker (host)
// ----------------------------------------------------

void check(cudaError_t err)
{
    if (err != cudaSuccess)
    {
        printf("CUDA ERROR: %s\n", cudaGetErrorString(err));
        exit(1);
    }
}

// ----------------------------------------------------
// Update Ez kernel
//
// ∂Ez/∂t = (1/ε0) * ( ∂Hy/∂x - ∂Hx/∂y )
//
// Uses backward differences for Hx and Hy.
// ----------------------------------------------------

__global__
void updateEz(
    float *Ez,
    const float *Hx,
    const float *Hy
)
{
    int x = blockIdx.x * blockDim.x + threadIdx.x;
    int y = blockIdx.y * blockDim.y + threadIdx.y;

    // Interior points only (need neighbors in -x and -y)
    if (x < 1 || x >= NX-1) return;
    if (y < 1 || y >= NY-1) return;

    int i = IDX(x, y);

    float dHydx =
        (Hy[i] - Hy[IDX(x-1, y)]) / dx_d;

    float dHxdy =
        (Hx[i] - Hx[IDX(x, y-1)]) / dy_d;

    Ez[i] += (dt_d / eps0_d) * (dHydx - dHxdy);
}

// ----------------------------------------------------
// MAIN (host)
// ----------------------------------------------------

int main()
{
    // Total number of cells
    size_t bytes = NX * NY * sizeof(float);

    // Device pointers
    float *dEz;
    float *dHx;
    float *dHy;

    // Allocate device memory
    check(cudaMalloc(&dEz, bytes));
    check(cudaMalloc(&dHx, bytes));
    check(cudaMalloc(&dHy, bytes));

    // Initialize fields to zero
    check(cudaMemset(dEz, 0, bytes));
    check(cudaMemset(dHx, 0, bytes));
    check(cudaMemset(dHy, 0, bytes));

    // Copy physical constants to device constant memory
    check(cudaMemcpyToSymbol(dt_d,   &dt,   sizeof(float)));
    check(cudaMemcpyToSymbol(dx_d,   &dx,   sizeof(float)));
    check(cudaMemcpyToSymbol(dy_d,   &dy,   sizeof(float)));
    check(cudaMemcpyToSymbol(eps0_d, &eps0, sizeof(float)));
    check(cudaMemcpyToSymbol(mu0_d,  &mu0,  sizeof(float)));

    // Thread/block configuration for field updates
    dim3 threads(BX, BY);
    dim3 blocks(
        (NX + BX - 1) / BX,
        (NY + BY - 1) / BY
    );

    // 1D configuration for PEC kernel
    int nBlocks = (NX > NY ? NX : NY) + 255;
    nBlocks /= 256;

    printf("\n");
    printf("CUDA 2D TEz FDTD\n");
    printf("---------------------------\n");
    printf("Grid : %d x %d\n", NX, NY);
    printf("Steps: %d\n", NSTEPS);
    printf("dx   : %e m\n", dx);
    printf("dy   : %e m\n", dy);
    printf("dt   : %e s\n", dt);
    printf("\n");

    // Timing events
    cudaEvent_t start, stop;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);

    cudaEventRecord(start);

    // --------------------------------------------------
    // Time stepping loop
    // --------------------------------------------------

    for (int n = 0; n < NSTEPS; n++)
    {
        // Update magnetic fields
        updateHx<<<blocks, threads>>>(dHx, dEz);
        updateHy<<<blocks, threads>>>(dHy, dEz);

        // Update electric field
        updateEz<<<blocks, threads>>>(dEz, dHx, dHy);

        // Inject source
        sourceKernel<<<1, 1>>>(dEz, n);

        // Apply PEC boundaries
        applyPEC<<<nBlocks, 256>>>(dEz);

        if ((n % 500) == 0)
            printf("Step %d\n", n);
    }

    cudaEventRecord(stop);
    cudaEventSynchronize(stop);

    float ms;
    cudaEventElapsedTime(&ms, start, stop);

    printf("\nRuntime %.3f ms\n", ms);

    // --------------------------------------------------
    // Copy Ez field back to host
    // --------------------------------------------------

    float *Ez = (float*)malloc(bytes);

    check(cudaMemcpy(
        Ez,
        dEz,
        bytes,
        cudaMemcpyDeviceToHost
    ));

    // Free device memory and destroy events
    cudaFree(dEz);
    cudaFree(dHx);
    cudaFree(dHy);

    cudaEventDestroy(start);
    cudaEventDestroy(stop);

    // ----------------------------------------------------
    // Save Ez field to text file
    // ----------------------------------------------------

    FILE *fp = fopen("Ez.dat", "w");

    if (fp == NULL)
    {
        printf("Cannot open output file Ez.dat.\n");
        free(Ez);
        return 1;
    }

    for (int j = 0; j < NY; j++)
    {
        for (int i = 0; i < NX; i++)
        {
            fprintf(fp, "%e ", Ez[j * NX + i]);
        }
        fprintf(fp, "\n");
    }

    fclose(fp);

    printf("Saved field to Ez.dat\n");

    // ----------------------------------------------------
    // Optional PGM image of Ez
    // ----------------------------------------------------

    FILE *img = fopen("Ez.pgm", "wb");

    if (img)
    {
        fprintf(img, "P5\n%d %d\n255\n", NX, NY);

        // Find max absolute Ez for normalization
        float vmax = 0.0f;
        for (int i = 0; i < NX * NY; i++)
        {
            float a = fabsf(Ez[i]);
            if (a > vmax) vmax = a;
        }
        if (vmax == 0.0f) vmax = 1.0f;

        // Map Ez to grayscale [0,255]
        for (int i = 0; i < NX * NY; i++)
        {
            float x = Ez[i] / vmax;
            x = 0.5f * (x + 1.0f);   // shift to [0,1]

            if (x < 0.0f) x = 0.0f;
            if (x > 1.0f) x = 1.0f;

            unsigned char c =
                (unsigned char)(255.0f * x);

            fwrite(&c, 1, 1, img);
        }

        fclose(img);
        printf("Saved image Ez.pgm\n");
    }

    // ----------------------------------------------------
    // Cleanup
    // ----------------------------------------------------

    free(Ez);

    printf("\nFinished successfully.\n");

    return 0;
}

