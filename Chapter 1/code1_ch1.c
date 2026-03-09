//To complile: gcc -O3 -o aef code1.c -lm
//Ro run: ./aef
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define NX 200			 // Number of spatial cells
#define NY 200			 // Number of spatial cells	
#define NSTEPS 500		 // Number of time steps
#define DX 0.01f		 // Spatial step (meters)
#define DY 0.01f		 // Spatial step (meters)
#define EPS0 8.854e-12f  // Permittivity of free space
#define MU0 1.257e-6f    // Permeability of free space
#define C 299792458.0f   // Speed of light


#define CYL_X 100
#define CYL_Y 100
#define CYL_RADIUS 20

int main() {
    int i, j, n;

    float Ez[NX][NY] = {0.0f};
    float Hx[NX][NY] = {0.0f};
    float Hy[NX][NY] = {0.0f};

    // CFL stability: 2D
    // Time step (seconds)
	float dt = 1.0f / (C * sqrt(1.0f/(DX*DX) + 1.0f/(DY*DY)));
	printf("Time step (CFL-limited) = %e s\n", dt);

    float ce = dt / (EPS0 * DX);
    float chx = dt / (MU0 * DY);
    float chy = dt / (MU0 * DX);

    char filename[256];
    

    for(n=0; n<NSTEPS; n++) {

        // Update Hx
        for(i=0; i<NX; i++)
            for(j=0; j<NY-1; j++)
                Hx[i][j] -= chx * (Ez[i][j+1] - Ez[i][j]);

        // Update Hy
        for(i=0; i<NX-1; i++)
            for(j=0; j<NY; j++)
                Hy[i][j] += chy * (Ez[i+1][j] - Ez[i][j]);

        // Source
        Ez[10][NY/2] += expf(-0.5f * ((n-30)/10.0f)*((n-30)/10.0f));

        // Update Ez (TMz)
        for(i=1; i<NX; i++) {
            for(j=1; j<NY; j++) {
                int dx = i - CYL_X;
                int dy = j - CYL_Y;
                if(dx*dx + dy*dy <= CYL_RADIUS*CYL_RADIUS)
                    Ez[i][j] = 0.0f;  // PEC cylinder
                else
                    Ez[i][j] += ce * ((Hy[i][j] - Hy[i-1][j]) - (Hx[i][j] - Hx[i][j-1]));
            }
        }

        // Mur ABC
        for(i=1; i<NX-1; i++) {
            Ez[i][0] = Ez[i][1];
            Ez[i][NY-1] = Ez[i][NY-2];
        }
        for(j=1; j<NY-1; j++) {
            Ez[0][j] = Ez[1][j];
            Ez[NX-1][j] = Ez[NX-2][j];
        }

        // Output every 10 steps
        if(n % 10 == 0) {
            sprintf(filename, "Ez_step_%04d.dat", n);
            FILE *fp = fopen(filename, "w");
            for(j=0; j<NY; j++) {
                for(i=0; i<NX; i++) {
                    fprintf(fp, "%e ", Ez[i][j]);
                }
                fprintf(fp, "\n");
            }
            fclose(fp);
        }
    }

    return 0;
}
