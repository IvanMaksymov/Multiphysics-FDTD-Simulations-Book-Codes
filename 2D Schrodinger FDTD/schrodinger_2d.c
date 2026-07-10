// gcc schrodinger_2d.c -O3 -lm -o schrodinger_2d
// ./schrodinger_2d

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define L 8.0
#define DY 0.02
#define N 401     
#define NSTEPS 1000

static inline int IDX(int x, int y){ return y*N + x; }

int main(void){
    // Dimensionless time step
    const double dt = DY * DY / 4.0;
    
    // Physics constants 
    // H = - (1/2)*Laplacian + (1/2)*V
    const double c = dt / (2.0 * DY * DY); 
    const double dt_half = dt / 2.0;

    double *R = calloc(N*N, sizeof(double));
    double *I = calloc(N*N, sizeof(double));
    double *Rn = calloc(N*N, sizeof(double));
    double *In = calloc(N*N, sizeof(double));
    double *V = calloc(N*N, sizeof(double));

    // 1. Setup Potential 
    double w = 0.6, s = 0.5, a = 0.3, v0 = 600.0;
    int j0 = round((L - w) / (2 * DY));
    int j1 = round((L + w) / (2 * DY));
    int i0 = round((L + s) / (2 * DY) + a / DY);
    int i1 = round((L + s) / (2 * DY));
    int i2 = round((L - s) / (2 * DY));
    int i3 = round((L - s) / (2 * DY) - a / DY);

    // Note: With a=0.0, this creates a perfectly solid barrier.
    for(int y=0; y<i3; y++)  for(int x=j0; x<j1; x++) V[IDX(x,y)] = v0;
    for(int y=i2; y<i1; y++) for(int x=j0; x<j1; x++) V[IDX(x,y)] = v0;
    for(int y=i0; y<N; y++)  for(int x=j0; x<j1; x++) V[IDX(x,y)] = v0;

    // 2. Wavepacket Initialization (Gaussian with momentum k)
    double x0_pos = L / 5.0;
    double y0_pos = L / 2.0;
    double sigma = 0.5;
    double k = 15.0 * M_PI;

    for(int y=0; y<N; y++){
        for(int x=0; x<N; x++){
            double X_pos = x * DY;
            double Y_pos = y * DY;
            double env = exp(-0.5 * (pow(X_pos - x0_pos, 2) + pow(Y_pos - y0_pos, 2)) / (sigma * sigma));
            
            R[IDX(x,y)] = env * cos(k * (X_pos - x0_pos));
            I[IDX(x,y)] = env * sin(k * (X_pos - x0_pos));
        }
    }

    // 3. Time Evolution 
    for(int step=1; step<NSTEPS; step++){
        
        // Update Real Component
        for(int y=1; y<N-1; y++){
            for(int x=1; x<N-1; x++){
                int idx = IDX(x,y);
                double lapI = (I[IDX(x+1,y)] - 2*I[idx] + I[IDX(x-1,y)]) +
                              (I[IDX(x,y+1)] - 2*I[idx] + I[IDX(x,y-1)]);
                Rn[idx] = R[idx] + (-c * lapI + dt_half * V[idx] * I[idx]);
            }
        }
        
        // Update Imaginary Component
        for(int y=1; y<N-1; y++){
            for(int x=1; x<N-1; x++){
                int idx = IDX(x,y);
                double lapR = (Rn[IDX(x+1,y)] - 2*Rn[idx] + Rn[IDX(x-1,y)]) +
                              (Rn[IDX(x,y+1)] - 2*Rn[idx] + Rn[IDX(x,y-1)]);
                In[idx] = I[idx] - (-c * lapR + dt_half * V[idx] * Rn[idx]);
            }
        }

        // Swap Pointers
        double *tmp;
        tmp=R; R=Rn; Rn=tmp;
        tmp=I; I=In; In=tmp;

        // Output specific frames 
        if(step==5 || step==315 || step==665 || step==995){
            char filename[32];
            sprintf(filename, "psi2_step_%d.dat", step);
            FILE *fp = fopen(filename, "w");
            for(int y=0; y<N; y++){
                for(int x=0; x<N; x++){
                    int idx = IDX(x,y);
                    fprintf(fp, "%e ", R[idx]*R[idx] + I[idx]*I[idx]);
                }
                fprintf(fp, "\n");
            }
            fclose(fp);
            printf("Saved %s\n", filename);
        }
    }

    free(R); free(I); free(Rn); free(In); free(V);
    return 0;
}
