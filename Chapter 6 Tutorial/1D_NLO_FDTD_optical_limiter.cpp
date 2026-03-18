#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <complex>

using namespace std;

/* ---------- types ---------- */
typedef float MyType;
typedef complex<float> DCMPLX;

/* ---------- constants ---------- */
#define Sqr(x) ((x)*(x))
#define LIGHT_SPEED 299792458.0f
#define MU_0 1.25663706e-6f
#define EPSILON_0 8.85418782e-12f
#define MAXTIMESTEP (16384*32)
#define PI 3.1415926f

/* ---------- global fields ---------- */
DCMPLX *hy,*ez,*jz,*sigma_nl_z,*sigma_l_z,*eps_rmn,*omega_rmn,*G_rmn;
DCMPLX *chi_z,*alpha,*bbeta,*c_epsz,*pz_kerr,*pz_kerr_old,*ez_old,*c_x3z;
DCMPLX **SavedField;
MyType dx;

/* ---------- allocation helpers ---------- */
DCMPLX* CArray1d(int n){
    DCMPLX* a=(DCMPLX*)malloc(n*sizeof(DCMPLX));
    for(int i=0;i<n;i++) a[i]=DCMPLX(0.0f,0.0f);
    return a;
}

DCMPLX** CArray2d(int nx,int ny){
    DCMPLX** a=(DCMPLX**)malloc(nx*sizeof(DCMPLX*));
    for(int i=0;i<nx;i++){
        a[i]=(DCMPLX*)malloc(ny*sizeof(DCMPLX));
        for(int j=0;j<ny;j++) a[i][j]=DCMPLX(0.0f,0.0f);
    }
    return a;
}

/* ---------- FDTD solver ---------- */
void calc_radiation(int timesteps,MyType dt,int nx){

    FILE *out;
    DCMPLX eps_eff,c1,c2;
    DCMPLX da(dt/dx,0.0f);

    DCMPLX ez_low_m1=0.0f,ez_low_m2=0.0f;
    DCMPLX ez_high_m1=0.0f,ez_high_m2=0.0f;

    MyType tw=120*dt;
    MyType t0=2.1f*tw;

    for(int i=0;i<timesteps;i++){

        /* debug output */
        if(i%100==0){
            printf("Timestep %d\n",i);
            char filename[64];
            sprintf(filename,"step%d.dat",i);
            out=fopen(filename,"w");
            for(int l=0;l<=nx;l++) fprintf(out,"%g\n",real(ez[l]));
            fclose(out);
        }

        /* update H field */
        for(int l=1;l<=nx;l++)
            hy[l]+= (da/MU_0)*(ez[l]-ez[l-1]);

        /* update E field */
        for(int l=0;l<=nx;l++){
            eps_eff=EPSILON_0*(c_epsz[l]+chi_z[l]);

            c1=1.0f-(dt*sigma_l_z[l])/(2.0f*eps_eff);
            c2=1.0f+(dt*sigma_l_z[l])/(2.0f*eps_eff);

            ez_old[l]=ez[l];

            ez[l]=(c1/c2)*ez[l]
                 +(da/(eps_eff*c2))*(hy[l+1]-hy[l])
                 -(dt/(eps_eff*c2))*jz[l]
                 -(pz_kerr[l]-pz_kerr_old[l])/(eps_eff*c2);
        }

        /* excitation pulse */
        //ez[3]+=1e9f*exp(-pow((i*dt-t0)/tw,2.0f))*sin(2.0f*PI*6.4620e14*(i*dt-t0));
        ez[3]+=4e9f*exp(-pow((i*dt-t0)/tw,2.0f))*sin(2.0f*PI*2.8837e14f*(i*dt-t0));
        //ez[3]+=1e0f*sin(2.0f*PI*6.4620e14*i*dt);
        //broadband delta pulse
        //if(i==0) ez[3] = 3.0e10f; else ez[3] = ez[3];

        /* absorbing boundary */
        ez[0]=ez_low_m2;
        ez_low_m2=ez_low_m1;
        ez_low_m1=ez[1];

        ez[nx]=ez_high_m2;
        ez_high_m2=ez_high_m1;
        ez_high_m1=ez[nx-1];

        /* nonlinear material update */
        for(int l=0;l<=nx;l++){

            jz[l]=((1.0f+alpha[l]*dt/2.0f)/(1.0f-alpha[l]*dt/2.0f))*jz[l]
                 +(bbeta[l]*dt/(1.0f-alpha[l]*dt/2.0f))*ez[l];

            sigma_nl_z[l]=((1.0f-G_rmn[l]*dt/2.0f)/(1.0f+G_rmn[l]*dt/2.0f))*sigma_nl_z[l]
                         +(Sqr(omega_rmn[l])*dt/(1.0f+G_rmn[l]*dt/2.0f))
                         *(eps_rmn[l]*Sqr(abs(ez[l]))-chi_z[l]);

            chi_z[l]+=dt*sigma_nl_z[l];

            pz_kerr_old[l]=pz_kerr[l];
            pz_kerr[l]=EPSILON_0*c_x3z[l]*Sqr(ez[l])*ez[l];
        }

        /* probes */
        SavedField[0][i]=ez[150];
        SavedField[1][i]=ez[850];
        SavedField[2][i]=hy[150];
        SavedField[3][i]=hy[850];
    }

    /* save output */
    out=fopen("spectrum.dat","w");
    for(int i=0;i<MAXTIMESTEP;i++)
        fprintf(out,"%g %g\n",SavedField[1][i].real(),SavedField[3][i].real());
    fclose(out);
}

/* ---------- material profile ---------- */
void calc_epsilon(int nx){

    FILE *out=fopen("Eps.dat","w");
    FILE *out1=fopen("bbeta.dat","w");

    for(int i=0;i<=nx;i++){
        c_epsz[i]=1.0f;
        c_x3z[i]=alpha[i]=bbeta[i]=0.0f;
        eps_rmn[i]=omega_rmn[i]=G_rmn[i]=sigma_l_z[i]=0.0f;
    }

	int j_center = (160 + 840)/2;

	for(int i=0;i<=nx;i++){
	    for(int j=160;j<=840;j+=20){
		if(j==j_center) continue;//comment to remove the 'defect'

			if(i>j && i<=j+10){
				c_epsz[i]=3.4f*3.4f;
				c_x3z[i]=4.4e-19f;
			}
	    }

        fprintf(out,"%g\n",real(c_epsz[i]));
        fprintf(out1,"%g\n",real(bbeta[i]));
    }

    fclose(out);
    fclose(out1);
}

/* ---------- main ---------- */
int main(){

    int nx=1000;

    dx=20e-9f;
    MyType dt=dx/(2.0f*LIGHT_SPEED);

    printf("dt=%e\n",dt);

    /* allocate memory */
    ez=CArray1d(nx+1);
    ez_old=CArray1d(nx+1);
    pz_kerr=CArray1d(nx+1);
    pz_kerr_old=CArray1d(nx+1);
    jz=CArray1d(nx+1);

    sigma_l_z=CArray1d(nx+1);
    sigma_nl_z=CArray1d(nx+1);
    omega_rmn=CArray1d(nx+1);
    eps_rmn=CArray1d(nx+1);
    G_rmn=CArray1d(nx+1);
    chi_z=CArray1d(nx+1);

    hy=CArray1d(nx+2);

    c_epsz=CArray1d(nx+1);
    c_x3z=CArray1d(nx+1);
    alpha=CArray1d(nx+1);
    bbeta=CArray1d(nx+1);

    SavedField=CArray2d(4,MAXTIMESTEP);

    /* run simulation */
    calc_epsilon(nx);
    calc_radiation(MAXTIMESTEP,dt,nx);

    return 0;
}
