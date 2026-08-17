#include <jni.h>
#include <string.h>
#include <math.h>

#define NX 120
#define NY 120
#define NSTEPS 200

extern "C"
JNIEXPORT jstring JNICALL
Java_com_example_myfdtdapp_MainActivity_runFDTD(JNIEnv* env, jobject thiz)
{
    static float Ez[NX][NY];
    static float Hx[NX][NY];
    static float Hy[NX][NY];

    const float ce = 0.5f;
    const float chx = 0.5f;
    const float chy = 0.5f;

    memset(Ez, 0, sizeof(Ez));
    memset(Hx, 0, sizeof(Hx));
    memset(Hy, 0, sizeof(Hy));

    for (int n = 0; n < NSTEPS; n++) {

        for (int i = 0; i < NX; i++)
            for (int j = 0; j < NY - 1; j++)
                Hx[i][j] -= chx * (Ez[i][j+1] - Ez[i][j]);

        for (int i = 0; i < NX - 1; i++)
            for (int j = 0; j < NY; j++)
                Hy[i][j] += chy * (Ez[i+1][j] - Ez[i][j]);

        Ez[NX/2][NY/2] += expf(-0.5f * powf((n - 40)/10.0f, 2.0f));

        for (int i = 1; i < NX; i++)
            for (int j = 1; j < NY; j++)
                Ez[i][j] += ce * ((Hy[i][j] - Hy[i-1][j]) -
                                  (Hx[i][j] - Hx[i][j-1]));
    }

    float centre = Ez[NX/2][NY/2];
    char buf[128];
    snprintf(buf, sizeof(buf), "FDTD OK. Ez centre = %.5f", centre);

    return env->NewStringUTF(buf);
}
