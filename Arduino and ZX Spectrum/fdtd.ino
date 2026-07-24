#include <Arduino.h>
#include <stdint.h>
#include <math.h>

//--------------------------------------------------
// 1D Acoustic FDTD (Q15 Fixed Point)
//--------------------------------------------------

#define KE      200
#define STEPS   400

//--------------------------------------------------
// Q15 constants
//--------------------------------------------------

#define Q15_ONE 32767

// Courant number S = c dt / dz
// Use S = 0.5 for stability.
const int16_t S = 16384;      // 0.5 in Q15

//--------------------------------------------------
// Q15 helpers
//--------------------------------------------------

inline int16_t float_to_q15(float x)
{
    if (x >= 0.999969f) x = 0.999969f;
    if (x < -1.0f)      x = -1.0f;
    return (int16_t)(x * 32768.0f);
}

inline int16_t q15_mul(int16_t a, int16_t b)
{
    return (int16_t)(((int32_t)a * (int32_t)b) >> 15);
}

inline int16_t q15_add(int16_t a, int16_t b)
{
    int32_t t = (int32_t)a + b;

    if (t > 32767)  t = 32767;
    if (t < -32768) t = -32768;

    return (int16_t)t;
}

inline int16_t q15_sub(int16_t a, int16_t b)
{
    int32_t t = (int32_t)a - b;

    if (t > 32767)  t = 32767;
    if (t < -32768) t = -32768;

    return (int16_t)t;
}

//--------------------------------------------------
// Fields
//--------------------------------------------------

int16_t p[KE];
int16_t u[KE];

// local sound speed (normalised to cmax)
int16_t c[KE];

//--------------------------------------------------
// Mur boundaries
//--------------------------------------------------

int16_t p_left_1  = 0;
int16_t p_left_2  = 0;

int16_t p_right_1 = 0;
int16_t p_right_2 = 0;

//--------------------------------------------------
// Setup
//--------------------------------------------------

void setup()
{
    Serial.begin(115200);

    while (!Serial)
        ;

    //--------------------------------------------------
    // Initialise fields
    //--------------------------------------------------

    for (int k = 0; k < KE; k++)
    {
        p[k] = 0;
        u[k] = 0;

        // Left medium
        c[k] = float_to_q15(1.0f);
    }

    //--------------------------------------------------
    // Second material
    //--------------------------------------------------

    int interface = KE / 2;

    for (int k = interface; k < KE; k++)
    {
        // Slower medium
        c[k] = float_to_q15(1000.0f / 1500.0f);
    }

    Serial.println();
    Serial.println("--------------------------------");
    Serial.println("1D Acoustic FDTD (Q15)");
    Serial.println("--------------------------------");
    Serial.println();
}

//--------------------------------------------------
// ASCII waveform display
//--------------------------------------------------

void display_waveform(int step)
{
    Serial.print("Step ");
    Serial.println(step);

    for (int k = 0; k < KE; k++)
    {
        int16_t value = p[k];

        // Convert Q15 value to display position
        int level = value / 800;

        // Display centre line
        if (level == 0)
        {
            Serial.println("|");
            continue;
        }

        if (level > 0)
        {
            Serial.print("|");

            for (int i = 0; i < level; i++)
                Serial.print("*");

            Serial.println();
        }
        else
        {
            for (int i = 0; i < -level; i++)
                Serial.print(" ");

            Serial.println("*");
        }
    }

    Serial.println();
}

//--------------------------------------------------
// Time stepping
//--------------------------------------------------

void loop()
{
    static int n = 0;
    static int t = 0;

    if (n >= STEPS)
    {
        Serial.println("Simulation finished.");
        while (1)
            ;
    }

    t++;

    //--------------------------------------------------
    // Pressure update
    //--------------------------------------------------

    for (int k = 1; k < KE; k++)
    {
        int16_t coeff = q15_mul(S, c[k]);

        int16_t du = q15_sub(u[k - 1], u[k]);

        p[k] = q15_add(
                    p[k],
                    q15_mul(coeff, du)
                );
    }

    //--------------------------------------------------
    // Gaussian source
    //--------------------------------------------------

    {
        float t0 = 80.0f;
        float spread = 20.0f;

        float x = ((float)t - t0) / spread;

        int16_t src = float_to_q15(0.8f * expf(-x * x));

        p[15] = q15_add(p[15], src);
    }

    //--------------------------------------------------
    // Simple Mur absorbing boundaries
    //--------------------------------------------------

    p[0] = p_left_2;
    p_left_2 = p_left_1;
    p_left_1 = p[1];

    p[KE - 1] = p_right_2;
    p_right_2 = p_right_1;
    p_right_1 = p[KE - 2];

    //--------------------------------------------------
    // Velocity update
    //--------------------------------------------------

    for (int k = 0; k < KE - 1; k++)
    {
        int16_t coeff = q15_mul(S, c[k]);

        int16_t dp = q15_sub(p[k], p[k + 1]);

        u[k] = q15_add(
                    u[k],
                    q15_mul(coeff, dp)
                );
    }

    //--------------------------------------------------
    // Display every 50 time steps
    //--------------------------------------------------

    if ((n % 50) == 0)
    {
        //display_waveform(n);

        for (int k = 0; k < KE; k++)
        {
            Serial.print((float)p[k] / 32768.0f);

            if (k < KE-1)
                Serial.print(",");
        }

        Serial.println();

    }

    n++;
}

