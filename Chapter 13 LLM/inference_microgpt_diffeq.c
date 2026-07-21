/*
 * GPT with Subword Tokenization - Differential Equation Solver Version (Inference Only)
 * Uses RPN (Reverse Polish Notation) for mathematical expressions
 * * Compilation:
  gcc -Ofast -mcmodel=medium -march=native -ffast-math -funroll-loops -flto \
  -I./llama.cpp/include -I./llama.cpp/ggml/include \
  -L./llama.cpp/build/bin -Wl,-rpath,./llama.cpp/build/bin \
  inference_microgpt_diffeq_V6.c -lggml-base -lggml -lllama -lopenblas -lm -lpthread -o aef
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <ctype.h>
#include <cblas.h>
#include "gguf.h"
#include "ggml.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* Configuration */
#define MAX_VOCAB 512
#define MAX_TOKEN_LEN 64
#define MAX_DOC_LEN 512

/* Model hyper-parameters */
#define N_EMBD 512
#define N_HEAD 8
#define N_LAYER 12
#define BLOCK_SIZE 128
#define HEAD_DIM (N_EMBD / N_HEAD)
#define MLP_DIM (4 * N_EMBD)
#define BATCH_SIZE 16

/* Global Tokenizer Variables */
static char vocab_tokens[MAX_VOCAB][MAX_TOKEN_LEN];
static int vocab_count = 0;
static int vocab_size, BOS_id, EOS_id, UNK_id;
static int space_id = -1;
static int sol_id = -1;

/* Active Inference Parameters */
static float *wte;
static float *wpe;
static float *lm_head;

static float *attn_wq[N_LAYER];
static float *attn_wk[N_LAYER];
static float *attn_wv[N_LAYER];
static float *attn_wo[N_LAYER];
static float *mlp_fc1[N_LAYER];
static float *mlp_fc2[N_LAYER];

/* Forward Activations Struct */
typedef struct {
    float x_embed[N_EMBD];
    float rms_scale_init;
    float x_in[N_LAYER][N_EMBD];
    float xn_attn[N_LAYER][N_EMBD];
    float rms_scale_attn[N_LAYER];
    float q[N_LAYER][N_EMBD];
    float aw[N_LAYER][N_HEAD][BLOCK_SIZE];
    float attn_out[N_LAYER][N_EMBD];
    float x_mid[N_LAYER][N_EMBD];
    float xn_mlp[N_LAYER][N_EMBD];
    float rms_scale_mlp[N_LAYER];
    float mlp_pre[N_LAYER][MLP_DIM];
    float mlp_post[N_LAYER][MLP_DIM];
    float x_out[N_EMBD];
} PosActs;

/* Static KV Cache Context */
static float kv_keys[BATCH_SIZE][N_LAYER][BLOCK_SIZE][N_EMBD];
static float kv_vals[BATCH_SIZE][N_LAYER][BLOCK_SIZE][N_EMBD];

/* GGUF Global Context */
static struct gguf_context *g_gguf = NULL;

/* ------------------------------------------------------------------ */
/* RNG Utilities (Inference Sampling)                                 */
/* ------------------------------------------------------------------ */
static unsigned long long rng_state = 42;

static unsigned long long rng_next(void) {
    rng_state ^= rng_state << 13;
    rng_state ^= rng_state >> 7;
    rng_state ^= rng_state << 17;
    return rng_state;
}

static double rng_uniform(void) {
    return (rng_next() >> 11) * (1.0 / 9007199254740992.0);
}

/* ------------------------------------------------------------------ */
/* Tokenizer Functions                                                */
/* ------------------------------------------------------------------ */
static int find_token_id(const char *token_str) {
    for (int i = 0; i < vocab_count; i++)
        if (strcmp(vocab_tokens[i], token_str) == 0) return i;
    return -1;
}

static void encode_de_string(const char *input, int *output_tokens, int *count) {
    int n = 0;
    output_tokens[n++] = BOS_id;

    char buffer[MAX_DOC_LEN];
    strcpy(buffer, input);
    char *token = strtok(buffer, " \n");

    while (token && n < MAX_DOC_LEN - 1) {
        int id = find_token_id(token);
        output_tokens[n++] = (id != -1) ? id : UNK_id;
        token = strtok(NULL, " \n");
    }
    *count = n;
}

static const char* id_to_token(int id) {
    if (id < 0 || id >= vocab_count) return "[INVALID]";
    if (id == BOS_id) return "";
    if (id == EOS_id) return "\n";
    if (id == UNK_id) return "";
    if (id == sol_id) return "[SOL]";
    return vocab_tokens[id];
}

static void build_tokenizer(void) {
    FILE *f = fopen("math_vocab.bin", "rb");
    if (!f) {
        printf("Error: math_vocab.bin not found\n");
        exit(1);
    }

    vocab_count = 0;
    char buffer[MAX_TOKEN_LEN];
    int idx = 0;
    int c;

    while ((c = fgetc(f)) != EOF) {
        if (c == '\0' || c == '\n' || c == ' ') {
            if (idx > 0) {
                buffer[idx] = '\0';
                if (vocab_count >= MAX_VOCAB) break;
                strncpy(vocab_tokens[vocab_count], buffer, MAX_TOKEN_LEN - 1);
                vocab_tokens[vocab_count][MAX_TOKEN_LEN - 1] = '\0';
                vocab_count++;
                idx = 0;
            }
        } else {
            if (idx < MAX_TOKEN_LEN - 1) buffer[idx++] = (char)c;
        }
    }
    fclose(f);

    UNK_id = 3;
    BOS_id = 1;
    EOS_id = 2;
    vocab_size = vocab_count;
    space_id = find_token_id(" ");
    sol_id   = find_token_id("[SOL]");

    printf("Loaded vocab=%d space=%d sol=%d\n", vocab_count, space_id, sol_id);
}

/* ------------------------------------------------------------------ */
/* Mathematical Processing Layers                                     */
/* ------------------------------------------------------------------ */
static inline void linear_fwd(const float *restrict x, const float *restrict w, 
                             int nout, int nin, float *restrict out) {
    cblas_sgemv(CblasRowMajor, CblasNoTrans, nout, nin, 1.0f, w, nin, x, 1, 0.0f, out, 1);
}

static inline float rmsnorm_fwd(const float *x, int n, float *out) {
    float ms = 0;
    for (int i = 0; i < n; i++) ms += x[i] * x[i];
    ms /= n;
    float scale = 1.0f / sqrtf(ms + 1e-5f);
    for (int i = 0; i < n; i++) out[i] = x[i] * scale;
    return scale;
}

static inline void softmax_fwd(const float *logits, int n, float *probs) {
    float mx = logits[0];
    for (int i = 1; i < n; i++) if (logits[i] > mx) mx = logits[i];
    
    float sum = 0;
    for (int i = 0; i < n; i++) {
        probs[i] = expf(logits[i] - mx);
        sum += probs[i];
    }
    float inv = 1.0f / sum;
    for (int i = 0; i < n; i++) probs[i] *= inv;
}

/* ------------------------------------------------------------------ */
/* Core Network Evaluation Pipeline                                   */
/* ------------------------------------------------------------------ */
static void gpt_forward(int batch_idx, int token_id, int pos_id, float *logits_out, PosActs *act) {
    float x[N_EMBD], tmp[MLP_DIM > N_EMBD ? MLP_DIM : N_EMBD];

    for (int i = 0; i < N_EMBD; i++)
        x[i] = wte[token_id * N_EMBD + i] + wpe[pos_id * N_EMBD + i];
    memcpy(act->x_embed, x, sizeof(x));

    act->rms_scale_init = rmsnorm_fwd(x, N_EMBD, x);

    for (int li = 0; li < N_LAYER; li++) {
        memcpy(act->x_in[li], x, sizeof(x));

        float xn[N_EMBD];
        act->rms_scale_attn[li] = rmsnorm_fwd(x, N_EMBD, xn);
        memcpy(act->xn_attn[li], xn, sizeof(xn));

        float q[N_EMBD], k[N_EMBD], v[N_EMBD];
        linear_fwd(xn, attn_wq[li], N_EMBD, N_EMBD, q);
        linear_fwd(xn, attn_wk[li], N_EMBD, N_EMBD, k);
        linear_fwd(xn, attn_wv[li], N_EMBD, N_EMBD, v);
        memcpy(act->q[li], q, sizeof(q));

        memcpy(kv_keys[batch_idx][li][pos_id], k, sizeof(k));
        memcpy(kv_vals[batch_idx][li][pos_id], v, sizeof(v));
        
        int seq_len = pos_id + 1;
        float scale = 1.0f / sqrtf((float)N_EMBD / (float)N_HEAD);
        float ao[N_EMBD];

        for (int h = 0; h < N_HEAD; h++) {
            int hs = h * HEAD_DIM;
            float al[BLOCK_SIZE];
            
            for (int tt = 0; tt < seq_len; tt++) {
                float dot = 0;
                for (int j = 0; j < HEAD_DIM; j++)
                    dot += q[hs + j] * kv_keys[batch_idx][li][tt][hs + j];
                al[tt] = dot * scale;
            }
            
            float mx = al[0];
            for (int tt = 1; tt < seq_len; tt++)
                if (al[tt] > mx) mx = al[tt];
            
            float sm = 0;
            for (int tt = 0; tt < seq_len; tt++) {
                al[tt] = expf(al[tt] - mx);
                sm += al[tt];
            }
            
            float inv = 1.0f / sm;
            for (int tt = 0; tt < seq_len; tt++) {
                al[tt] *= inv;
                act->aw[li][h][tt] = al[tt];
            }
            
            for (int j = 0; j < HEAD_DIM; j++) {
                float s = 0;
                for (int tt = 0; tt < seq_len; tt++)
                    s += al[tt] * kv_vals[batch_idx][li][tt][hs + j];
                ao[hs + j] = s;
            }
        }
        
        memcpy(act->attn_out[li], ao, sizeof(ao));
        linear_fwd(ao, attn_wo[li], N_EMBD, N_EMBD, tmp);
        
        for (int i = 0; i < N_EMBD; i++)
            x[i] = tmp[i] + act->x_in[li][i];
        memcpy(act->x_mid[li], x, sizeof(x));

        float xn_m[N_EMBD];
        act->rms_scale_mlp[li] = rmsnorm_fwd(x, N_EMBD, xn_m);
        memcpy(act->xn_mlp[li], xn_m, sizeof(xn_m));

        float h1[MLP_DIM];
        linear_fwd(xn_m, mlp_fc1[li], MLP_DIM, N_EMBD, h1);
        memcpy(act->mlp_pre[li], h1, MLP_DIM * sizeof(float));

        float h2[MLP_DIM];
        for (int i = 0; i < MLP_DIM; i++)
            h2[i] = h1[i] > 0 ? h1[i] * h1[i] : 0;
        memcpy(act->mlp_post[li], h2, MLP_DIM * sizeof(float));

        linear_fwd(h2, mlp_fc2[li], N_EMBD, MLP_DIM, tmp);
        for (int i = 0; i < N_EMBD; i++)
            x[i] = tmp[i] + act->x_mid[li][i];
    }

    memcpy(act->x_out, x, sizeof(x));
    linear_fwd(x, lm_head, vocab_size, N_EMBD, logits_out);
}

static int weighted_choice(const float *w, int n) {
    float total = 0;
    for (int i = 0; i < n; i++) total += w[i];
    float r = (float)rng_uniform() * total, cum = 0;
    for (int i = 0; i < n; i++) {
        cum += w[i];
        if (r < cum) return i;
    }
    return n - 1;
}

/* ------------------------------------------------------------------ */
/* Autoregressive Generative Execution Token Loop                     */
/* ------------------------------------------------------------------ */
static void solve_de(const char *de_problem_rpn) {
    printf("\n--- Solving DE: %s ---\n", de_problem_rpn);
    
    int problem_tokens[BLOCK_SIZE];
    int problem_len = 0;
    encode_de_string(de_problem_rpn, problem_tokens, &problem_len);
    
    memset(kv_keys, 0, sizeof(kv_keys));
    memset(kv_vals, 0, sizeof(kv_vals));
    
    float logits_solve[MAX_VOCAB];
    
    for (int pos = 0; pos < problem_len; pos++) {
        PosActs act;
        gpt_forward(0, problem_tokens[pos], pos, logits_solve, &act);
    }
    
    printf("Generated solution (RPN): ");
    int current_token = sol_id;  
    int pos = problem_len;
    
    if (sol_id != -1) {
        printf("[SOL] ");
        fflush(stdout);
    }
    
    int eq_id = find_token_id("=");
    int y_id = find_token_id("y");
    int dy_id = find_token_id("y'");
    int ddy_id = find_token_id("y''");
    
    for (int gen_pos = 0; gen_pos < 100; gen_pos++) {  
        PosActs act;
        gpt_forward(0, current_token, pos, logits_solve, &act);
        
        logits_solve[eq_id] = -1e10f;
        logits_solve[y_id] = -1e10f;
        logits_solve[dy_id] = -1e10f;
        logits_solve[ddy_id] = -1e10f;
        
        for (int i = 0; i < vocab_size; i++)
            logits_solve[i] /= 0.8f;
        
        float probs_solve[MAX_VOCAB];
        softmax_fwd(logits_solve, vocab_size, probs_solve);
        current_token = weighted_choice(probs_solve, vocab_size);
        
        if (current_token == EOS_id || current_token == BOS_id) break;
        
        const char *token_str = id_to_token(current_token);
        if (token_str && token_str[0] != '\0' && token_str[0] != '[') {
            printf("%s ", token_str);
            fflush(stdout);
        } else if (token_str && strcmp(token_str, "\n") != 0 && token_str[0] != '\0') {
            if (strcmp(token_str, "[SOL]") != 0) {
                printf("%s ", token_str);
                fflush(stdout);
            }
        }
        
        pos++;
        if (pos >= BLOCK_SIZE - 1) break;
    }
    printf("\n\n");
}

/* ------------------------------------------------------------------ */
/* WEIGHT ENGINE: REFACTORED GGUF LOADER                              */
/* ------------------------------------------------------------------ */
static void copy_tensor_f32(struct ggml_tensor *t, float *dst) {
    memcpy(dst, t->data, ggml_nbytes(t));
}

static void load_weights_gguf(const char *fname) {
    printf("Loading weights from GGUF...\n");
    struct ggml_context *ggml_ctx = NULL;
    struct gguf_init_params params = { .no_alloc = false, .ctx = &ggml_ctx };

    g_gguf = gguf_init_from_file(fname, params);
    if (!g_gguf || !ggml_ctx) { exit(1); }

    struct ggml_tensor *t;

    t = ggml_get_tensor(ggml_ctx, "tok_embd.weight");
    if (!t) { exit(1); }
    wte = malloc(ggml_nbytes(t));
    copy_tensor_f32(t, wte);

    t = ggml_get_tensor(ggml_ctx, "pos_embd.weight");
    if (t) {
        wpe = malloc(ggml_nbytes(t));
        copy_tensor_f32(t, wpe);
    } else {
        wpe = calloc(BLOCK_SIZE * N_EMBD, sizeof(float));
    }

    t = ggml_get_tensor(ggml_ctx, "lm_head.weight");
    if (!t) {
        lm_head = wte;
    } else {
        lm_head = malloc(ggml_nbytes(t));
        copy_tensor_f32(t, lm_head);
    }

    int as = N_EMBD * N_EMBD;
    int ms = MLP_DIM * N_EMBD; 
    char name[256];

    for (int i = 0; i < N_LAYER; i++) {
        snprintf(name, sizeof(name), "blk.%d.attn_q.weight", i);
        t = ggml_get_tensor(ggml_ctx, name);
        attn_wq[i] = malloc(as * sizeof(float));
        if (t) copy_tensor_f32(t, attn_wq[i]); else memset(attn_wq[i], 0, as * sizeof(float));

        snprintf(name, sizeof(name), "blk.%d.attn_k.weight", i);
        t = ggml_get_tensor(ggml_ctx, name);
        attn_wk[i] = malloc(as * sizeof(float));
        if (t) copy_tensor_f32(t, attn_wk[i]); else memset(attn_wk[i], 0, as * sizeof(float));

        snprintf(name, sizeof(name), "blk.%d.attn_v.weight", i);
        t = ggml_get_tensor(ggml_ctx, name);
        attn_wv[i] = malloc(as * sizeof(float));
        if (t) copy_tensor_f32(t, attn_wv[i]); else memset(attn_wv[i], 0, as * sizeof(float));

        snprintf(name, sizeof(name), "blk.%d.attn_wo.weight", i);
        t = ggml_get_tensor(ggml_ctx, name);
        attn_wo[i] = malloc(as * sizeof(float));
        if (t) copy_tensor_f32(t, attn_wo[i]); else memset(attn_wo[i], 0, as * sizeof(float));

        snprintf(name, sizeof(name), "blk.%d.ffn.fc1.weight", i);
        t = ggml_get_tensor(ggml_ctx, name);
        mlp_fc1[i] = malloc(ms * sizeof(float));
        if (t) copy_tensor_f32(t, mlp_fc1[i]); else memset(mlp_fc1[i], 0, ms * sizeof(float));

        snprintf(name, sizeof(name), "blk.%d.ffn.fc2.weight", i);
        t = ggml_get_tensor(ggml_ctx, name);
        mlp_fc2[i] = malloc(ms * sizeof(float));
        if (t) copy_tensor_f32(t, mlp_fc2[i]); else memset(mlp_fc2[i], 0, ms * sizeof(float));
    }
    printf("Loaded all %d layers cleanly from GGUF.\n", N_LAYER);
}

/* ------------------------------------------------------------------ */
/* Execution Target Entrypoint                                        */
/* ------------------------------------------------------------------ */
int main(void) {
    openblas_set_num_threads(1);
    printf("DE Inference Engine Initialized.\n");

    build_tokenizer();
    
    /* Interchangeable Weight Load Engine Directives */
    // load_weights();
    load_weights_gguf("vector_model_complete.gguf");
    
	printf("\n[Standard Tests]\n");
	solve_de("y' 1 y * + 0 =");      // Expected: [SOL] C 1 NEG x * exp *
	solve_de("y' 2 y * + 0 =");      // Expected: [SOL] C 2 NEG x * exp *
	solve_de("y'' 1 y * + 0 =");     // Expected: [SOL] C 1 x * cos * D 1 x * sin * +

	printf("\n[Algebraic Logic - Sign Flipping]\n");
	// These test if the Sign Penalty is working
	solve_de("y' 1 y * - 0 =");      // Expected: [SOL] C 1 POS x * exp *
	solve_de("y' 3 NEG y * + 0 =");  // Expected: [SOL] C 3 POS x * exp *

	printf("\n[Coefficient Stress Test - Primes & Magnitude]\n");
	solve_de("y' 7 y * + 0 =");      // Expected: [SOL] C 7 NEG x * exp *
	solve_de("y' 10 y * + 0 =");     // Expected: [SOL] C 10 NEG x * exp *

	printf("\n[Frequency Logic - Harmonic Oscillator]\n");
	solve_de("y'' 4 y * + 0 =");     // Checks if sqrt(4)=2: [SOL] C 2 x * cos * D 2 x * sin * +
	solve_de("y'' 9 y * + 0 =");     // Checks if sqrt(9)=3: [SOL] C 3 x * cos * D 3 x * sin * +

	printf("\n[Order & Structure Tests]\n");
	solve_de("1 y * y' + 0 =");      // Commutative check: [SOL] C 1 NEG x * exp *
	solve_de("y' y 1 * + 0 =");      // Variable placement check

	printf("\n[Decay vs Oscillation - Structural Switch]\n");
	solve_de("y'' 1 y * + 0 ="); // Should stay Trig: C 1 x * cos * D 1 x * sin * +
	solve_de("y'' 1 y * - 0 ="); // Needs to switch to Exp: C 1 POS x * exp * D 1 NEG x * exp * +

	printf("\n[Advanced Phase/Frequency Stress]\n");
	solve_de("y'' 16 y * + 0 =");   // Tests the limit of square root logic (sqrt 16 = 4)            
	solve_de("y'' 25 y * + 0 =");   // Tests the limit of square root logic (sqrt 25 = 5)            

	printf("\n[Advanced - Second Order Decay/Growth]\n");
	// Tests if the model recognizes the difference between y'' + ky and y'' - ky
	solve_de("y'' 1 y * - 0 =");     // Expected: [SOL] C 1 POS x * exp * D 1 NEG x * exp * + 
	solve_de("y'' 4 y * - 0 =");     // Expected: [SOL] C 2 POS x * exp * D 2 NEG x * exp * +

	printf("\n[Structural Logic & Noise Robustness]\n");            
	// 1. Implicit Identity: y' + y = 0
	// Logic: coeff is 1, so exponent is 1 NEG.
	solve_de("y' y + 0 =");  
	// Expected: [SOL] C 1 NEG x * exp *

	// 2. Algebraic Distraction: y' + 2y = 5 - 5 (which is y' + 2y = 0)
	// Logic: coeff is 2, so exponent is 2 NEG.
	solve_de("y' 2 y * + 5 5 - ="); 
	// Expected: [SOL] C 2 NEG x * exp *

	// 3. Commutative Coefficients: y' + (y * -3) = 0
	// Note: input changed to use the new layout "3 NEG" instead of string "-3"
	solve_de("y' y 3 NEG * + 0 =");
	// Expected: [SOL] C 3 POS x * exp *

	// 4. Inverted Physics: 0 = 4y + y''
	// Logic: y'' + 4y = 0. sqrt(4) = 2. Trig solution.
	solve_de("0 y 4 * + y'' + =");
	// Expected: [SOL] C 2 x * cos * D 2 x * sin * +

    return 0;
}





