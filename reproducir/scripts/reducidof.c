/* Versión rápida de campo() en precisión simple (se compila con -Ofast: raíz inversa aproximada
   con una iteración de Newton). Las sumas se acumulan en simple dentro de cada fila. */
#include <math.h>

void campof(int n, const float *x, const float *y, const float *m, float *sx, float *sy)
{
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < n; i++) {
        float ax = 0.0f, ay = 0.0f;
        const float xi = x[i], yi = y[i];
        #pragma omp simd reduction(+:ax,ay)
        for (int j = 0; j < n; j++) {
            float dx = xi - x[j], dy = yi - y[j];
            float inv = m[j]/sqrtf(dx*dx + dy*dy + 1e-30f);
            ax += dx*inv;
            ay += dy*inv;
        }
        sx[i] = ax;
        sy[i] = ay;
    }
}

/* Momento angular distribuido: u_i = sqrt(kappa_0) R(L_i) es el radio de la órbita circular de la
   partícula i en unidades del plano de fase, y d = u_i - u_j. El núcleo de interacción promediado es
   G(rho, d) = sqrt(rho^2 - d^2) + d asin(d/rho) si |d| < rho, (pi/2)|d| si no, y
   s_i = sum_j m_j sqrt(max(rho^2 - d^2, 0))/rho^2 (Z_i - Z_j). */
void campofl(int n, const float *x, const float *y, const float *u, const float *m, float *sx, float *sy)
{
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < n; i++) {
        float ax = 0.0f, ay = 0.0f;
        const float xi = x[i], yi = y[i], ui = u[i];
        #pragma omp simd reduction(+:ax,ay)
        for (int j = 0; j < n; j++) {
            float dx = xi - x[j], dy = yi - y[j], du = ui - u[j];
            float r2 = dx*dx + dy*dy;
            float q = fmaxf(r2 - du*du, 0.0f);
            float w = m[j]*sqrtf(q)/(r2 + 1e-30f);
            ax += dx*w;
            ay += dy*w;
        }
        sx[i] = ax;
        sy[i] = ay;
    }
}
