/* Modelo reducido del bulto: suma directa sobre pares en el plano de fase del epiciclo.
   campo:   s_i = sum_j m_j (Z_i - Z_j)/|Z_i - Z_j|     (dos componentes)
   energia: u_i = sum_j m_j |Z_i - Z_j|                                              */
#include <math.h>

void campo(int n, const double *x, const double *y, const double *m, double *sx, double *sy)
{
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < n; i++) {
        double ax = 0.0, ay = 0.0;
        const double xi = x[i], yi = y[i];
        #pragma omp simd reduction(+:ax,ay)
        for (int j = 0; j < n; j++) {
            double dx = xi - x[j], dy = yi - y[j];
            double inv = m[j]/sqrt(dx*dx + dy*dy + 1e-300);
            ax += dx*inv;
            ay += dy*inv;
        }
        sx[i] = ax;
        sy[i] = ay;
    }
}

void energia(int n, const double *x, const double *y, const double *m, double *u)
{
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < n; i++) {
        double a = 0.0;
        const double xi = x[i], yi = y[i];
        #pragma omp simd reduction(+:a)
        for (int j = 0; j < n; j++) {
            double dx = xi - x[j], dy = yi - y[j];
            a += m[j]*sqrt(dx*dx + dy*dy);
        }
        u[i] = a;
    }
}

/* Energía con momento angular distribuido: e_i = sum_j m_j G(rho_ij, d_ij), d_ij = u_i - u_j. */
void energial(int n, const double *x, const double *y, const double *u, const double *m, double *e)
{
    const double pi2 = 1.5707963267948966;
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < n; i++) {
        double a = 0.0;
        const double xi = x[i], yi = y[i], ui = u[i];
        for (int j = 0; j < n; j++) {
            if (j == i) continue;
            double dx = xi - x[j], dy = yi - y[j], d = ui - u[j];
            double r = sqrt(dx*dx + dy*dy);
            a += m[j]*(fabs(d) < r ? sqrt(r*r - d*d) + d*asin(d/r) : pi2*fabs(d));
        }
        e[i] = a;
    }
}
