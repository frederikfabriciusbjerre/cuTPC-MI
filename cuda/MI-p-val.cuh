#ifndef CALC_STATISTICS_CUH
#define CALC_STATISTICS_CUH

#include <cmath>
#include <cuda_runtime.h>
// #include "pt_old.cuh"  
#include "pt.cuh"

__device__ double df_reiter(double B, double U, double m, double dfcom) {
    // Intermediate calculations
    double t = m - 1.0;
    double r = (1.0 + 1.0 / m) * B / U;
    double a = r * t / (t - 2.0);
    double vstar = ((dfcom + 1.0) / (dfcom + 3.0)) * dfcom;

    double c0 = 1.0 / (t - 4.0);
    double c1 = vstar - 2.0 * (1.0 + a);
    double c2 = vstar - 4.0 * (1.0 + a);

    // Calculations for z
    double z = 1.0 / c2 +
               c0 * (pow(a, 2) * c1 / (pow(1.0 + a, 2) * c2)) +
               c0 * (8.0 * pow(a, 2) * c1 / ((1.0 + a) * pow(c2, 2)) + 4.0 * pow(a, 2) / ((1.0 + a) * c2)) +
               c0 * (4.0 * pow(a, 2) / (c2 * c1) + 16.0 * pow(a, 2) * c1 / pow(c2, 3)) +
               c0 * (8.0 * pow(a, 2) / pow(c2, 2));

    // Final calculation for v
    double v = 4.0 + 1.0 / z;

    return v;
}

__device__ double compute_MI_p_value(const double* z_m, int M, int nrows, int ord, int df_method) {
    // if m = 1
    // calculate avgz
    double z_sum = 0.0;
    for (int m = 0; m < M; m++) {
        z_sum += z_m[m];
    }
    double avgz = z_sum / M;
    // calculate within-imputation variance
    double df_com = nrows - 3 - ord;
    double W = 1.0 / df_com;

    // calculate between-imputation variance
    double B = 0.0;
    if (M > 1) {
        double B_sum = 0.0;
        for (int m = 0; m < M; m++) {
            double diff = z_m[m] - avgz;
            B_sum += diff * diff;
        }
        B = B_sum / (M - 1.0);
    }

    // calculate total variance
    double temp = (1.0 + 1.0 / M) * B;
    double TV = W + temp;

    // calculate test stat
    double ts = avgz / sqrt(TV);

    // calculate degrees of freedom
    double df;

    // calculate lambda
    double lambda = temp / TV;
    double lambda_sq = lambda * lambda;
    if (M > 1){
        // choose the degrees of freedom method
        if (df_method == 0) {
            // rubin's original approximation
            double df_old = (M - 1) * (1 + (W / B) * (M/(M + 1))) * (1 + (W / B) * (M/(M + 1)));
            df = df_old;
        } else if (df_method == 1) {
            // barnard and rubin's approximation (1999)
            // based on rewrite that can be found in master thesis by Frederik Fabricius-Bjerre
            double br_const = (1.0 - lambda) * (1.0 + df_com) * df_com;
            double df_br =  (M - 1) * br_const / ((df_com + 3.0) * (M - 1) + lambda_sq * br_const); 
            df = df_br;   
            if (isnan(df) || isinf(df)) {
                df = ((df_com + 1.0) / (df_com + 3.0)) * df_com;
            }
        } else if (df_method == 2) {
            df = df_reiter(B, W, M, df_com);
            if (isnan(df) || isinf(df)) {
                df = ((df_com + 1.0) / (df_com + 3.0)) * df_com;
            }
        } else {
            df = INFINITY; // fallback for invalid input
        }
    } else {
        df = INFINITY;
    }

    // return p-value
    // double p_val = 2.0 * (1.0 - pt(fabs(ts), df)); // if using old
    double p_val = 2.0 * pt(fabs(ts), df, 0, 0);
    return p_val;
}

#endif // CALC_STATISTICS_CUH
