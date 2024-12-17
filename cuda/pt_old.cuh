#ifndef PT_DEVICE_CUH
#define PT_DEVICE_CUH

#include <cuda_runtime.h>
#include <math_constants.h>


// Device function to compute the continued fraction for the incomplete beta function
__device__ double betacf(double a, double b, double x) {
    const int MAX_ITER = 1000;
    const double EPS = 1.0e-15;
    const double FPMIN = 1.0e-30;

    double qab = a + b;
    double qap = a + 1.0;
    double qam = a - 1.0;

    double c = 1.0;
    double d = 1.0 - qab * x / qap;
    if (fabs(d) < FPMIN) d = FPMIN;
    d = 1.0 / d;
    double h = d;

    for (int m = 1; m <= MAX_ITER; m++) {
        int m2 = 2 * m;
        // Even step
        double aa = m * (b - m) * x / ((qam + m2) * (a + m2));
        d = 1.0 + aa * d;
        if (fabs(d) < FPMIN) d = FPMIN;
        c = 1.0 + aa / c;
        if (fabs(c) < FPMIN) c = FPMIN;
        d = 1.0 / d;
        h *= d * c;
        // Odd step
        aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2));
        d = 1.0 + aa * d;
        if (fabs(d) < FPMIN) d = FPMIN;
        c = 1.0 + aa / c;
        if (fabs(c) < FPMIN) c = FPMIN;
        d = 1.0 / d;
        double del = d * c;
        h *= del;
        if (fabs(del - 1.0) < EPS) break;
    }

    return h;
}

// Device function to compute the incomplete beta function
__device__ double betai(double a, double b, double x) {
    if (x < 0.0 || x > 1.0) return 0.0; // Invalid input
    if (x == 0.0 || x == 1.0) return x; // Edge cases

    // Compute ln(Beta(a, b))
    double ln_beta = lgamma(a) + lgamma(b) - lgamma(a + b);

    // Compute front factor
    double front = exp(log(x) * a + log(1.0 - x) * b - ln_beta) / a;

    // Compute continued fraction
    double cf = betacf(a, b, x);

    return front * cf;
}

// Device function to compute the CDF of the Student's t-distribution
__device__ double pt(double t, double df) {
    if (isinf(df)) {
        // Use the standard normal distribution CDF for infinite degrees of freedom
        return 0.5 * (1.0 + erf(t / sqrt(2.0)));  
    }
    double x = df / (df + t * t);
    double a = df / 2.0;
    double b = 0.5;

    double ibeta = betai(a, b, x);

    double cdf;
    if (t >= 0) {
        cdf = 1.0 - 0.5 * ibeta;
    } else {
        cdf = 0.5 * ibeta;
    }

    return cdf;
}

// float versions
__device__ float betacff(float a, float b, float x) {
    const int MAX_ITER = 1000000;
    const float EPS = 1.0e-10;
    const float FPMIN = 1.0e-30;

    float qab = a + b;
    float qap = a + 1.0;
    float qam = a - 1.0;

    float c = 1.0;
    float d = 1.0 - qab * x / qap;
    if (fabs(d) < FPMIN) d = FPMIN;
    d = 1.0 / d;
    float h = d;

    for (int m = 1; m <= MAX_ITER; m++) {
        int m2 = 2 * m;
        // Even step
        float aa = m * (b - m) * x / ((qam + m2) * (a + m2));
        d = 1.0 + aa * d;
        if (fabs(d) < FPMIN) d = FPMIN;
        c = 1.0 + aa / c;
        if (fabs(c) < FPMIN) c = FPMIN;
        d = 1.0 / d;
        h *= d * c;
        // Odd step
        aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2));
        d = 1.0 + aa * d;
        if (fabs(d) < FPMIN) d = FPMIN;
        c = 1.0 + aa / c;
        if (fabs(c) < FPMIN) c = FPMIN;
        d = 1.0 / d;
        float del = d * c;
        h *= del;
        if (fabs(del - 1.0) < EPS) break;
    }

    return h;
}

// Device function to compute the incomplete beta function
__device__ float betaif(float a, float b, float x) {
    if (x < 0.0 || x > 1.0) return 0.0; // Invalid input
    if (x == 0.0 || x == 1.0) return x; // Edge cases

    // Compute ln(Beta(a, b))
    float ln_beta = lgammaf(a) + lgammaf(b) - lgammaf(a + b);

    // Compute front factor
    float front = exp(log(x) * a + log(1.0 - x) * b - ln_beta) / a;

    // Compute continued fraction
    float cf = betacff(a, b, x);

    return front * cf;
}

// Device function to compute the CDF of the Student's t-distribution
__device__ float ptf(double t, double df) {
    // Convert input doubles to floats for internal calculations
    float tf = static_cast<float>(t);
    float dff = static_cast<float>(df);

    if (isinf(dff)) {
        // Use the standard normal distribution CDF for infinite degrees of freedom
        return 0.5f * (1.0f + erf(tf / sqrt(2.0f)));  
    }
    float x = dff / (dff + tf * tf);
    float a = dff / 2.0f;
    float b = 0.5f;

    float ibeta = betaif(a, b, x);

    float cdf;
    if (tf >= 0) {
        cdf = 1.0f - 0.5f * ibeta;
    } else {
        cdf = 0.5f * ibeta;
    }

    return static_cast<double>(cdf);
}

#endif // PT_DEVICE_CUH
