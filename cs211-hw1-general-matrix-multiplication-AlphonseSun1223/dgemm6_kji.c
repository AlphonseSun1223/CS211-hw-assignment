void dgemm6_kji(double *C, double *A, double *B, int n) {
    int k, j, i;
    for (k = 0; k < n; k++) {
        for (j = 0; j < n; j++) {
            register double r = B[k * n + j];
            for (i = 0; i < n; i++) {
                C[i * n + j] += A[i * n + k] * r;
            }
        }
    }
}


void dgemm6_kji2(double *C, double *A, double *B, int n) {
    int k, jj, j, ii, i, kk;
    int b = 10;  

    for (k = 0; k < n; k += b) {
        for (j = 0; j < n; j += b) {
            for (i = 0; i < n; i += b) {
                for (ii = i; ii < i + b; ii++) {
                    for (jj = j; jj < j + b; jj++) {
                        register double r = C[ii * n + jj];
                        for (kk = k; kk < k + b; kk++) {
                            r += A[ii * n + kk] * B[kk * n + jj];
                        }
                        C[ii * n + jj] = r;
                    }
                }
            }
        }
    }
}
