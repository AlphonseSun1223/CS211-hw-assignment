void dgemm3(double *C, double *A, double *B, int n) {
    int i, j, k;
    
    for (i = 0; i < n; i += 3) {
        for (j = 0; j < n; j += 3) {
            
            register int t = i * n + j, t1 = t + n, t2 = t1 + n;
            register double c00 = C[t], c01 = C[t + 1], c02 = C[t + 2];
            register double c10 = C[t1], c11 = C[t1 + 1], c12 = C[t1 + 2];
            register double c20 = C[t2], c21 = C[t2 + 1], c22 = C[t2 + 2];

            for (k = 0; k < n; k += 3) {
                
                register int ta = i * n + k, ta1 = ta + n, ta2 = ta1 + n;
                register int tb = k * n + j, tb1 = tb + n, tb2 = tb1 + n;

                
                register double a00 = A[ta], a01 = A[ta1], a10 = A[ta2];
                
                
                register double b00 = B[tb], b01 = B[tb + 1], b10 = B[tb + 2];

                
                c00 += a00 * b00; c01 += a00 * b01; c02 += a00 * b10;
                c10 += a01 * b00; c11 += a01 * b01; c12 += a01 * b10;
                c20 += a10 * b00; c21 += a10 * b01; c22 += a10 * b10;

                
                a00 = A[ta + 1]; a01 = A[ta1 + 1]; a10 = A[ta2 + 1];
                b00 = B[tb1]; b01 = B[tb1 + 1]; b10 = B[tb1 + 2];

                
                c00 += a00 * b00; c01 += a00 * b01; c02 += a00 * b10;
                c10 += a01 * b00; c11 += a01 * b01; c12 += a01 * b10;
                c20 += a10 * b00; c21 += a10 * b01; c22 += a10 * b10;

                
                a00 = A[ta + 2]; a01 = A[ta1 + 2]; a10 = A[ta2 + 2];
                b00 = B[tb2]; b01 = B[tb2 + 1]; b10 = B[tb2 + 2];

                
                c00 += a00 * b00; c01 += a00 * b01; c02 += a00 * b10;
                c10 += a01 * b00; c11 += a01 * b01; c12 += a01 * b10;
                c20 += a10 * b00; c21 += a10 * b01; c22 += a10 * b10;
            }

            
            C[t] = c00; C[t + 1] = c01; C[t + 2] = c02;
            C[t1] = c10; C[t1 + 1] = c11; C[t1 + 2] = c12;
            C[t2] = c20; C[t2 + 1] = c21; C[t2 + 2] = c22;
        }
    }
}
