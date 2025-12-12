#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// Declaración de función LAPACK para calcular eigenvalores seleccionados
extern void dsyevr_(char* jobz, char* range, char* uplo, int* n, 
                    double* a, int* lda, double* vl, double* vu,
                    int* il, int* iu, double* abstol, int* m,
                    double* w, double* z, int* ldz, int* isuppz,
                    double* work, int* lwork, int* iwork, int* liwork, int* info);

extern double dlamch_(char* cmach);

/**
 * Read covariance matrix from binary file
 * Format: [n (int32), m (int32), data (float32 array)]
 */
void read_covariance_binary(const char* filename, double** matrix, int* n) {
    FILE* fp = fopen(filename, "rb");
    if (!fp) {
        fprintf(stderr, "Error opening file %s\n", filename);
        exit(1);
    }
    
    // Read dimensions
    int dims[2];
    if (fread(dims, sizeof(int), 2, fp) != 2) {
        fprintf(stderr, "Error reading dimensions\n");
        exit(1);
    }
    
    *n = dims[0];
    int m = dims[1];
    
    if (*n != m) {
        fprintf(stderr, "Error: matrix is not square (%d x %d)\n", *n, m);
        exit(1);
    }
    
    printf("Matrix size detected: %d x %d\n", *n, *n);
    
    // Allocate memory
    size_t matrix_size = (size_t)(*n) * (*n);
    *matrix = (double*)malloc(matrix_size * sizeof(double));
    if (!*matrix) {
        fprintf(stderr, "Error allocating memory for matrix\n");
        exit(1);
    }
    
    // Read as float32, convert to double
    float* temp = (float*)malloc(matrix_size * sizeof(float));
    if (!temp) {
        fprintf(stderr, "Error allocating temporary buffer\n");
        exit(1);
    }
    
    size_t read_count = fread(temp, sizeof(float), matrix_size, fp);
    if (read_count != matrix_size) {
        fprintf(stderr, "Error reading matrix data (got %zu, expected %zu)\n", 
                read_count, matrix_size);
        exit(1);
    }
    
    // Convert to double
    for (size_t i = 0; i < matrix_size; i++) {
        (*matrix)[i] = (double)temp[i];
    }
    
    free(temp);
    fclose(fp);
    printf("Binary covariance matrix loaded successfully\n");
}

/**
 * Read covariance matrix from CSV file
 * 
 * Format: First line is header, then n x n matrix
 */
void read_covariance_matrix(const char* filename, double** matrix, int* n) {
    FILE* fp = fopen(filename, "r");
    if (!fp) {
        fprintf(stderr, "Error opening file %s\n", filename);
        exit(1);
    }
    
    // Buffer for reading lines
    size_t line_size = 50000000;  // 50MB buffer
    char* line = (char*)malloc(line_size);
    if (!line) {
        fprintf(stderr, "Error allocating memory for line buffer\n");
        exit(1);
    }
    
    // Read header and count variables
    if (!fgets(line, line_size, fp)) {
        fprintf(stderr, "Error reading header\n");
        exit(1);
    }
    
    *n = 1;
    for (size_t i = 0; line[i]; i++) {
        if (line[i] == ',') (*n)++;
    }
    
    printf("Matrix size detected: %d x %d\n", *n, *n);
    
    // Allocate memory for matrix
    size_t matrix_size = (size_t)(*n) * (*n);
    *matrix = (double*)malloc(matrix_size * sizeof(double));
    if (!*matrix) {
        fprintf(stderr, "Error allocating memory for matrix\n");
        exit(1);
    }
    
    // Read matrix data
    for (int i = 0; i < *n; i++) {
        if (!fgets(line, line_size, fp)) {
            fprintf(stderr, "Error reading matrix row %d\n", i);
            exit(1);
        }
        
        char* token = strtok(line, ",");
        for (int j = 0; j < *n && token != NULL; j++) {
            (*matrix)[i * (*n) + j] = atof(token);
            token = strtok(NULL, ",");
        }
        
        if ((i + 1) % 500 == 0) {
            printf("Reading matrix... %d/%d rows\n", i + 1, *n);
        }
    }
    
    fclose(fp);
    free(line);
    printf("Covariance matrix loaded successfully\n");
}

int main(int argc, char* argv[]) {
    if (argc < 5) {
        printf("Usage: %s covariance.(csv|bin) K eigenvalues.csv eigenvectors.csv\n", argv[0]);
        printf("\nArguments:\n");
        printf("  covariance.(csv|bin) : Input covariance matrix (CSV or binary format)\n");
        printf("  K                    : Number of eigenvalues to compute (the K largest)\n");
        printf("  eigenvalues.csv      : Output file for eigenvalues\n");
        printf("  eigenvectors.csv     : Output file for eigenvectors\n");
        printf("\nFormat detection:\n");
        printf("  - Files ending in .bin are treated as binary format\n");
        printf("  - All other files are treated as CSV format\n");
        return 1;
    }
    
    char* cov_file = argv[1];
    int K = atoi(argv[2]);
    char* eigenval_file = argv[3];
    char* eigenvec_file = argv[4];
    
    double* cov_matrix;
    int n;
    
    // Detect file format
    size_t len = strlen(cov_file);
    int is_binary = (len > 4 && strcmp(cov_file + len - 4, ".bin") == 0);
    
    printf("Reading covariance matrix...\n");
    if (is_binary) {
        printf("Detected binary format (.bin)\n");
        read_covariance_binary(cov_file, &cov_matrix, &n);
    } else {
        printf("Detected CSV format\n");
        read_covariance_matrix(cov_file, &cov_matrix, &n);
    }
    
    // Validate K
    if (K > n) {
        printf("Warning: K=%d is greater than n=%d. Using K=%d\n", K, n, n);
        K = n;
    }
    if (K <= 0) {
        fprintf(stderr, "Error: K must be greater than 0\n");
        return 1;
    }
    
    printf("Computing the %d largest eigenvalues...\n", K);
    
    // Prepare for LAPACK dsyevr
    char jobz = 'V';      // Compute eigenvalues and eigenvectors
    char range = 'I';     // Select by index
    char uplo = 'U';      // Upper triangle
    int lda = n;
    
    double vl = 0.0, vu = 0.0;  // Not used with range='I'
    int il = n - K + 1;          // Initial index (the K largest)
    int iu = n;                  // Final index
    
    // Tolerance for convergence
    char cmach = 'S';
    double abstol = 2.0 * dlamch_(&cmach);
    
    int m;  // Number of eigenvalues found
    double* w = (double*)malloc(n * sizeof(double));
    double* z = (double*)malloc((size_t)n * K * sizeof(double));
    int ldz = n;
    int* isuppz = (int*)malloc(2 * K * sizeof(int));
    
    if (!w || !z || !isuppz) {
        fprintf(stderr, "Error allocating memory for eigenvalues/eigenvectors\n");
        return 1;
    }
    
    // Query optimal workspace size
    double work_query;
    int lwork = -1;
    int iwork_query;
    int liwork = -1;
    int info;
    
    printf("Querying optimal workspace size...\n");
    dsyevr_(&jobz, &range, &uplo, &n, cov_matrix, &lda, &vl, &vu, &il, &iu,
            &abstol, &m, w, z, &ldz, isuppz, &work_query, &lwork, 
            &iwork_query, &liwork, &info);
    
    lwork = (int)work_query;
    liwork = iwork_query;
    printf("Workspace required: %d doubles (%.2f MB), %d ints (%.2f MB)\n", 
           lwork, lwork * sizeof(double) / 1024.0 / 1024.0,
           liwork, liwork * sizeof(int) / 1024.0 / 1024.0);
    
    double* work = (double*)malloc(lwork * sizeof(double));
    int* iwork = (int*)malloc(liwork * sizeof(int));
    
    if (!work || !iwork) {
        fprintf(stderr, "Error allocating memory for workspace\n");
        return 1;
    }
    
    // Compute eigenvalues and eigenvectors
    printf("Computing eigenvalues and eigenvectors (this may take a while)...\n");
    dsyevr_(&jobz, &range, &uplo, &n, cov_matrix, &lda, &vl, &vu, &il, &iu,
            &abstol, &m, w, z, &ldz, isuppz, work, &lwork, iwork, &liwork, &info);
    
    if (info != 0) {
        fprintf(stderr, "Error in dsyevr: info = %d\n", info);
        return 1;
    }
    
    printf("Eigenvalues found: %d\n", m);
    printf("Saving results...\n");
    
    // Save eigenvalues (in descending order)
    FILE* fp_val = fopen(eigenval_file, "w");
    if (!fp_val) {
        fprintf(stderr, "Error opening %s\n", eigenval_file);
        return 1;
    }
    
    fprintf(fp_val, "index,eigenvalue\n");
    for (int i = m - 1; i >= 0; i--) {
        fprintf(fp_val, "%d,%.6e\n", m - i, w[i]);
    }
    fclose(fp_val);
    
    // Save eigenvectors in ascending order (ev_1 = largest eigenvalue)
    FILE* fp_vec = fopen(eigenvec_file, "w");
    if (!fp_vec) {
        fprintf(stderr, "Error opening %s\n", eigenvec_file);
        return 1;
    }
    
    // Header: ev_1, ev_2, ..., ev_m
    for (int i = 1; i <= m; i++) {
        fprintf(fp_vec, "ev_%d", i);
        if (i < m) fprintf(fp_vec, ",");
    }
    fprintf(fp_vec, "\n");
    
    // Data: each row is a component, each column an eigenvector
    for (int i = 0; i < n; i++) {
        for (int j = m - 1; j >= 0; j--) {
            fprintf(fp_vec, "%.6e", z[j * n + i]);
            if (j > 0) fprintf(fp_vec, ",");
        }
        fprintf(fp_vec, "\n");
        
        if ((i + 1) % 500 == 0) {
            printf("Saving eigenvectors... %d/%d\n", i + 1, n);
        }
    }
    fclose(fp_vec);
    
    printf("\nResults saved:\n");
    printf("  - Eigenvalues: %s (%d values)\n", eigenval_file, m);
    printf("  - Eigenvectors: %s (%d vectors of dimension %d)\n", 
           eigenvec_file, m, n);
    printf("  - ev_1 corresponds to the largest eigenvalue\n");
    
    free(cov_matrix);
    free(w);
    free(z);
    free(isuppz);
    free(work);
    free(iwork);
    
    return 0;
}