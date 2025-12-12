git remote add origin https://github.com/shirelepfl/CMT-Project.git
git branch -M main
git push -u origin main

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

// This is because we store the 3D field u in a 1D array, each point in our domain will be stored at a certain index in the 1d array; 
// the index of each number is the IDX, with k being slices, y being rows and x being columns
#define IDX(i, j, k) ((k) * Ny * Nx + (j) * Nx + (i))

// Forward declarations of functions that we are gonna use in main()
static void thomas_solve(int N, const double *a, const double *b, const double *c, double *d);
static void build_tridiag_matrix_1D(int N, double alpha, double *a, double *b, double *c);
void IMEX_ADI_3D(double *u, double *tmp1, double *tmp2,
                     int Nx, int Ny, int Nz,
                     double dx, double dy, double dz,
                     double dt, double D, double r, double K, double lambda);

// Structure used to store one set of biological parameters
typedef struct
{
    char name[32];
    double r, D, K, lambda;
} Case;

int main(void)
{
    // Baseline parameters (same for both profiles)
    const double r_0 = 0.03;      // average net proliferation rate days^-1
    const double D_0 = 0.5;       // average effective diffusion mm^2/days
    const double lambda_0 = 0.01; // average net immune clearance term days^-1
    const double K_0 = 1e6;       // average carrying capacity cells/mm^3

    // Total simulation time (days)
    const double T = 120.0; 

    // Number of grid points in each dimension
    int Nx = 60, Ny = 60, Nz = 60;

    // Physical domain size (mm)
    const double L = 120.0;

    // Spatial grid spacing
    const double dx = L / (Nx - 1);
    const double dy = L / (Ny - 1);
    const double dz = L / (Nz - 1);

    // Time step and number of iterations
    double dt = 5;  // every 5 days
    int Nt = (int)round(T / dt);

    // Create the two cases: healthy and smoker
    Case cases[2];

    // Healthy lung parameters → slower tumor growth, stronger immune clearance
    strcpy(cases[0].name, "healthy");
    cases[0].r = 0.8 * r_0;
    cases[0].D = 0.9 * D_0;
    cases[0].K = 0.85 * K_0;
    cases[0].lambda = 1.4 * lambda_0;

    // Smoker lung parameters → faster growth, weaker immune clearance
    strcpy(cases[1].name, "smoker");
    cases[1].r = 1.5 * r_0;
    cases[1].D = 1.4 * D_0;
    cases[1].K = 1.3 * K_0;
    cases[1].lambda = 1.7 * lambda_0;

    const size_t Nxyz = (size_t)Nx * Ny * Nz;

    // Allocate the 3D solution arrays (stored as 1D)
    // u = current tumor cell density
    // tmp1, tmp2 = work arrays used in ADI steps
    double *u = malloc(Nxyz * sizeof(double));
    double *tmp1 = malloc(Nxyz * sizeof(double));
    double *tmp2 = malloc(Nxyz * sizeof(double));

    // Precompute coordinates to avoid recomputing i*dx many times
    double *x = malloc(Nx * sizeof(double));
    double *y = malloc(Ny * sizeof(double));
    double *z = malloc(Nz * sizeof(double));

    for (int i = 0; i < Nx; ++i) x[i] = i * dx;
    for (int j = 0; j < Ny; ++j) y[j] = j * dy;
    for (int k = 0; k < Nz; ++k) z[k] = k * dz;
  
    // Center of the initial tumor
    const double xc = L / 2.0, yc = L / 2.0, zc = L / 2.0;

    // Parameters for the initial 3D Gaussian tumor
    const double A0 = 1e5;
    const double sigma = 2.0;
    const double inv_2sigma2 = 1.0 / (2.0 * sigma * sigma);

    // Loop over the two scenarios
    for (int ic = 0; ic < 2; ++ic)
    {
        double r = cases[ic].r;
        double D = cases[ic].D;
        double K = cases[ic].K;
        double lambda = cases[ic].lambda;

        // Build the initial Gaussian tumor profile
        for (int k = 0; k < Nz; ++k)
        {
            for (int j = 0; j < Ny; ++j)
            {
                for (int i = 0; i < Nx; ++i)
                {
                    double dx_ = x[i] - xc;
                    double dy_ = y[j] - yc;
                    double dz_ = z[k] - zc;
                    double d2 = dx_ * dx_ + dy_ * dy_ + dz_ * dz_;
                    u[IDX(i, j, k)] = A0 * exp(-d2 * inv_2sigma2);
                }
            }
        }

        printf("=== %s case ===\n", cases[ic].name);

        // Arrays to store the results for the CSV output
        double *t_vals = malloc((Nt + 1) * sizeof(double));
        long *count_vals = malloc((Nt + 1) * sizeof(long));

        // Time-stepping loop
        for (int n = 0; n <= Nt; ++n)
        {
            // Perform one IMEX-ADI time step
            IMEX_ADI_3D(u, tmp1, tmp2,
                            Nx, Ny, Nz,
                            dx, dy, dz,
                            dt, D, r, K, lambda);
            double t = n * dt;

            // Count number of voxels with high density (u > 1 here)
            long count = 0;
            for (int k = 0; k < Nz; ++k)
            {
                for (int j = 0; j < Ny; ++j)
                {
                    for (int i = 0; i < Nx; ++i)
                    {
                        double val = u[IDX(i, j, k)];
                        if (val > 1)
                        { // ou > 1.0, ou > K*0.5, etc.
                            count++;
                        }
                    }
                }
            }

            // Store the results for this time step in CSV
            t_vals[n] = t;
            count_vals[n] = count;
        }

        // Save results into a CSV file named "<case>_results.csv"
        char filename[64];
        snprintf(filename, sizeof(filename), "%s_results.csv", cases[ic].name);

        FILE *f = fopen(filename, "w");

        // Write header
        fprintf(f, "t,count\n");

        // Write values
        for (int n = 0; n <= Nt; ++n)
        {
            fprintf(f, "%f,%ld\n", t_vals[n], count_vals[n]);
        }

        fclose(f);

        free(t_vals);
        free(count_vals);
    }

    // Free all allocated memory
    free(u);
    free(tmp1);
    free(tmp2);
    free(x);
    free(y);
    free(z);

    return 0;
}

/* 
Thomas algorithm solver
This solves a tridiagonal linear system A x = d,
where A has three diagonals: a (lower), b (main), c (upper).
We overwrite d[] with the solution. 
*/
static void thomas_solve(int N, const double *a, const double *b, const double *c, double *d)
{    
    // Temporary arrays (forward-sweep coefficients)
    double cp[N - 1];
    double dp[N]; // memory allocation function to save our tables

    // First element
    cp[0] = c[0] / b[0];
    dp[0] = d[0] / b[0];

    // Forward sweep: eliminate the lower diagonal
    for (int i = 1; i < N - 1; i++)
    {
        double denom = b[i] - a[i] * cp[i - 1];
        cp[i] = c[i] / denom;
        dp[i] = (d[i] - a[i] * dp[i - 1]) / denom;
    } 
    // cp[i] and dp[i] are the modified coefficients of the superior triangular system

    // Last element of forward sweep
    dp[N - 1] = (d[N - 1] - a[N - 1] * dp[N - 2]) / (b[N - 1] - a[N - 1] * cp[N - 2]);

    // Backward substitution
    d[N - 1] = dp[N - 1]; // solving the superior triangular system
    for (int i = N - 2; i >= 0; i--)
        d[i] = dp[i] - cp[i] * d[i + 1];   // results are kept in d which is the solution matrix
}

/*
Build the coefficients of a 1D tridiagonal matrix with Neumann BCs.
The matrix corresponds to (I - alpha * Laplacian_1D) under Neumann BCs.
a = lower diagonal, b = main diagonal, c = upper diagonal
*/

static void build_tridiag_matrix_1D(int N, double alpha, double *a, double *b, double *c)
{    
    // Initialize everything to identity matrix
    for (int i = 0; i < N; i++){
        a[i] = 0;
        b[i] = 1;
        c[i] = 0;
    }

    // Special cases for Neumann boundaries
    if (N == 1){
        b[0] = 1;
        return;
    }

    // First row (left Neumann BC)
    b[0] = 1 + 2 * alpha;
    c[0] = -2 * alpha;

    // Interior rows
    for (int i = 1; i < N - 1; i++){
        a[i] = -alpha;
        b[i] = 1 + 2 * alpha;
        c[i] = -alpha;
    }

    // Last row (right Neumann BC)
    a[N - 1] = -2 * alpha;
    b[N - 1] = 1 + 2 * alpha;
}

/*
One IMEX-ADI time step in 3D.

IMEX = IMplicit-EXplicit: 
treat diffusion implicitly, treat reaction explicitly.

ADI = Alternating Direction Implicit: 
instead of solving a big 3D implicit system (costly), we split the implicit diffusion solve into 
three 1-D implicit solves: first X, then Y, then Z. 
Each of those is a tridiagonal system solvable efficiently with Thomas.

u    = input/output array (3D field)
tmp1 = first temporary array
tmp2 = second temporary array
*/
void IMEX_ADI_3D(double *restrict u, double *restrict tmp1, double *restrict tmp2,
                     int Nx, int Ny, int Nz,
                     double dx, double dy, double dz,
                     double dt, double D, double r, double K, double lambda)
{   
    // Diffusion coefficients in each direction
    double ax = dt * D / (dx * dx);
    double ay = dt * D / (dy * dy);
    double az = dt * D / (dz * dz);

    // 1. Explicit reaction + explicit Y and Z diffusion
    for (int k = 0; k < Nz; k++)
    {
        for (int j = 0; j < Ny; j++) // for loop on each axis going to the border of our domain
        {
            for (int i = 0; i < Nx; i++)
            {
                int p = IDX(i, j, k);  // flattened index
                double un = u[p];     // current value

                // Logistic reaction term
                double f = r * un * (1.0 - un / K) - lambda * un; // this is the term without the Lagrangian, it is the net local growth rate of tumor cells

                // Find neighbors in y and z with Neumann BC
                // boundaries for the Neumann condition over y(j) and z(k)
                int jm = (j == 0) ? 1 : j - 1;
                int jp = (j == Ny - 1) ? Ny - 2 : j + 1;
                int km = (k == 0) ? 1 : k - 1;
                int kp = (k == Nz - 1) ? Nz - 2 : k + 1;

                double uy = u[IDX(i, jp, k)] - 2.0 * un + u[IDX(i, jm, k)];
                double uz = u[IDX(i, j, kp)] - 2.0 * un + u[IDX(i, j, km)];
                // approximation to solve the first part of FK

                tmp1[p] = un + dt * f + ay * uy + az * uz; // the new density of the local tumor cell from the last one (un) because it is a suite
            }
        }
    }

  /* 
  2. Implicit solve in X direction
 since 3d solving is difficult, the IMEX ADI method does the thomas solve in each direction then recomputes 
 */
    {
        double ax_a[Nx], ax_b[Nx], ax_c[Nx]; // the three diagonals: a is low, b is main and c is upper 
        build_tridiag_matrix_1D(Nx, ax, ax_a, ax_b, ax_c);

        for (int k = 0; k < Nz; k++)
        {
            for (int j = 0; j < Ny; j++)
            {
                // Extracting the row along x
                for (int i = 0; i < Nx; i++) // solve i for every possible j and k 
                    tmp2[i] = tmp1[IDX(i, j, k)]; // take the current value in the 3d compressed array and store it in a 1d array

                // Solve the tridiagonal system
                thomas_solve(Nx, ax_a, ax_b, ax_c, tmp2);

                // Inject solution back
                for (int i = 0; i < Nx; i++)
                    tmp1[IDX(i, j, k)] = tmp2[i];
            }
        }
    }

    // 3. Implicit solve in Y direction
    {
        double ay_a[Ny], ay_b[Ny], ay_c[Ny];
        build_tridiag_matrix_1D(Ny, ay, ay_a, ay_b, ay_c);

        for (int k = 0; k < Nz; k++)
        {
            for (int i = 0; i < Nx; i++)
            {
                for (int j = 0; j < Ny; j++)
                    tmp2[j] = tmp1[IDX(i, j, k)];

                thomas_solve(Ny, ay_a, ay_b, ay_c, tmp2);

                for (int j = 0; j < Ny; j++)
                    tmp1[IDX(i, j, k)] = tmp2[j];
            }
        }
    }

    // 4. Implicit solve in Z direction
    {
        double az_a[Nz], az_b[Nz], az_c[Nz];
        build_tridiag_matrix_1D(Nz, az, az_a, az_b, az_c);

        for (int j = 0; j < Ny; j++)
        {
            for (int i = 0; i < Nx; i++)
            {
                for (int k = 0; k < Nz; k++)
                    tmp2[k] = tmp1[IDX(i, j, k)];

                thomas_solve(Nz, az_a, az_b, az_c, tmp2);

                for (int k = 0; k < Nz; k++)
                    u[IDX(i, j, k)] = tmp2[k];
            }
        }
    }
}








