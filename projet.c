git remote add origin https : // github.com/shirelepfl/CMT-Project.git
                              git branch -
                              M main
                                  git push -
                              u origin main

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define IDX(i, j, k) ((k) * Ny * Nx + (j) * Nx + (i))

                              static void thomas_solve(int N, const double *a, const double *b, const double *c, double *d);
static void build_tridiag_neumann_1D(int N, double alpha, double *a, double *b, double *c);
void step_IMEX_ADI3D(double *u, double *tmp1, double *tmp2,
                     int Nx, int Ny, int Nz,
                     double dx, double dy, double dz,
                     double dt, double D, double r, double K, double lambda);
double compute_front_radius(const double *u, int Nx, int Ny, int Nz,
                            double dx, double dy, double dz, double K);

typedef struct
{
    char name[32];
    double r, D, K, lambda;
} Case;

long count_affected(const double *u, int Nx, int Ny, int Nz, double threshold)
{
    long count = 0;
    int Ntot = Nx * Ny * Nz;

    for (int p = 0; p < Ntot; ++p)
    {
        if (u[p] > threshold)
        {
            count++;
        }
    }
    return count;
}

int main(void)
{
    // Baseline parameters (same for both profiles)
    const double r_0 = 0.03;      // days^-1
    const double D_0 = 0.5;       // mm^2/days
    const double lambda_0 = 0.01; // days^-1
    const double K_0 = 1e6;       // cells/mm^3 // vrai valeur K_0 = const double K_0 = 1e6;

    const double T = 122.0; // days
    int Nx = 60, Ny = 60, Nz = 60;
    const double L = 120.0; // mm

    const double dx = L / (Nx - 1);
    const double dy = L / (Ny - 1);
    const double dz = L / (Nz - 1);
    double dt = 5;
    int Nt = (int)round(T / dt);

    // Create the two cases
    Case cases[2];

    // Healthy
    strcpy(cases[0].name, "healthy");
    cases[0].r = 0.8 * r_0;
    cases[0].D = 0.9 * D_0;
    cases[0].K = 0.85 * K_0;
    cases[0].lambda = 1.4 * lambda_0;

    // Smoker
    strcpy(cases[1].name, "smoker");
    cases[1].r = 1.5 * r_0;
    cases[1].D = 1.4 * D_0;
    cases[1].K = 1.3 * K_0;
    cases[1].lambda = 1.7 * lambda_0;

    const size_t Nxyz = (size_t)Nx * Ny * Nz;

    // Allocate 3D arrays in flattened 1D form
    double *u = malloc(Nxyz * sizeof(double));
    double *tmp1 = malloc(Nxyz * sizeof(double));
    double *tmp2 = malloc(Nxyz * sizeof(double));

    // Pré-calcul des coordonnées (évite i*dx etc. dans la triple boucle)
    double *x = malloc(Nx * sizeof(double));
    double *y = malloc(Ny * sizeof(double));
    double *z = malloc(Nz * sizeof(double));

    for (int i = 0; i < Nx; ++i)
        x[i] = i * dx;
    for (int j = 0; j < Ny; ++j)
        y[j] = j * dy;
    for (int k = 0; k < Nz; ++k)
        z[k] = k * dz;

    const double xc = L / 2.0, yc = L / 2.0, zc = L / 2.0;
    const double A0 = 1e5;
    const double sigma = 2.0;
    const double inv_2sigma2 = 1.0 / (2.0 * sigma * sigma);

    for (int ic = 0; ic < 2; ++ic)
    {
        double r = cases[ic].r;
        double D = cases[ic].D;
        double K = cases[ic].K;
        double lambda = cases[ic].lambda;

        // Initial Gaussian tumor (parallélisée)
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

        printf("=== Cas %s ===\n", cases[ic].name);

        double *t_vals = malloc((Nt + 1) * sizeof(double));
        long *count_vals = malloc((Nt + 1) * sizeof(long));

        // Time loop (le plus coûteux → à optimiser dans step_IMEX_ADI3D aussi)
        for (int n = 0; n <= Nt; ++n)
        {
            step_IMEX_ADI3D(u, tmp1, tmp2,
                            Nx, Ny, Nz,
                            dx, dy, dz,
                            dt, D, r, K, lambda);
            double t = n * dt;
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

            // store results for this time step in csv
            t_vals[n] = t;
            count_vals[n] = count;
        }

        char filename[64];
        snprintf(filename, sizeof(filename), "%s_results.csv", cases[ic].name);

        FILE *f = fopen(filename, "w");

        // write header
        fprintf(f, "t,count\n");

        // write values
        for (int n = 0; n <= Nt; ++n)
        {
            fprintf(f, "%f,%ld\n", t_vals[n], count_vals[n]);
        }

        fclose(f);

        free(t_vals);
        free(count_vals);
    }
    free(u);
    free(tmp1);
    free(tmp2);
    free(x);
    free(y);
    free(z);

    return 0;
}

// on résout Thomas en 3D, c'est une matrice tridiagonale A donc le triangle supérieur est a, puis la diagonale est b et puis dessous c'est c
// et d est la solution, forme est Ax = d

static void thomas_solve(int N, const double *a, const double *b, const double *c, double *d)
{
    double cp[N - 1];
    double dp[N]; // memory allocation function to save our tables

    cp[0] = c[0] / b[0];
    dp[0] = d[0] / b[0];
    for (int i = 1; i < N - 1; i++)
    {
        double denom = b[i] - a[i] * cp[i - 1];
        cp[i] = c[i] / denom;
        dp[i] = (d[i] - a[i] * dp[i - 1]) / denom;
    } // cp[i] and dp[i] are the modified coefficients of the superior triangular system
    dp[N - 1] = (d[N - 1] - a[N - 1] * dp[N - 2]) / (b[N - 1] - a[N - 1] * cp[N - 2]);
    d[N - 1] = dp[N - 1]; // solving the superior triangular system
    for (int i = N - 2; i >= 0; i--)
        d[i] = dp[i] - cp[i] * d[i + 1]; // results are kept in d which is the solution matrix
}

static void build_tridiag_neumann_1D(int N, double alpha, double *a, double *b, double *c)
{
    for (int i = 0; i < N; i++)
    {
        a[i] = 0;
        b[i] = 1;
        c[i] = 0;
    }
    if (N == 1)
    {
        b[0] = 1;
        return;
    }
    b[0] = 1 + 2 * alpha;
    c[0] = -2 * alpha;
    for (int i = 1; i < N - 1; i++)
    {
        a[i] = -alpha;
        b[i] = 1 + 2 * alpha;
        c[i] = -alpha;
    }
    a[N - 1] = -2 * alpha;
    b[N - 1] = 1 + 2 * alpha;
}

void step_IMEX_ADI3D(double *restrict u, double *restrict tmp1, double *restrict tmp2,
                     int Nx, int Ny, int Nz,
                     double dx, double dy, double dz,
                     double dt, double D, double r, double K, double lambda)
{
    double ax = dt * D / (dx * dx);
    double ay = dt * D / (dy * dy);
    double az = dt * D / (dz * dz);

    // 1. réaction + YZ explicites
    for (int k = 0; k < Nz; k++)
    {
        for (int j = 0; j < Ny; j++)
        {
            for (int i = 0; i < Nx; i++)
            {
                int p = IDX(i, j, k);
                double un = u[p];

                // Terme de réaction (logistique - lambda*u)
                double f = r * un * (1.0 - un / K) - lambda * un;

                // Indices avec conditions de Neumann en y et z
                int jm = (j == 0) ? 1 : j - 1;
                int jp = (j == Ny - 1) ? Ny - 2 : j + 1;
                int km = (k == 0) ? 1 : k - 1;
                int kp = (k == Nz - 1) ? Nz - 2 : k + 1;

                double uy = u[IDX(i, jp, k)] - 2.0 * un + u[IDX(i, jm, k)];
                double uz = u[IDX(i, j, kp)] - 2.0 * un + u[IDX(i, j, km)];

                tmp1[p] = un + dt * f + ay * uy + az * uz;
            }
        }
    }

    // 2. X-implicite : on construit la tridiagonale sur la pile
    {
        double ax_a[Nx], ax_b[Nx], ax_c[Nx];
        build_tridiag_neumann_1D(Nx, ax, ax_a, ax_b, ax_c);

        for (int k = 0; k < Nz; k++)
        {
            for (int j = 0; j < Ny; j++)
            {
                // extraire la ligne en x
                for (int i = 0; i < Nx; i++)
                    tmp2[i] = tmp1[IDX(i, j, k)];

                // résoudre
                thomas_solve(Nx, ax_a, ax_b, ax_c, tmp2);

                // réinjecter
                for (int i = 0; i < Nx; i++)
                    tmp1[IDX(i, j, k)] = tmp2[i];
            }
        }
    }

    // 3. Y-implicite
    {
        double ay_a[Ny], ay_b[Ny], ay_c[Ny];
        build_tridiag_neumann_1D(Ny, ay, ay_a, ay_b, ay_c);

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

    // 4. Z-implicite
    {
        double az_a[Nz], az_b[Nz], az_c[Nz];
        build_tridiag_neumann_1D(Nz, az, az_a, az_b, az_c);

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

// volume où u>=K/2  (utile pour front moyen)
// calcule un rayon moyen du front tumoral dans le domaine 3D
double compute_front_radius(const double *u, int Nx, int Ny, int Nz,
                            double dx, double dy, double dz, double K)
{
    double seuil = 0.75 * K; // seuil de densité -> moitié de la capacité maximale
    double V = 0;            // initialisation du volume
    for (int i = 0; i < Nx * Ny * Nz; i++)
        if (u[i] >= seuil)
            V += 1.0;
    V *= dx * dy * dz;
    double R = pow(3.0 * V / (4.0 * M_PI), 1.0 / 3.0); // Formule inversée de V = (4/3) π R³ pour obtenir R,rayon sphérique correspondant à ce volume
    return R;
}
