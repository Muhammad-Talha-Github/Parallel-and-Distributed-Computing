#include <iostream>
#include <cmath>
#include <vector>
#include <cstdlib>
#include <omp.h>

using namespace std;

const int MAX_ITER = 1000;
const float TOL = 1e-5;

//applying boundary conditions
void apply_boundary_conditions(vector<vector<float> >& grid, int N) {
    for (int j = 0; j < N; ++j) {
        grid[0][j] = 5.0f;          // top row
        grid[N - 1][j] = 5.0f;      // bottom row
    }
    for (int i = 0; i < N; ++i) {
        grid[i][0] = 0.0f;          // left
        grid[i][N - 1] = 0.0f;      // right
    }
}

int main(int argc, char** argv) {
    if (argc != 2) {
        cerr << "Usage: ./laplace_openmp <grid_size>\n";
        return 1;
    }

    int N = atoi(argv[1]);
    vector<vector<float> > V(N, vector<float>(N, 0.0f));
    vector<vector<float> > V_new = V;

    apply_boundary_conditions(V, N);
    apply_boundary_conditions(V_new, N);

    float error = 1.0f;
    int iter = 0;
    double start = omp_get_wtime();

    while (error > TOL && iter < MAX_ITER) {
        error = 0.0f;

#pragma omp parallel for reduction(+:error) collapse(2)
        for (int i = 1; i < N - 1; ++i) {
            for (int j = 1; j < N - 1; ++j) {
                V_new[i][j] = 0.25f * (
                    V[i + 1][j] + V[i - 1][j] + V[i][j + 1] + V[i][j - 1]);
                float diff = V_new[i][j] - V[i][j];
                error += diff * diff;
            }
        }

        swap(V, V_new);
        error = sqrt(error);
        iter++;
    }

    double end = omp_get_wtime();

    cout << "Grid size: " << N << " x " << N << endl;
    cout << "OpenMP Laplace Solver completed in " << iter << " iterations\n";
    cout << "Final error: " << error << "\n";
    cout << "Execution time: " << (end - start) * 1000 << " ms\n";

    return 0;
}
