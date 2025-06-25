#include <mpi.h>
#include <iostream>
#include <cmath>
#include <vector>
#include <cstdlib>  // for atoi

using namespace std;

const int MAX_ITER = 1000;
const float TOL = 1e-5;

//applying boundary conditions
void apply_boundary_conditions(vector<vector<float> >& grid, int rank, int size, int local_rows, int N) {
    if (rank == 0) {
        for (int j = 0; j < N; ++j)
            grid[1][j] = 5.0f;  // top row
    }
    if (rank == size - 1) {
        for (int j = 0; j < N; ++j)
            grid[local_rows][j] = 5.0f;  // bottom row
    }
    for (int i = 1; i <= local_rows; ++i) {
        grid[i][0] = 0.0f;
        grid[i][N - 1] = 0.0f;
    }
}

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (argc != 2) {
        if (rank == 0) cerr << "Usage: mpirun -np <p> ./laplace_mpi <grid_size>\n";
        MPI_Finalize();
        return 1;
    }

    int N = atoi(argv[1]);
    if (N % size != 0) {
        if (rank == 0) cerr << "Grid size N must be divisible by number of processes.\n";
        MPI_Finalize();
        return 1;
    }

    int local_rows = N / size;
    int padded_rows = local_rows + 2; // includes ghost rows

    vector<vector<float> > V(padded_rows, vector<float>(N, 0.0f));
    vector<vector<float> > V_new = V;

    apply_boundary_conditions(V, rank, size, local_rows, N);
    apply_boundary_conditions(V_new, rank, size, local_rows, N);

    float local_error = 1.0f, global_error = 1.0f;
    int iter = 0;
    double start = MPI_Wtime();

    while (global_error > TOL && iter < MAX_ITER) {
        // Halo exchange
        if (rank != 0)
            MPI_Sendrecv(V[1].data(), N, MPI_FLOAT, rank - 1, 0,
                V[0].data(), N, MPI_FLOAT, rank - 1, 0,
                MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        if (rank != size - 1)
            MPI_Sendrecv(V[local_rows].data(), N, MPI_FLOAT, rank + 1, 0,
                V[local_rows + 1].data(), N, MPI_FLOAT, rank + 1, 0,
                MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        local_error = 0.0f;
        for (int i = 1; i <= local_rows; ++i) {
            for (int j = 1; j < N - 1; ++j) {
                V_new[i][j] = 0.25f * (
                    V[i + 1][j] + V[i - 1][j] + V[i][j + 1] + V[i][j - 1]);
                float diff = V_new[i][j] - V[i][j];
                local_error += diff * diff;
            }
        }

        // Swap
        swap(V, V_new);

        // Reduce error
        MPI_Allreduce(&local_error, &global_error, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
        global_error = sqrt(global_error);

        iter++;
    }

    double end = MPI_Wtime();

    if (rank == 0) {
        cout << "Grid size: " << N << " x " << N << endl;
        cout << "MPI Laplace Solver completed in " << iter << " iterations\n";
        cout << "Final error: " << global_error << "\n";
        cout << "Execution time: " << (end - start) * 1000 << " ms\n";
    }

    MPI_Finalize();
    return 0;
}
