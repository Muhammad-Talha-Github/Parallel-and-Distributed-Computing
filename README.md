# 🔬 Laplace Solver Using OpenMP and MPI

This project implements a parallel **Laplace Equation Solver** using two different approaches:

* **OpenMP** for shared-memory systems
* **MPI (Message Passing Interface)** for distributed-memory systems

The solution is benchmarked and tested on a High-Performance Computing (HPC) environment provided by **ScREC (Supercomputing Research and Education Center), NUST**, using the supercomputer **afrit.rcms.nust.edu.pk**.

---

## 📌 Problem Statement

> **Implement a Laplace Solver on HPC using MPI.**
> Test for correctness using smaller 2D matrices, and compare results with a single-threaded implementation.
> Gradually scale the problem and evaluate the performance of the MPI implementation against the OpenMP version.
> For large matrix sizes and iterations, implement proper inter-process communication (boundary exchange) using MPI.
> Highlight any performance bottlenecks encountered.

---

## 📁 Project Structure

```
Laplace-Solver/
│
├── laplace_openmp.cpp   # Laplace solver using OpenMP
├── laplace_mpi.cpp      # Laplace solver using MPI
├── Project_PDP          # A word file containing results
└── README.md            # Project documentation
```

---

## 🛠️ How to Compile & Run

### 🔷 OpenMP Version

```bash
g++ -fopenmp -O3 -o laplace_openmp laplace_openmp.cpp
./laplace_openmp
```

### 🔶 MPI Version

```bash
mpic++ -O3 -o laplace_mpi laplace_mpi.cpp
mpirun -np <num_processes> ./laplace_mpi
```

Replace `<num_processes>` with the number of MPI processes you want to use.

---

## 📈 Performance Analysis

* **Matrix sizes tested**: 100×100 to 1000×1000
* **Epochs**: Scaled up with matrix size (100 to 10,000+)
* **Speedup measured** for both implementations by comparing against serial performance
* **Inter-process communication** in MPI ensured boundary data exchange per epoch

---

## ⚙️ HPC Execution Context

The MPI implementation was executed on the **ScREC Supercomputer** using:

* **MPI Implementation**: Open MPI 1.6.2 (default on ScREC)
* **Node Access**: Via SSH to compute nodes
* **Job Submission**: Through `qsub` and `mpirun` as per ScREC's job scheduling policy
* **Supercomputer Specs**:

  * 272 CPU cores
  * NVidia Tesla S1070 GPU (unused in this project)
  * 132 Teraflops peak performance
