# Bioinformatics Genome Similarity Tool

This project implements a bioinformatics application to compute the genetic similarity between bacterial genomes. The program takes genome sequences of multiple bacteria as input, analyzes 6-mers (substrings of length 6), and determines a correlation between them. It utilizes parallel computing to improve efficiency and scalability for large datasets.

---

## Features
- **Input Handling**: Reads bacterial genome sequences from files.
- **6-mer Analysis**: Computes frequency vectors for 6-mers.
- **Stochastic Model**: Estimates expected occurrences of k-mers.
- **Parallel Processing**: Uses OpenMP to speed up computations by distributing tasks across multiple threads.

---

## How It Works
### Key Components
1. **Initialization**: Sets parameters for k-mer sizes.
2. **File Reading**: Reads bacterial genome sequences from input files.
3. **Comparison**: Computes similarity scores between bacterial genomes using a stochastic model.

### Parallelization
The application leverages parallelization using OpenMP:
- **Nested Loops**: The primary computational workload lies in comparing bacterial pairs. These comparisons are distributed across threads.
- **Dynamic Scheduling**: Balances workload efficiently across threads.
- **Critical Sections**: Ensures thread-safe operations for output.

---

## Performance
### Results
Performance improvements with OpenMP are shown below:
| Threads | Execution Time (seconds) | Speedup |
|---------|--------------------------|---------|
| 1       | 944.67                   | 1.46    |
| 2       | 484                      | 2.85    |
| 4       | 247.67                   | 5.57    |
| 6       | 190                      | 7.26    |
| 8       | 173.33                   | 7.96    |
| 10      | 152.67                   | 9.03    |
| 12      | 145.33                   | 9.49    |
| 14      | 142.67                   | 9.67    |
| 16      | 142.33                   | 9.69    |

The results show a near-linear speedup up to 10 threads, with diminishing returns beyond that point.

### Profiling Tools
- **Visual Studio**: Used for enabling OpenMP and profiling performance.
- **Dynamic Scheduling**: Optimizes load balancing across threads.

---

## Code Modifications
1. **Memory Management**: Transitioned to `std::vector` for efficient and scalable memory handling.
2. **Optimized Functions**:
   - `CompareAllBacteria`: Reduced dynamic memory operations.
   - `CompareBacteria`: Replaced modular calls with inline calculations for reduced overhead.
3. **File Reading**: Buffered reads for improved efficiency.
4. **Parallel Directives**:
   - `#pragma omp parallel for`: Parallelized bacterial pair comparisons.
   - `#pragma omp critical`: Ensured thread-safe outputs.

---

## Requirements
- **Compiler**: C++ compiler with OpenMP support (e.g., Visual Studio, GCC).
- **Operating System**: 64-bit environment.
- **Dependencies**: OpenMP library.

---

## Installation
1. Clone the repository:
   ```bash
   git clone https://github.com/username/bioinformatics-genome-similarity.git
   ```
2. Build the project with OpenMP enabled:
   ```bash
   g++ -fopenmp -o genome_similarity main.cpp
   ```
3. Run the program:
   ```bash
   ./genome_similarity <input_file>
   ```

---

## Future Enhancements
- **Task-level Parallelism**: Explore more granular parallelization.
- **Hybrid Parallelization**: Combine OpenMP with MPI for distributed processing.

---

## Author
**Gustavo Amorim De Almeida**

For questions or suggestions, feel free to open an issue or contact the author directly.

