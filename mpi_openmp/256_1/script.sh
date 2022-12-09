module load SpectrumMPI

mpicc m_256_1.c -fopenmp -o m_256_1 -O3 -lm

bsub < start_256_1_1.lsf
bsub < start_256_2_1.lsf
bsub < start_256_4_1.lsf
bsub < start_256_8_1.lsf