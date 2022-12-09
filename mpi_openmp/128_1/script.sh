module load SpectrumMPI

mpicc m_128_1.c -fopenmp -o m_128_1 -O3 -lm

bsub < start_128_1_1.lsf
bsub < start_128_2_1.lsf
bsub < start_128_4_1.lsf
bsub < start_128_8_1.lsf