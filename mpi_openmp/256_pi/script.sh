module load SpectrumMPI

mpicc m_256_pi.c -fopenmp -o m_256_pi -O3 -lm

bsub < start_256_1_pi.lsf
bsub < start_256_2_pi.lsf
bsub < start_256_4_pi.lsf
bsub < start_256_8_pi.lsf