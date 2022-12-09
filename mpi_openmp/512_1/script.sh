module load SpectrumMPI

mpicc m_512_1.c -fopenmp -o m_512_1 -O3 -lm

bsub < start_512_1_1.lsf
bsub < start_512_2_1.lsf
bsub < start_512_4_1.lsf
bsub < start_512_8_1.lsf