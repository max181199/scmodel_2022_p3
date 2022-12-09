module load SpectrumMPI

mpicc m_512_pi.c -fopenmp -o m_512_pi -O3 -lm

bsub < start_512_1_pi.lsf
bsub < start_512_2_pi.lsf
bsub < start_512_4_pi.lsf
bsub < start_512_8_pi.lsf