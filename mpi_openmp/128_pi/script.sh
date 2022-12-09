module load SpectrumMPI

mpicc m_128_pi.c -fopenmp -o  m_128_pi -O3 -lm

bsub < start_128_1_pi.lsf
bsub < start_128_2_pi.lsf
bsub < start_128_4_pi.lsf
bsub < start_128_8_pi.lsf