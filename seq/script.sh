module load SpectrumMPI

mpicc s_128_1.c -o s_128_1 -O3 -lm
mpicc s_128_pi.c -o s_128_pi -O3 -lm
mpicc s_256_1.c -o s_256_1 -O3 -lm
mpicc s_256_pi.c -o s_256_pi -O3 -lm
mpicc s_512_1.c -o s_512_1 -O3 -lm
mpicc s_512_pi.c -o s_512_pi -O3 -lm

mpisubmit.pl -p 1 ./s_128_1
mpisubmit.pl -p 1 ./s_128_pi 
mpisubmit.pl -p 1 ./s_256_1 
mpisubmit.pl -p 1 ./s_256_pi 
mpisubmit.pl -p 1 ./s_512_1 
mpisubmit.pl -p 1 ./s_512_pi  