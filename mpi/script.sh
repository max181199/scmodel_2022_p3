module load SpectrumMPI

mpicc m_128_1.c -o m_128_1 -O3 -lm
mpicc m_128_pi.c -o m_128_pi -O3 -lm
mpicc m_256_1.c -o m_256_1 -O3 -lm
mpicc m_256_pi.c -o m_256_pi -O3 -lm
mpicc m_512_1.c -o m_512_1 -O3 -lm
mpicc m_512_pi.c -o m_512_pi -O3 -lm

mpisubmit.pl -p 1 ./m_128_1
mpisubmit.pl -p 1 ./m_128_pi 
mpisubmit.pl -p 1 ./m_256_1 
mpisubmit.pl -p 1 ./m_256_pi 
mpisubmit.pl -p 1 ./m_512_1 
mpisubmit.pl -p 1 ./m_512_pi  

mpisubmit.pl -p 4 ./m_128_1
mpisubmit.pl -p 4 ./m_128_pi 
mpisubmit.pl -p 4 ./m_256_1 
mpisubmit.pl -p 4 ./m_256_pi 
mpisubmit.pl -p 4 ./m_512_1 
mpisubmit.pl -p 4 ./m_512_pi  

mpisubmit.pl -p 8 ./m_128_1
mpisubmit.pl -p 8 ./m_128_pi 
mpisubmit.pl -p 8 ./m_256_1 
mpisubmit.pl -p 8 ./m_256_pi 
mpisubmit.pl -p 8 ./m_512_1 
mpisubmit.pl -p 8 ./m_512_pi  

mpisubmit.pl -p 16 ./m_128_1
mpisubmit.pl -p 16 ./m_128_pi 
mpisubmit.pl -p 16 ./m_256_1 
mpisubmit.pl -p 16 ./m_256_pi 
mpisubmit.pl -p 16 ./m_512_1 
mpisubmit.pl -p 16 ./m_512_pi  

mpisubmit.pl -p 32 ./m_128_1
mpisubmit.pl -p 32 ./m_128_pi 
mpisubmit.pl -p 32 ./m_256_1 
mpisubmit.pl -p 32 ./m_256_pi 
mpisubmit.pl -p 32 ./m_512_1 
mpisubmit.pl -p 32 ./m_512_pi  