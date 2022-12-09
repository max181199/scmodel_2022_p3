#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <mpi.h>

#define N 127
#define K 220
#define M_PI_L 3.14159265358

typedef struct Gin {
  double lx;
  double ly;
  double lz;
  double T;
  double hx;
  double hy;
  double hz;
  double tau;
} gin;


double fi (gin * g, int i, int j, int k){
  return sin((g->hx * (double)i) * M_PI_L / g->lx) * 
          sin((g->hy *  (double)j) * 2.0 * M_PI_L / g->ly) * 
            sin((g->hz *  (double)k) * 3.0 * M_PI_L / g->lz);
}

double u (gin * g, int i, int j, int k, int t){
  return fi(g, i, j ,k) * cos(
    M_PI_L *
    sqrt(
      1.0/(g->lx * g->lx) + 
      4.0/(g->ly * g->ly) + 
      9.0/(g->lz * g->lz)
    ) * (g->tau * (double)t)
  );
} 

void u_0_init (double (*a)[N + 1][N + 1], gin * g){
  for (int i = 0; i <= N; i++){
    for (int j = 0; j <= N; j++){
      for (int k = 0; k <= N; k++){
        a[i][j][k] = fi(g, i, j, k);
      }
    }
  }
}

void u_1_init (double (*a)[N + 1][N + 1], double (*pa)[N + 1][N + 1], gin * g){

  for (int i = 0; i <= N; i++){
    for (int j = 0; j <= N; j++){
      a[i][j][0] = 0.;
      a[i][j][N] = 0.;
    }
  }

  for (int j = 0; j <= N; j++){
    for (int k = 0; k <= N; k++){
      a[0][j][k] = 0.;
      a[N][j][k] = 0.;
    }
  }

  for (int i = 1; i < N; i++){
    for (int j = 1; j <= N; j++){
      for (int k = 1; k < N; k++){
        a[i][j][k] = pa[i][j][k] + (g->tau * g->tau / 2.0) * (
          ((fi(g, i-1, j, k) - 2. * (fi(g, i, j, k)) + fi(g, i+1, j, k)) / (g->hx * g->hx)) +
          ((fi(g, i, j-1, k) - 2. * (fi(g, i, j, k)) + fi(g, i, j+1, k)) / (g->hy * g->hy)) +
          ((fi(g, i, j, k-1) - 2. * (fi(g, i, j, k)) + fi(g, i, j, k+1)) / (g->hz * g->hz)) 
        );
      }
    }
  }
  
  for (int i = 1; i < N; i++){
    for (int k = 1; k < N; k++){
      a[i][0][k] = a[i][N][k];
    }
  }

  return;
}

void u_step (double (*a)[N + 1][N + 1], double (*pa)[N + 1][N + 1], double (*ppa)[N + 1][N + 1], gin * g){

  for (int i = 0; i <= N; i++){
    for (int j = 0; j <= N; j++){
      a[i][j][0] = 0;
      a[i][j][N] = 0;
    }
  }

  for (int j = 0; j <= N; j++){
    for (int k = 0; k <= N; k++){
      a[0][j][k] = 0;
      a[N][j][k] = 0;
    }
  }

  for (int i = 1; i < N; i++){
    for (int k = 1; k < N; k++){ 
      for (int j = 1; j < N; j++){
        a[i][j][k] = 2 * pa[i][j][k] - ppa[i][j][k] + (g->tau * g->tau) * (
          ((pa[i-1][j][k] - 2 * pa[i][j][k] + pa[i+1][j][k]) / (g->hx * g->hx)) +
          ((pa[i][j-1][k] - 2 * pa[i][j][k] + pa[i][j+1][k]) / (g->hy * g->hy)) +
          ((pa[i][j][k-1] - 2 * pa[i][j][k] + pa[i][j][k+1]) / (g->hz * g->hz))
        );
      }
      a[i][N][k] = 2 * pa[i][N][k] - ppa[i][N][k] + (g->tau * g->tau) * (
          ((pa[i-1][N][k] - 2 * pa[i][N][k] + pa[i+1][N][k]) / (g->hx * g->hx)) +
          ((pa[i][N-1][k] - 2 * pa[i][N][k] + pa[i][1][k]) / (g->hy * g->hy)) +
          ((pa[i][N][k-1] - 2 * pa[i][N][k] + pa[i][N][k+1]) / (g->hz * g->hz))
        );
      a[i][0][k] = a[i][N][k];
    }
  }

  return;
}

double u_err(double (*a)[N + 1][N + 1], int t, gin * g){
  double max = 0;
  double tmp;
   for (int i = 0; i <= N; i++){
    for (int j = 0; j <= N; j++){
      for (int k = 0; k <= N; k++){
        //printf("\t[%d][%d][%d]\t%lf\n",i,j,k,u(g,i,j,k,t));
        tmp = fabs(a[i][j][k] - u(g,i,j,k,t));
        if (tmp > max)
          max = tmp;
      }
    }
  }   
  return max;
}

int main(){

  clock_t start = clock();

  MPI_Init(NULL, NULL);

  int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

  double (*A)[N + 1][N + 1];
  double (*B)[N + 1][N + 1];
  double (*C)[N + 1][N + 1];
  A = malloc((N + 1) * sizeof(*A));
  B = malloc((N + 1) * sizeof(*B));
  C = malloc((N + 1)* sizeof(*C));


  if (world_size != 1){
    printf("BAD PROCESS COUNT MUST BE 1\n");
    return 0;
  }


  gin g;
  g.lx = 1.0;
  g.ly = 1.0;
  g.lz = 1.0;
  g.T  = 1.0;
  g.hx = g.lx / (double)N;
  g.hy = g.ly / (double)N;
  g.hz = g.lz / (double)N;
  g.tau = g.T / (double)K;

  double err = 0;
  double max_err = 0;

  // INIT A, t = 0
    u_0_init(A, &g);
    err = u_err(A, 0, &g);
    if (err > max_err)
        max_err = err;
    //print_grid(A);
    printf("ERROR[%d]::%10.10lf;\n",0,err);
  
  // INIT B, t = 1 
    u_1_init(B, A, &g);
    err = u_err(B, 1, &g);
    if (err > max_err)
        max_err = err;
    //print_grid(B);
    printf("ERROR[%d]::%10.10lf;\n",1, err);

  for(int s = 2; s <= 20 ; s++){
    if (s % 3 == 2){
      u_step(C, B, A, &g);
      err = u_err(C, s, &g);
      if (err > max_err)
        max_err = err;
      printf("ERROR[%d]::%10.10lf;\n",s,err);
    } else if (s % 3 == 1){
      u_step(B, A, C, &g);
      err = u_err(B, s, &g);
      if (err > max_err)
        max_err = err;
      printf("ERROR[%d]::%10.10lf;\n",s,err);
    } else if (s % 3 == 0){
      u_step(A, C, B, &g);
      err = u_err(A, s, &g);
      if (err > max_err)
        max_err = err;
      printf("ERROR[%d]::%10.10lf;\n",s,err);
    }
  }

  clock_t end = clock();
  double time = ((double)(end - start)) / CLOCKS_PER_SEC;

  printf("TIME=%lf;\n", time);
  printf("MAX_ERROR=%10.10lf;\n",max_err);

  MPI_Finalize();

  return 0;
}

