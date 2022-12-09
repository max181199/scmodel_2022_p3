#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <mpi.h>

#define N 511
#define K 886
#define FALSE 0  
#define TRUE !(FALSE)
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

typedef struct Lim {
  int r; // rank
  int world_size;
  MPI_Comm comm; // communicator
  int x; // x-node-coordinate
  int y; // y-node-coordinate
  int z; // z-node-coordinate
  int x_0; // x-inner-point-start
  int x_1; // x-inner-point-end
  int y_0; // y-inner-point-start
  int y_1; // y-inner-point-end
  int z_0; // z-inner-point-start
  int z_1; // z-inner-point-end
  int h_x_r; // x-halo-right
  int h_x_l; // x-halo-left
  int h_y_r; // y-halo-right
  int h_y_l; // y-halo-left
  int h_z_r; // z-halo-right
  int h_z_l; // z-halo-left
  int r_x; // x-right-neigborhood
  int l_x; // x-left-neigborhood
  int r_y; // y-right-neigborhood
  int l_y; // y-left-neigborhood
  int l_z; // z-left-neigborhood
  int r_z; // z-right-neigborhood
  int size_x;
  int size_y;
  int size_z;
  int from_N_to_0;
  int from_0_to_N;
} lim;


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

void u_0_init (double *a, gin * g, lim * l){

  for (int i = 1 ; i <= l->x_1 - l->x_0 + 1; i++){
    for (int j = 1; j <= l->y_1 - l->y_0 + 1; j++){
      for (int k = 1; k <= l->z_1 - l->z_0 + 1; k++){
        int index = i * l->size_y * l->size_z + j * l->size_z + k;
        a[index] = fi(g, l->x_0 + i - 1, l->y_0 + j - 1, l->z_0 + k - 1);
      }
    }
  }

}

double u_err(double (*a), int t, gin * g, lim * l){
  double max = 0;
  double tmp;

  for (int i = 1 ; i <= l->x_1 - l->x_0 + 1; i++){
    for (int j = 1; j <= l->y_1 - l->y_0 + 1; j++){
      for (int k = 1; k <= l->z_1 - l->z_0 + 1; k++){
        int index = i * l->size_y * l->size_z + j * l->size_z + k;
        tmp = fabs(a[index] - u(g, l->x_0 + i - 1, l->y_0 + j - 1, l->z_0 + k - 1, t));
        if (tmp > max)
          max = tmp;
      }
    }
  }  

  return max;
}

double u_full_err(double (*a), int t, gin * g, lim * l){
  double max = 0;
  double max_r_x = 0;
  double max_l_x = 0;
  double max_r_y = 0;
  double max_l_y = 0;
  double max_r_z = 0;
  double max_l_z = 0;
  double tmp;

  if (l->r_x != -1)
    for (int j = 1; j <= l->y_1 - l->y_0 + 1; j++){
      for (int k = 1; k <= l->z_1 - l->z_0 + 1; k++){
        int index = (l->x_1 - l->x_0 + 2) * l->size_y * l->size_z + j * l->size_z + k;
        tmp = fabs(a[index] - u(g, l->x_0 + (l->x_1 - l->x_0 + 2) - 1, l->y_0 + j - 1, l->z_0 + k - 1, t));
        if (tmp > max){
          max = tmp;
        }
        if (tmp > max_r_x){
          max_r_x = tmp;
        }
      }
  } 

  //printf("ER_RX[%d, %d, %d]--- = %lf\n", l->x, l->y, l->z, max_r_x);

  if (l->l_x != -1)
    for (int j = 1; j <= l->y_1 - l->y_0 + 1; j++){
      for (int k = 1; k <= l->z_1 - l->z_0 + 1; k++){
        int index = 0 * l->size_y * l->size_z + j * l->size_z + k;
        tmp = fabs(a[index] - u(g, l->x_0 + 0 - 1, l->y_0 + j - 1, l->z_0 + k - 1, t));
        if (tmp > max)
          max = tmp;
        if (tmp > max_l_x){
          max_l_x = tmp;
        }
      }
  } 

  //printf("ER_LX[%d, %d, %d]--- = %lf\n", l->x, l->y, l->z, max_l_x);

  if(l->r_y != -1)
    for (int i = 1 ; i <= l->x_1 - l->x_0 + 1; i++){
      for (int k = 1; k <= l->z_1 - l->z_0 + 1; k++){
        int index = i * l->size_y * l->size_z + (l->y_1 - l->y_0 + 2) * l->size_z + k;
        tmp = fabs(a[index] - u(g, l->x_0 + i - 1, l->y_0 + (l->y_1 - l->y_0 + 2) - 1, l->z_0 + k - 1, t));
        if (tmp > max)
          max = tmp;
        if (tmp > max_r_y){
          max_r_y = tmp;
        }
      }
    }
  
   //printf("ER_RY[%d, %d, %d][%d][L:%d, R:%d]--- = %lf\n", l->x, l->y, l->z,l->r,l->l_y,l->r_y, max_r_y);
  

  if(l->l_y != -1)
    for (int i = 1 ; i <= l->x_1 - l->x_0 + 1; i++){
      for (int k = 1; k <= l->z_1 - l->z_0 + 1; k++){
        int index = i * l->size_y * l->size_z + 0 * l->size_z + k;
        tmp = fabs(a[index] - u(g, l->x_0 + i - 1, l->y_0 + 0 - 1, l->z_0 + k - 1, t));
        if (tmp > max)
          max = tmp;
        if (tmp > max_l_y){
          max_l_y = tmp;
        }
      }
    }

     //printf("ER_LY[%d, %d, %d][%d][L:%d, R:%d]--- = %lf\n", l->x, l->y, l->z,l->r,l->l_y,l->r_y, max_l_y);

  if(l->r_z != -1)
  for (int i = 1 ; i <= l->x_1 - l->x_0 + 1; i++){
    for (int j = 1; j <= l->y_1 - l->y_0 + 1; j++){
        int index = i * l->size_y * l->size_z + j * l->size_z + (l->z_1 - l->z_0 + 2);
        tmp = fabs(a[index] - u(g, l->x_0 + i - 1, l->y_0 + j - 1, l->z_0 + (l->z_1 - l->z_0 + 2) - 1, t));
        if (tmp > max)
          max = tmp;
        if (tmp > max_r_z){
          max_r_z = tmp;
        }
    }
  } 

  //printf("ER_RZ[%d, %d, %d]--- = %lf\n", l->x, l->y, l->z, max_r_z);

  if(l->l_z != -1)
  for (int i = 1 ; i <= l->x_1 - l->x_0 + 1; i++){
    for (int j = 1; j <= l->y_1 - l->y_0 + 1; j++){
        int index = i * l->size_y * l->size_z + j * l->size_z + 0;
        tmp = fabs(a[index] - u(g, l->x_0 + i - 1, l->y_0 + j - 1, l->z_0 + 0 - 1, t));
        if (tmp > max)
          max = tmp;
        if (tmp > max_l_z){
          max_l_z = tmp;
        }
    }
  }

  //printf("ER_LZ[%d, %d, %d]--- = %lf\n", l->x, l->y, l->z, max_l_z);

  return max;
}

void swap_halo (double (*a), lim * l, gin * g){
  // create requests pull
  MPI_Request req_x_to_left;
  MPI_Request req_y_to_left;
  MPI_Request req_z_to_left;
  MPI_Request req_x_to_right;
  MPI_Request req_y_to_right;
  MPI_Request req_z_to_right;
  MPI_Request req_x_from_left;
  MPI_Request req_y_from_left;
  MPI_Request req_z_from_left;
  MPI_Request req_x_from_right;
  MPI_Request req_y_from_right;
  MPI_Request req_z_from_right;
  // create buffers
  double (*send_x_to_left)[l->z_1 - l->z_0 + 1];
  double (*recive_x_from_right)[l->z_1 - l->z_0 + 1];
  double (*send_x_to_right)[l->z_1 - l->z_0 + 1];
  double (*recive_x_from_left)[l->z_1 - l->z_0 + 1];
  double (*send_y_to_left)[l->z_1 - l->z_0 + 1];
  double (*send_y_to_right)[l->z_1 - l->z_0 + 1];
  double (*recive_y_from_left)[l->z_1 - l->z_0 + 1];
  double (*recive_y_from_right)[l->z_1 - l->z_0 + 1];
  double (*send_z_to_left)[l->y_1 - l->y_0 + 1];
  double (*send_z_to_right)[l->y_1 - l->y_0 + 1];
  double (*recive_z_from_left)[l->y_1 - l->y_0 + 1];
  double (*recive_z_from_right)[l->y_1 - l->y_0 + 1];
  send_x_to_left = malloc((l->y_1 - l->y_0 + 1) * sizeof (*send_x_to_left));
  recive_x_from_right = malloc((l->y_1 - l->y_0 + 1) * sizeof (*recive_x_from_right));
  send_x_to_right = malloc((l->y_1 - l->y_0 + 1) * sizeof (*send_x_to_right));
  recive_x_from_left = malloc((l->y_1 - l->y_0 + 1) * sizeof (*recive_x_from_left));
  send_y_to_left = malloc((l->x_1 - l->x_0 + 1) * sizeof (*send_y_to_left));
  send_y_to_right = malloc((l->x_1 - l->x_0 + 1) * sizeof (*send_y_to_right));
  recive_y_from_left = malloc((l->x_1 - l->x_0 + 1) * sizeof (*recive_y_from_left));
  recive_y_from_right = malloc((l->x_1 - l->x_0 + 1) * sizeof (*recive_y_from_right));
  send_z_to_left = malloc((l->x_1 - l->x_0 + 1) * sizeof (*send_z_to_left));
  send_z_to_right = malloc((l->x_1 - l->x_0 + 1) * sizeof (*send_z_to_right));
  recive_z_from_left = malloc((l->x_1 - l->x_0 + 1) * sizeof (*recive_z_from_left));
  recive_z_from_right = malloc((l->x_1 - l->x_0 + 1) * sizeof (*recive_z_from_right));

  for (int j = 1; j <= l->y_1 - l->y_0 + 1; j++){
    for (int k = 1; k <= l->z_1 - l->z_0 + 1; k++){
      send_x_to_right[j - 1][k - 1] = 0;
      send_x_to_left[j - 1][k - 1] = 0;
    }
  }
   
  for (int i = 1; i <= l->x_1 - l->x_0 + 1; i++){
    for (int k = 1; k <= l->z_1 - l->z_0 + 1; k++){
      send_y_to_right[i - 1][k - 1] = 0;
      send_y_to_left[i - 1][k - 1] = 0;
    }
  }

  for (int i = 1; i <= l->x_1 - l->x_0 + 1; i++){
    for (int j = 1; j <= l->y_1 - l->y_0 + 1; j++){
      send_z_to_right[i - 1][j - 1] = 0;
      send_z_to_left[i - 1][j - 1] = 0;
    }
  }

  // fill buffers   
  if(l->r_x != -1)  
    for (int j = 1; j <= l->y_1 - l->y_0 + 1; j++){
      for (int k = 1; k <= l->z_1 - l->z_0 + 1; k++){
        int index = (l->x_1 - l->x_0 + 1) * l->size_y * l->size_z + j * l->size_z + k;
        send_x_to_right[j - 1][k - 1] = a[index];
      }
    }
  if(l->l_x != -1) 
    for (int j = 1; j <= l->y_1 - l->y_0 + 1; j++){
      for (int k = 1; k <= l->z_1 - l->z_0 + 1; k++){
        int index = 1 * l->size_y * l->size_z + j * l->size_z + k;
        send_x_to_left[j - 1][k - 1] = a[index];
      }
    }
  if(l->r_y != -1)     
    for (int i = 1; i <= l->x_1 - l->x_0 + 1; i++){
      for (int k = 1; k <= l->z_1 - l->z_0 + 1; k++){
        int index = i * l->size_y * l->size_z + (l->y_1 - l->y_0 + 1) * l->size_z + k;
        send_y_to_right[i - 1][k - 1] = a[index];
        // if ((N >= l->y_0) && (N <= l->y_1))
        // printf("[%d on %d to  %d ]::%lf\n", l->y_1, l->r, l->r_y, recive_y_from_left[i - 1][k - 1]);
      }
    }
  if(l->l_y != -1) 
    for (int i = 1; i <= l->x_1 - l->x_0 + 1; i++){
      for (int k = 1; k <= l->z_1 - l->z_0 + 1; k++){
        int index = i * l->size_y * l->size_z + (1) * l->size_z + k;
        send_y_to_left[i - 1][k - 1] = a[index];
      }
    }
  if(l->r_z != -1) 
    for (int i = 1; i <= l->x_1 - l->x_0 + 1; i++){
      for (int j = 1; j <= l->y_1 - l->y_0 + 1; j++){
        int index = i * l->size_y * l->size_z + j * l->size_z + (l->z_1 - l->z_0 + 1);
        send_z_to_right[i - 1][j - 1] = a[index];
      }
    }
  if(l->l_z != -1) 
    for (int i = 1; i <= l->x_1 - l->x_0 + 1; i++){
      for (int j = 1; j <= l->y_1 - l->y_0 + 1; j++){
        int index = i * l->size_y * l->size_z + j * l->size_z + 1;
        send_z_to_left[i - 1][j - 1] = a[index];
      }
    }

  // SEND RECIVE X-R
  if (l->x % 2 == 0){
    if(l->r_x != -1){
      MPI_Send(send_x_to_right, (l->y_1 - l->y_0 + 1) * (l->z_1 - l->z_0 + 1), MPI_DOUBLE, l->r_x, 0, l->comm);
      MPI_Recv(recive_x_from_right, (l->y_1 - l->y_0 + 1) * (l->z_1 - l->z_0 + 1), MPI_DOUBLE, l->r_x, 0, l->comm, MPI_STATUS_IGNORE);
    }  
  } else {
    if(l->l_x != -1){
      MPI_Recv(recive_x_from_left, (l->y_1 - l->y_0 + 1) * (l->z_1 - l->z_0 + 1), MPI_DOUBLE, l->l_x, 0, l->comm, MPI_STATUS_IGNORE);
      MPI_Send(send_x_to_left, (l->y_1 - l->y_0 + 1) * (l->z_1 - l->z_0 + 1), MPI_DOUBLE, l->l_x, 0, l->comm);
    }  
  }
  // SEND RECIVE X-R
  
  // SEND RECIVE X-L
  if (l->x % 2 == 1){
    if(l->r_x != -1){
      MPI_Send(send_x_to_right, (l->y_1 - l->y_0 + 1) * (l->z_1 - l->z_0 + 1), MPI_DOUBLE, l->r_x, 0, l->comm);
      MPI_Recv(recive_x_from_right, (l->y_1 - l->y_0 + 1) * (l->z_1 - l->z_0 + 1), MPI_DOUBLE, l->r_x, 0, l->comm, MPI_STATUS_IGNORE);
    }  
  } else {
    if(l->l_x != -1){
      MPI_Recv(recive_x_from_left, (l->y_1 - l->y_0 + 1) * (l->z_1 - l->z_0 + 1), MPI_DOUBLE, l->l_x, 0, l->comm, MPI_STATUS_IGNORE);
      MPI_Send(send_x_to_left, (l->y_1 - l->y_0 + 1) * (l->z_1 - l->z_0 + 1), MPI_DOUBLE, l->l_x, 0, l->comm);
    }  
  }
  // // SEND RECIVE X-L
  
  // SEND RECIVE Y-R
  if (l->y % 2 == 0){
    if(l->r_y != -1){
      MPI_Send(send_y_to_right, (l->x_1 - l->x_0 + 1) * (l->z_1 - l->z_0 + 1), MPI_DOUBLE, l->r_y, 0, l->comm);
      MPI_Recv(recive_y_from_right, (l->x_1 - l->x_0 + 1) * (l->z_1 - l->z_0 + 1), MPI_DOUBLE, l->r_y, 0, l->comm, MPI_STATUS_IGNORE);
    } 
  } else {
    if(l->l_y != -1){
      MPI_Recv(recive_y_from_left, (l->x_1 - l->x_0 + 1) * (l->z_1 - l->z_0 + 1), MPI_DOUBLE, l->l_y, 0, l->comm, MPI_STATUS_IGNORE);
      MPI_Send(send_y_to_left, (l->x_1 - l->x_0 + 1) * (l->z_1 - l->z_0 + 1), MPI_DOUBLE, l->l_y, 0, l->comm);
    }  
  }
  // SEND RECIVE Y-R

  // SEND RECIVE Y-L
  if (l->y % 2 == 1){
    if(l->r_y != -1) {
      MPI_Send(send_y_to_right, (l->x_1 - l->x_0 + 1) * (l->z_1 - l->z_0 + 1), MPI_DOUBLE, l->r_y, 0, l->comm);
      MPI_Recv(recive_y_from_right, (l->x_1 - l->x_0 + 1) * (l->z_1 - l->z_0 + 1), MPI_DOUBLE, l->r_y, 0, l->comm, MPI_STATUS_IGNORE);
    }
  } else {
    if(l->l_y != -1){
      MPI_Recv(recive_y_from_left, (l->x_1 - l->x_0 + 1) * (l->z_1 - l->z_0 + 1), MPI_DOUBLE, l->l_y, 0, l->comm, MPI_STATUS_IGNORE);
      MPI_Send(send_y_to_left, (l->x_1 - l->x_0 + 1) * (l->z_1 - l->z_0 + 1), MPI_DOUBLE, l->l_y, 0, l->comm);
    }
  }
  // SEND RECIVE Y-L

  // SEND RECIVE Z-R
  if (l->z % 2 == 0){
    if(l->r_z != -1){
      MPI_Send(send_z_to_right, (l->x_1 - l->x_0 + 1) * (l->y_1 - l->y_0 + 1), MPI_DOUBLE, l->r_z, 0, l->comm);
      MPI_Recv(recive_z_from_right, (l->x_1 - l->x_0 + 1) * (l->y_1 - l->y_0 + 1), MPI_DOUBLE, l->r_z, 0, l->comm, MPI_STATUS_IGNORE);
    }  
  } else {
    if(l->l_z != -1){
      MPI_Recv(recive_z_from_left, (l->x_1 - l->x_0 + 1) * (l->y_1 - l->y_0 + 1), MPI_DOUBLE, l->l_z, 0, l->comm, MPI_STATUS_IGNORE);
      MPI_Send(send_z_to_left, (l->x_1 - l->x_0 + 1) * (l->y_1 - l->y_0 + 1), MPI_DOUBLE, l->l_z, 0, l->comm);
    }  
  }
  // SEND RECIVE Z-R

  // SEND RECIVE Z-L
  if (l->z % 2 == 1){
    if(l->r_z != -1){
      MPI_Send(send_z_to_right, (l->x_1 - l->x_0 + 1) * (l->y_1 - l->y_0 + 1), MPI_DOUBLE, l->r_z, 0, l->comm);
      MPI_Recv(recive_z_from_right, (l->x_1 - l->x_0 + 1) * (l->y_1 - l->y_0 + 1), MPI_DOUBLE, l->r_z, 0, l->comm, MPI_STATUS_IGNORE);
    }  
  } else {
    if(l->l_z != -1){
      MPI_Recv(recive_z_from_left, (l->x_1 - l->x_0 + 1) * (l->y_1 - l->y_0 + 1), MPI_DOUBLE, l->l_z, 0, l->comm, MPI_STATUS_IGNORE);
      MPI_Send(send_z_to_left, (l->x_1 - l->x_0 + 1) * (l->y_1 - l->y_0 + 1), MPI_DOUBLE, l->l_z, 0, l->comm);
    }  
  }
  // SEND RECIVE Z-L

  //send buffers without waiting
  // if(l->r_x != -1)  
  //   MPI_Isend(send_x_to_right, (l->y_1 - l->y_0 + 1) * (l->z_1 - l->z_0 + 1), MPI_DOUBLE, l->r_x, 1, l->comm, &req_x_to_right);
  // if(l->l_x != -1)  
  //   MPI_Isend(send_x_to_left, (l->y_1 - l->y_0 + 1) * (l->z_1 - l->z_0 + 1), MPI_DOUBLE, l->l_x, 2, l->comm, &req_x_to_left);
  // if(l->r_y != -1)
  //   MPI_Isend(send_y_to_right, (l->x_1 - l->x_0 + 1) * (l->z_1 - l->z_0 + 1), MPI_DOUBLE, l->r_y, 3, l->comm, &req_y_to_right);
  // if(l->l_y != -1)
  //   MPI_Isend(send_y_to_left, (l->x_1 - l->x_0 + 1) * (l->z_1 - l->z_0 + 1), MPI_DOUBLE, l->l_y, 4, l->comm, &req_y_to_left);
  // if(l->r_z != -1)
  //   MPI_Isend(send_z_to_right, (l->y_1 - l->y_0 + 1) * (l->x_1 - l->x_0 + 1), MPI_DOUBLE, l->r_z, 5, l->comm, &req_z_to_right);
  // if(l->l_z != -1)
  //   MPI_Isend(send_z_to_left, (l->y_1 - l->y_0 + 1) * (l->x_1 - l->x_0 + 1), MPI_DOUBLE, l->l_z, 6, l->comm, &req_z_to_left);
  // //recive buffers withoud waiting
  // if(l->r_x != -1)  
  //   MPI_Irecv(recive_x_from_right, (l->y_1 - l->y_0 + 1) * (l->z_1 - l->z_0 + 1), MPI_DOUBLE, l->r_x, 2, l->comm, &req_x_from_right);
  // if(l->l_x != -1)  
  //   MPI_Irecv(recive_x_from_left, (l->y_1 - l->y_0 + 1) * (l->z_1 - l->z_0 + 1), MPI_DOUBLE, l->l_x, 1, l->comm, &req_x_from_left);
  // if(l->r_y != -1)
  //   MPI_Irecv(recive_y_from_right, (l->x_1 - l->x_0 + 1) * (l->z_1 - l->z_0 + 1), MPI_DOUBLE, l->r_y, 4, l->comm, &req_y_from_right);
  // if(l->l_y != -1)
  //   MPI_Irecv(recive_y_from_left, (l->x_1 - l->x_0 + 1) * (l->z_1 - l->z_0 + 1), MPI_DOUBLE, l->l_y, 3, l->comm, &req_y_from_left);
  // if(l->r_z != -1)
  //   MPI_Irecv(recive_z_from_right, (l->y_1 - l->y_0 + 1) * (l->x_1 - l->x_0 + 1), MPI_DOUBLE, l->r_z, 6, l->comm, &req_z_from_right);
  // if(l->l_z != -1)
  //   MPI_Irecv(recive_z_from_left, (l->y_1 - l->y_0 + 1) * (l->x_1 - l->x_0 + 1), MPI_DOUBLE, l->l_z, 5, l->comm, &req_z_from_left);
  // //waiting all 
  // if(l->r_x != -1)  
  //   MPI_Wait(&req_x_from_right, MPI_STATUS_IGNORE);
  // if(l->l_x != -1)
  //   MPI_Wait(&req_x_from_left, MPI_STATUS_IGNORE);
  // if(l->r_y != -1)
  //   MPI_Wait(&req_y_from_right, MPI_STATUS_IGNORE);
  // if(l->l_y != -1)
  //   MPI_Wait(&req_y_from_left, MPI_STATUS_IGNORE);
  // if(l->r_z != -1)
  //   MPI_Wait(&req_z_from_right, MPI_STATUS_IGNORE); 
  // if(l->l_z != -1)
  //   MPI_Wait(&req_z_from_left, MPI_STATUS_IGNORE);
  //change a matrix
  if(l->r_x != -1)  
    for (int j = 1; j <= l->y_1 - l->y_0 + 1; j++){
      for (int k = 1; k <= l->z_1 - l->z_0 + 1; k++){
        int index = (l->x_1 - l->x_0 + 2) * l->size_y * l->size_z + j * l->size_z + k;
        a[index] = recive_x_from_right[j - 1][k - 1];
        //printf("{%d, %d, %d} = %lf\n", (l->x_1 - l->x_0 + 2), j, k, recive_x_from_right[j - 1][k - 1] - u(g, l->x_1 + 1, l->y_0 + j - 1, l->z_0 + k - 1, 0));
      }
    }
  if(l->l_x != -1)
    for (int j = 1; j <= l->y_1 - l->y_0 + 1; j++){
      for (int k = 1; k <= l->z_1 - l->z_0 + 1; k++){
        int index = (0) * l->size_y * l->size_z + j * l->size_z + k;
        a[index] = recive_x_from_left[j - 1][k - 1];
        //printf("{%d, %d, %d} = %lf\n", 0, j, k, recive_x_from_left[j - 1][k - 1] - u(g, l->x_0 - 1, l->y_0 + j - 1, l->z_0 + k - 1, 0));
      }
    }
  if(l->r_y != -1)
    for (int i = 1; i <= l->x_1 - l->x_0 + 1; i++){
      for (int k = 1; k <= l->z_1 - l->z_0 + 1; k++){
        int index = i * l->size_y * l->size_z + (l->y_1 - l->y_0 + 2) * l->size_z + k;
        a[index] = recive_y_from_right[i - 1][k - 1];
      }
    }
  if(l->l_y != -1)
    for (int i = 1; i <= l->x_1 - l->x_0 + 1; i++){
      for (int k = 1; k <= l->z_1 - l->z_0 + 1; k++){
        int index = i * l->size_y * l->size_z + 0 * l->size_z + k;
        a[index] = recive_y_from_left[i - 1][k - 1];
      }
    }
  if(l->r_z != -1)
    for (int i = 1; i <= l->x_1 - l->x_0 + 1; i++){
      for (int j = 1; j <= l->y_1 - l->y_0 + 1; j++){
        int index = i * l->size_y * l->size_z + j * l->size_z + (l->z_1 - l->z_0 + 2);
        a[index] = recive_z_from_right[i - 1][j - 1];
      }
    }
  if(l->l_z != -1)
    for (int i = 1; i <= l->x_1 - l->x_0 + 1; i++){
      for (int j = 1; j <= l->y_1 - l->y_0 + 1; j++){
        int index = i * l->size_y * l->size_z + j * l->size_z + 0;
        a[index] = recive_z_from_left[i - 1][j - 1];
      }
    }

  free(send_x_to_left);
  free(recive_x_from_right);
  free(send_x_to_right);
  free(recive_x_from_left);
  free(send_y_to_left);
  free(send_y_to_right);
  free(recive_y_from_left);
  free(recive_y_from_right);
  free(send_z_to_left);
  free(send_z_to_right);
  free(recive_z_from_left);
  free(recive_z_from_right);

  return;
}

void after_swap_copy_j0_to_j1(double * a, lim * l){
  if ((0 >= l->y_0) && (0 <= l->y_1)){
    for (int i = 1; i <= l->x_1 - l->x_0 + 1; i++){
      for (int k = 1; k <= l->z_1 - l->z_0 + 1; k++){
        int index_a = i * l->size_y * l->size_z + 0 * l->size_z + k;
        int index_b = i * l->size_y * l->size_z + 1 * l->size_z + k;
        //printf("[%d]::{%lf}||{%lf}\n",l->y_0, a[index_b], a[index_a]);
        a[index_b] = a[index_a];
      }
    }
  }
  return;
}


void u_1_init (double (*a), double (*pa), gin * g, lim * l){
  
  // for (int i = 1 ; i <= l->x_1 - l->x_0 + 1; i++){
  //   for (int j = 1; j <= l->y_1 - l->y_0 + 1; j++){
  //     for (int k = 1; k <= l->z_1 - l->z_0 + 1; k++){
  //       int index = i * l->size_y * l->size_z + j * l->size_z + k;
  //       a[index] = u(g, l->x_0 + i - 1, l->y_0 + j - 1, l->z_0 + k - 1, 1);
  //     }
  //   }
  // }
  
  if ((l->z_0 <= 0) && (0 <= l->z_1))
    for (int i = 1; i <= l->x_1 - l->x_0 + 1; i++){
      for (int j = 1; j <= l->y_1 - l->y_0 + 1; j++){
        int index = i * l->size_y * l->size_z + j * l->size_z + 1;
        a[index] = 0.;
      }
    }
  if((l->z_0 <= N) && (N <= l->z_1))
    for (int i = 1; i <= l->x_1 - l->x_0 + 1; i++){
      for (int j = 1; j <= l->y_1 - l->y_0 + 1; j++){
        int index = i * l->size_y * l->size_z + j * l->size_z + (l->z_1 - l->z_0 + 1);
        a[index] = 0.;
      }
    }


  if ((l->x_0 <= 0) && (0 <= l->x_1))
    for (int j = 1; j <= l->y_1 - l->y_0 + 1; j++){
      for (int k = 1; k <= l->z_1 - l->z_0 + 1; k++){
        int index = 0 * l->size_y * l->size_z + j * l->size_z + k;
        a[index] = 0.;
      }
    }
  if((l->x_0 <= N) && (N <= l->x_1))
    for (int j = 1; j <= l->y_1 - l->y_0 + 1; j++){
      for (int k = 1; k <= l->z_1 - l->z_0 + 1; k++){
        int index = (l->x_1 - l->x_0 + 1) * l->size_y * l->size_z + j * l->size_z + k;
        a[index] = 0.;
      }
    }

  for (int i = 1 ; i <= l->x_1 - l->x_0 + 1; i++){
    for (int j = 1; j <= l->y_1 - l->y_0 + 1; j++){
      for (int k = 1; k <= l->z_1 - l->z_0 + 1; k++){
        int index = i * l->size_y * l->size_z + j * l->size_z + k;
        int il = l->x_0 + i - 1;
        int jl = l->y_0 + j - 1;
        int kl = l->z_0 + k - 1;
        a[index] = pa[index] + (g->tau * g->tau / 2.0) * (
          ((fi(g, il - 1, jl, kl) - 2. * (fi(g, il, jl, kl)) + fi(g, il + 1, jl, kl)) / (g->hx * g->hx)) +
          ((fi(g, il, jl - 1, kl) - 2. * (fi(g, il, jl, kl)) + fi(g, il, jl + 1, kl)) / (g->hy * g->hy)) +
          ((fi(g, il, jl, kl - 1) - 2. * (fi(g, il, jl, kl)) + fi(g, il, jl, kl + 1)) / (g->hz * g->hz)) 
        );
      }
    }
  }

  return;
}

void u_step (double (*a), double (*pa), double (*ppa), gin * g, lim * l, int s ){
  
  // for (int i = 1 ; i <= l->x_1 - l->x_0 + 1; i++){
  //   for (int j = 1; j <= l->y_1 - l->y_0 + 1; j++){
  //     for (int k = 1; k <= l->z_1 - l->z_0 + 1; k++){
  //       int index = i * l->size_y * l->size_z + j * l->size_z + k;
  //       a[index] = u(g, l->x_0 + i - 1, l->y_0 + j - 1, l->z_0 + k - 1, s);
  //     }
  //   }
  // }

  if ((l->z_0 <= 0) && (0 <= l->z_1))
    for (int i = 1; i <= l->x_1 - l->x_0 + 1; i++){
      for (int j = 1; j <= l->y_1 - l->y_0 + 1; j++){
        int index = i * l->size_y * l->size_z + j * l->size_z + 1;
        a[index] = 0.;
      }
    }
  if((l->z_0 <= N) && (N <= l->z_1))
    for (int i = 1; i <= l->x_1 - l->x_0 + 1; i++){
      for (int j = 1; j <= l->y_1 - l->y_0 + 1; j++){
        int index = i * l->size_y * l->size_z + j * l->size_z + (l->z_1 - l->z_0 + 1);
        a[index] = 0.;
      }
    }


  if ((l->x_0 <= 0) && (0 <= l->x_1))
    for (int j = 1; j <= l->y_1 - l->y_0 + 1; j++){
      for (int k = 1; k <= l->z_1 - l->z_0 + 1; k++){
        int index = 0 * l->size_y * l->size_z + j * l->size_z + k;
        a[index] = 0.;
      }
    }
  if((l->x_0 <= N) && (N <= l->x_1))
    for (int j = 1; j <= l->y_1 - l->y_0 + 1; j++){
      for (int k = 1; k <= l->z_1 - l->z_0 + 1; k++){
        int index = (l->x_1 - l->x_0 + 1) * l->size_y * l->size_z + j * l->size_z + k;
        a[index] = 0.;
      }
    }


  int i_start = l->x_0 == 0 ? 2 : 1;
  int i_stop =  l->x_1 == N ? l->x_1 - l->x_0 : l->x_1 - l->x_0 + 1;
  int j_start = l->y_0 == 0 ? 2 : 1;
  int j_stop =  l->y_1 == N ? l->y_1 - l->y_0 : l->y_1 - l->y_0 + 1;
  int z_start = l->z_0 == 0 ? 2 : 1;
  int z_stop =  l->z_1 == N ? l->z_1 - l->z_0 : l->z_1 - l->z_0 + 1;

  for (int i = i_start; i <= i_stop; i++){
    for (int j = j_start; j <= j_stop; j++){ 
      for (int k = z_start; k <= z_stop; k++){
        int index = i * l->size_y * l->size_z + j * l->size_z + k;
        int index_xl = (i - 1) * l->size_y * l->size_z + j * l->size_z + k;
        int index_xr = (i + 1) * l->size_y * l->size_z + j * l->size_z + k;
        int index_yl = i * l->size_y * l->size_z + (j - 1) * l->size_z + k;
        int index_yr = i * l->size_y * l->size_z + (j + 1) * l->size_z + k;
        int index_zl = i * l->size_y * l->size_z + j * l->size_z + (k - 1);
        int index_zr = i * l->size_y * l->size_z + j * l->size_z + (k + 1);
        a[index] = 2.0 * pa[index] - ppa[index] + (g->tau * g->tau) * (
          ((pa[index_xl] - 2. * pa[index] + pa[index_xr]) / (g->hx * g->hx)) +
          ((pa[index_yl] - 2. * pa[index] + pa[index_yr]) / (g->hy * g->hy)) +
          ((pa[index_zl] - 2. * pa[index] + pa[index_zr]) / (g->hz * g->hz))
        );
      }
    }
  }

  if (l->world_size == 1 || l->from_0_to_N == l->r || l->from_N_to_0 == l->r){
    for (int i = 1; i <= l->x_1 - l->x_0 + 1; i++){
        for (int k = 1; k <= l->z_1 - l->z_0 + 1; k++){
          int index = i * l->size_y * l->size_z + (l->y_1 - l->y_0 + 1) * l->size_z + k;
          int index_xl = (i - 1) * l->size_y * l->size_z + (l->y_1 - l->y_0 + 1) * l->size_z + k;
          int index_xr = (i + 1) * l->size_y * l->size_z + (l->y_1 - l->y_0 + 1) * l->size_z + k;
          int index_yl = i * l->size_y * l->size_z + ((l->y_1 - l->y_0 + 1) - 1) * l->size_z + k;
          int index_zl = i * l->size_y * l->size_z + (l->y_1 - l->y_0 + 1) * l->size_z + (k - 1);
          int index_zr = i * l->size_y * l->size_z + (l->y_1 - l->y_0 + 1) * l->size_z + (k + 1);
          int index_s = i * l->size_y * l->size_z + (2) * l->size_z + k;
          a[index] = 2 * pa[index] - ppa[index] + (g->tau * g->tau) * (
            ((pa[index_xl] - 2 * pa[index] + pa[index_xr]) / (g->hx * g->hx)) +
            ((pa[index_yl] - 2 * pa[index] + pa[index_s]) / (g->hy * g->hy)) +
            ((pa[index_zl] - 2 * pa[index] + pa[index_zr]) / (g->hz * g->hz))
          );
        }
      }  
  } else {
    if (l->from_0_to_N != -1){
      double (*send)[l->z_1 - l->z_0 + 1];
      double (*recive)[l->z_1 - l->z_0 + 1];
      send = malloc((l->x_1 - l->x_0 + 1) * sizeof(*send));
      recive = malloc((l->x_1 - l->x_0 + 1) * sizeof(*recive));
      for (int i = 1; i <= l->x_1 - l->x_0 + 1; i++){
        for (int k = 1; k <= l->z_1 - l->z_0 + 1; k++){
          int index = i * l->size_y * l->size_z + (2) * l->size_z + k;
          send[i - 1][k - 1] = pa[index];
        }
      } 
      MPI_Send(send, (l->x_1 - l->x_0 + 1) * (l->z_1 - l->z_0 + 1), MPI_DOUBLE, l->from_0_to_N, l->from_0_to_N,  l->comm);  
      MPI_Recv(recive, (l->x_1 - l->x_0 + 1) * (l->z_1 - l->z_0 + 1), MPI_DOUBLE, l->from_0_to_N, l->r,  l->comm, MPI_STATUS_IGNORE);
      for (int i = 1; i <= l->x_1 - l->x_0 + 1; i++){
        for (int k = 1; k <= l->z_1 - l->z_0 + 1; k++){
          int index = i * l->size_y * l->size_z + (1) * l->size_z + k;
          pa[index] = recive[i - 1][k - 1];
        }
      } 
      free(send);
    }

    if (l->from_N_to_0 != -1){
      double (*recive)[l->z_1 - l->z_0 + 1];
      double (*send)[l->z_1 - l->z_0 + 1];
      recive = malloc((l->x_1 - l->x_0 + 1) * sizeof(*recive));
      send = malloc((l->x_1 - l->x_0 + 1) * sizeof(*send));
      for (int i = 1; i <= l->x_1 - l->x_0 + 1; i++){
        for (int k = 1; k <= l->z_1 - l->z_0 + 1; k++){
          int index = i * l->size_y * l->size_z + (l->y_1 - l->y_0 + 1) * l->size_z + k;
          send[i - 1][k - 1] = pa[index];
        }
      } 
      MPI_Recv(recive, (l->x_1 - l->x_0 + 1) * (l->z_1 - l->z_0 + 1), MPI_DOUBLE, l->from_N_to_0, l->r,  l->comm, MPI_STATUS_IGNORE);
      MPI_Send(send, (l->x_1 - l->x_0 + 1) * (l->z_1 - l->z_0 + 1), MPI_DOUBLE, l->from_N_to_0, l->from_N_to_0,  l->comm);  
      for (int i = i_start; i <= i_stop; i++){
        for (int k = z_start; k <= z_stop; k++){
          int index = i * l->size_y * l->size_z + (l->y_1 - l->y_0 + 1) * l->size_z + k;
          int index_xl = (i - 1) * l->size_y * l->size_z + (l->y_1 - l->y_0 + 1) * l->size_z + k;
          int index_xr = (i + 1) * l->size_y * l->size_z + (l->y_1 - l->y_0 + 1) * l->size_z + k;
          int index_yl = i * l->size_y * l->size_z + ((l->y_1 - l->y_0 + 1) - 1) * l->size_z + k;
          int index_zl = i * l->size_y * l->size_z + (l->y_1 - l->y_0 + 1) * l->size_z + (k - 1);
          int index_zr = i * l->size_y * l->size_z + (l->y_1 - l->y_0 + 1) * l->size_z + (k + 1);
          a[index] = 2 * pa[index] - ppa[index] + (g->tau * g->tau) * (
            ((pa[index_xl] - 2 * pa[index] + pa[index_xr]) / (g->hx * g->hx)) +
            ((pa[index_yl] - 2 * pa[index] + recive[i - 1][k - 1]) / (g->hy * g->hy)) +
            ((pa[index_zl] - 2 * pa[index] + pa[index_zr]) / (g->hz * g->hz))
          );
        }
      }   
      free(recive);
    }
  }
 


  return;
}


int main(){

  clock_t start = clock();

  MPI_Init(NULL, NULL);
  int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

  int dim[3];
  dim[0] = 0;
  dim[1] = 0;
  dim[2] = 0;

  MPI_Dims_create(world_size, 3, dim);

  if (dim[0] * dim[1] * dim[2] > world_size){
    printf("MPI_Dims_create_error\n");  
    return -1;
  }

  MPI_Comm MPI_COMM_CART;
  
  int period[3];
  period[1] = TRUE; period[2] = FALSE; period[0] = FALSE;
  int reorder;
  reorder = FALSE;
  int coords[3];
  MPI_Cart_create(MPI_COMM_WORLD, 3, dim, period, reorder, &MPI_COMM_CART);
  MPI_Cart_coords(MPI_COMM_CART, world_rank, 3, coords);
  //printf("DIM[%d]:{%d, %d, %d}\n", world_rank, dim[0], dim[1], dim[2]);
  lim l;
  l.world_size = world_size;
  l.r = world_rank;
  l.comm = MPI_COMM_CART;
  l.x = coords[0];
  l.y = coords[1];
  l.z = coords[2];
  //printf("COORDS[%d]:{%d, %d, %d}\n", world_rank, l.x, l.y, l.z);
  l.x_0 = coords[0] * ((N + 1) / dim[0]) + (coords[0] < (N + 1) % dim[0] ? coords[0] : (N + 1) % dim[0]);
  l.x_1 = l.x_0 + (N + 1) / dim[0] - 1 + (coords[0] < (N + 1) % dim[0] ?  1 : 0);
  l.y_0 = coords[1] * ((N + 1) / dim[1]) + (coords[1] < (N + 1) % dim[1] ? coords[1] : (N + 1) % dim[1]);
  l.y_1 = l.y_0 + (N + 1) / dim[1] - 1 + (coords[1] < (N + 1) % dim[1] ?  1 : 0);      
  l.z_0 = coords[2] * ((N + 1) / dim[2]) + (coords[2] < (N + 1) % dim[2] ? coords[2] : (N + 1) % dim[2]);
  l.z_1 = l.z_0 + (N + 1) / dim[2] - 1 + (coords[2] < (N + 1) % dim[2] ?  1 : 0);   
  //printf("LIMITATION[%d]:x{%d, %d}, y{%d, %d}, z{%d, %d} :: \n", world_rank, l.x_0, l.x_1, l.y_0, l.y_1, l.z_0, l.z_1);
  l.h_x_l = l.x_0 - 1;
  l.h_x_r = l.x_1 + 1;
  l.h_y_l = l.y_0 - 1;
  l.h_y_r = l.y_1 + 1;
  l.h_z_l = l.z_0 - 1;
  l.h_z_r = l.z_1 + 1;
  //printf("HALO[%d]:x{%d, %d}, y{%d, %d}, z{%d, %d}\n", world_rank, l.h_x_l, l.h_x_r, l.h_y_l, l.h_y_r, l.h_z_l, l.h_z_r);
  int src, dst;
  MPI_Cart_shift(MPI_COMM_CART, 0, 1, &src, &dst);
  l.r_x = dst;
  MPI_Cart_shift(MPI_COMM_CART, 0, -1, &src, &dst);
  l.l_x = dst;
  MPI_Cart_shift(MPI_COMM_CART, 1, 1, &src, &dst);
  l.r_y = l.y_1 == N ? -1 : dst;
  l.from_N_to_0 = l.y_1 != N ? -1 : dst;
  MPI_Cart_shift(MPI_COMM_CART, 1, -1, &src, &dst);
  l.l_y = l.y_0 == 0 ? -1 : dst;
  l.from_0_to_N = l.y_0 != 0 ? -1 : dst;
  MPI_Cart_shift(MPI_COMM_CART, 2, 1, &src, &dst);
  l.r_z = dst;
  MPI_Cart_shift(MPI_COMM_CART, 2, -1, &src, &dst);
  l.l_z = dst;
  //printf("NEIGBORHOODS[%d]:x{%d, %d}, y{%d, %d}, z{%d, %d}\n", world_rank, l.l_x, l.r_x, l.l_y, l.r_y, l.l_z, l.r_z);
  
  // Matrix Map
  //  0       1      ...     N        N + 1     
  //  h_x_l   x_0    ...     x_1      h_x_r
  //

  double (*A);
  double (*B);
  double (*C);
  A = malloc((l.x_1 - l.x_0 + 3) * (l.y_1 - l.y_0 + 3) * (l.z_1 - l.z_0 + 3) * sizeof(double));
  B = malloc((l.x_1 - l.x_0 + 3) * (l.y_1 - l.y_0 + 3) * (l.z_1 - l.z_0 + 3) * sizeof(double));
  C = malloc((l.x_1 - l.x_0 + 3) * (l.y_1 - l.y_0 + 3) * (l.z_1 - l.z_0 + 3) * sizeof(double));

  l.size_x = l.x_1 - l.x_0 + 3;
  l.size_y = l.y_1 - l.y_0 + 3;
  l.size_z = l.z_1 - l.z_0 + 3;

  for (int i = 0 ; i <= l.x_1 - l.x_0 + 2; i++){
    for (int j = 0; j <= l.y_1 - l.y_0 + 2; j++){ 
      for (int k = 0; k <= l.z_1 - l.z_0 + 2; k++){
        int index = i * l.size_y * l.size_z + j * l.size_z + k;
        A[index] = 0;
        B[index] = 0;
        C[index] = 0;
      }
    }
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
  double comm_err = 0;

  MPI_Barrier(MPI_COMM_CART);

  u_0_init(A, &g, &l);
  swap_halo(A, &l, &g);
  // err = u_full_err(A, 0, &g, &l);
  // printf("FE[%d, %d, %d]::%lf\n",l.x,l.y,l.z,err);
  after_swap_copy_j0_to_j1(A, &l);
  err = u_err(A, 0, &g, &l);
  MPI_Reduce(&err, &comm_err, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_CART);
  if (world_rank == 0)
    printf("ERROR[%d]::%10.10lf;\n",0,comm_err);

  MPI_Barrier(MPI_COMM_CART);

  u_1_init(B, A, &g, &l);
  swap_halo(B, &l,  &g);
  //  err = u_full_err(B, 1, &g, &l);
  //  printf("FE[%d, %d, %d]::%lf\n",l.x,l.y,l.z,err);
  after_swap_copy_j0_to_j1(B, &l);
  err = u_err(B, 1, &g, &l);
  MPI_Reduce(&err, &comm_err, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_CART);
  if (world_rank == 0)
    printf("ERROR[%d]::%10.10lf;\n",1,comm_err);

  MPI_Barrier(MPI_COMM_CART);

  // u_step(C, B, A, &g, &l, 2);
  // swap_halo(C, &l, &g);
  // after_swap_copy_j0_to_j1(C, &l);
  // err = u_err(C, 2, &g, &l);
  // MPI_Reduce(&err, &comm_err, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_CART);
  // if (world_rank == 0)
  //   printf("ERROR[%d]::%10.10lf;\n",2,comm_err);


  for(int s = 2; s <= 20 ; s++){
    if (s % 3 == 2){
      u_step(C, B, A, &g, &l, s);
      swap_halo(C, &l, &g);
      after_swap_copy_j0_to_j1(C, &l);
      //err = u_full_err(C, s, &g, &l);
      //printf("FE[%d, %d, %d][%d]::%lf\n",l.x,l.y,l.z,s,err);
      err = u_err(C, s, &g, &l);
      MPI_Reduce(&err, &comm_err, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_CART);
      if (world_rank == 0)
        printf("ERROR[%d]::%10.10lf;\n",s,comm_err);
    } else if (s % 3 == 1){
      u_step(B, A, C, &g, &l, s);
      swap_halo(B, &l, &g);
      after_swap_copy_j0_to_j1(B, &l);
      // err = u_full_err(B, s, &g, &l);
      // printf("FE[%d, %d, %d][%d]::%lf\n",l.x,l.y,l.z,s,err);
      err = u_err(B, s, &g, &l);
      MPI_Reduce(&err, &comm_err, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_CART);
      if (world_rank == 0)
        printf("ERROR[%d]::%10.10lf;\n",s,comm_err);
    } else if (s % 3 == 0){
      u_step(A, C, B, &g, &l, s);
      swap_halo(A, &l, &g);
      after_swap_copy_j0_to_j1(A, &l);
      // err = u_full_err(A, s, &g, &l);
      // printf("FE[%d, %d, %d][%d]::%lf\n",l.x,l.y,l.z,s,err);
      err = u_err(A, s, &g, &l);
      MPI_Reduce(&err, &comm_err, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_CART);
      if (world_rank == 0)
        printf("ERROR[%d]::%10.10lf;\n",s,comm_err);
    }
    MPI_Barrier(MPI_COMM_CART);
  }

  free(A);
  free(B);
  free(C);

  clock_t end = clock();
  double time = ((double)(end - start)) / CLOCKS_PER_SEC;

  if (world_rank == 0){
    double time_max;
    MPI_Reduce(&time, &time_max, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    printf("TIME=%lf;\n ", time_max);
  } else {
    MPI_Reduce(&time, NULL, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);  
  }

  MPI_Finalize();


  return 0;
}

