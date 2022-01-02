#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>

int **alloc_2d_int(int n) 
{
    int i;
    int *data = (int *)malloc(n*n*sizeof(int));
    int **array= (int **)malloc(n*sizeof(int*));
    for (i=0; i<n; i++)
        array[i] = &(data[n*i]);

    return array;
}

int main(int argc, char** argv)
{
    int n, P, p, rank, nprocs, root, x, y, shift, dest;
    double commTime = 0.0;
    double compTime = 0.0;
    MPI_Init(&argc, &argv);  
    MPI_Comm_size(MPI_COMM_WORLD,&nprocs);  
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  
    if(argc != 3)
    {
        printf("Incorrect amount of arguments");
        return;
    }
      
    n = atoi(argv[1]);
    P = nprocs;
    p = sqrt(P);
      
    int A[n][n];
    int B[n][n];
      
    int recvA[(n/p)*(n/p)];
    int recvB[(n/p)*(n/p)];
    int i, j, k, l, h;
    int bufA[n*n];
    int bufB[n*n];
    int lsum[n/p][n/p];
    int z = 0;
    
    if(rank == 0)
    {
        for(i = 0; i < n; i++)
        {
            for(j = 0; j < n; j++)
            {
                A[i][j] = 1;
                B[i][j] = 1;
            }
        }
        

        int f = 0;
        for(i = 0; i < p; i++)
        {
          for(j = 0; j < p; j++)
          {
            for(k = 0; k < n/p; k++)
            {
              for(l = 0; l < n/p; l++)
              {
                bufA[f] = A[i*(n/p)+k][j*(n/p)+l];
                bufB[f] = B[i*(n/p)+k][j*(n/p)+l];
                f++;
              }
            }
          }
        }
        
    }
    
    if(rank ==0)
    {
        commTime -= MPI_Wtime()*1000;    
    }
    
    MPI_Scatter(bufA, (n*n)/P, MPI_INT, recvA, (n*n)/P, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Scatter(bufB, (n*n)/P, MPI_INT, recvB, (n*n)/P, MPI_INT, 0, MPI_COMM_WORLD);
    if(rank ==0)
    {
        commTime += MPI_Wtime()*1000;    
    }
    int **localA, **localB, **newA, **newB;
    localA = alloc_2d_int(n/p);
    localB = alloc_2d_int(n/p);
    newA = alloc_2d_int(n/p);
    newB = alloc_2d_int(n/p);

    if(nprocs == 1)
    {
        compTime -= MPI_Wtime()*1000; 
        
        for(i = 0; i < n/p; i++)
        {
            for(j = 0; j < n/p; j++)
            {
                for(k = 0; k < n/p; k++)
                {
                    lsum[i][j] += A[i][k] * B[k][j];
                }
               
            }
            
        }

        compTime += MPI_Wtime()*1000;    
        printf("\nTotalTime: %f\n", compTime);
        return 0;
        
    }
    for(i = 0; i < n/p; i++)
    {
        for(j = 0; j < n/p; j++)
        {
            
            localA[i][j] = recvA[(n/p) * i + j];
            localB[i][j] = recvB[(n/p) * i + j];
            lsum[i][j] = 0;
        }
        
    }
    
    if(rank ==0)
    {
        commTime -= MPI_Wtime()*1000;    
    }
    if(rank >= p)
    {
        shift = rank/p;
        root = shift * p;
        
        if((rank-shift)/p != shift)
        {
            x = (rank-shift) % p;
            //printf("ranko {%d} : %d\n", rank, shift);
            if(rank + shift >= p * (shift + 1))
            {
                MPI_Sendrecv(&(localA[0][0]), (n*n)/P, MPI_INT, (root+x), 0, &(newA[0][0]), (n*n)/P, MPI_INT, ((rank + shift) % p) + root, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
            else
            {
                MPI_Sendrecv(&(localA[0][0]), (n*n)/P, MPI_INT, (root+x), 0, &(newA[0][0]), (n*n)/P, MPI_INT, rank + shift, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
        }
        else
        {
            if(rank + shift >= p * (shift + 1))
            {
                MPI_Sendrecv(&(localA[0][0]), (n*n)/P, MPI_INT, (rank-shift), 0, &(newA[0][0]), (n*n)/P, MPI_INT, ((rank + shift) % p) + root, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
            else
            {
                MPI_Sendrecv(&(localA[0][0]), (n*n)/P, MPI_INT, (rank-shift), 0, &(newA[0][0]), (n*n)/P, MPI_INT, (rank+shift), 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
            
        }
        localA = newA;
    }   
    
    
    if(rank % p != 0)
    {
        shift = p *(rank % p);
        if((rank-shift) < 1)
        {
            y = ((shift-rank) / p);
            root = (P - p + (shift/p));
            if(rank + shift > P)
            {
                int loop = (rank + shift - P)/p;
                
                MPI_Sendrecv(&(localB[0][0]), (n*n)/P, MPI_INT, (root-(y*p)), 0, &(newB[0][0]), (n*n)/P, MPI_INT, (loop * p) + (rank % p), 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
            else
            {
                MPI_Sendrecv(&(localB[0][0]), (n*n)/P, MPI_INT, (root-(y*p)), 0, &(newB[0][0]), (n*n)/P, MPI_INT, (rank + shift), 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
        }
        else
        {
            if(rank + shift > P)
            {
                int loop = (rank + shift - P)/p;
                
                MPI_Sendrecv(&(localB[0][0]), (n*n)/P, MPI_INT, (rank-shift), 0, &(newB[0][0]), (n*n)/P, MPI_INT, (loop * p) + (rank % p), 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
            else
            {
                MPI_Sendrecv(&(localB[0][0]), (n*n)/P, MPI_INT, (rank-shift), 0, &(newB[0][0]), (n*n)/P, MPI_INT, (rank + shift), 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
            
        }    
        localB = newB;
    }
    
    if(rank ==0)
    {
        commTime += MPI_Wtime()*1000;    
    }
    
    if(rank ==0)
    {
        compTime -= MPI_Wtime()*1000;    
    }
    
    for(i = 0; i < n/p; i++)
    {
        for(j = 0; j < n/p; j++)
        {
            for(k = 0; k < n/p; k++)
            {
                lsum[i][j] += localA[i][k] * localB[k][j];
            }    
        }
        
    }
    
    if(rank ==0)
    {
        compTime += MPI_Wtime()*1000;    
    }
    
    MPI_Barrier(MPI_COMM_WORLD);
    
    for(h = 1; h < p; h++)
    {
        
        if(rank ==0)
        {
            commTime -= MPI_Wtime()*1000;    
        }
        if(rank < p)
        {
            dest = P - p + (rank % p);
            MPI_Sendrecv(&(localB[0][0]), (n*n)/P, MPI_INT, dest, 0, &(newB[0][0]), (n*n)/P, MPI_INT, rank + p, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        else if(rank >= (P-p))
        {
            MPI_Sendrecv(&(localB[0][0]), (n*n)/P, MPI_INT, rank-p, 0, &(newB[0][0]), (n*n)/P, MPI_INT, rank % p, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        else
        {
            MPI_Sendrecv(&(localB[0][0]), (n*n)/P, MPI_INT, rank-p, 0, &(newB[0][0]), (n*n)/P, MPI_INT, rank + p, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        
        if(rank % p == 0)
        {
            dest = rank-1+p;
            MPI_Sendrecv(&(localA[0][0]), (n*n)/P, MPI_INT, dest, 0, &(newA[0][0]), (n*n)/P, MPI_INT, rank + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        else if(rank % p == p-1)
        {
            MPI_Sendrecv(&(localA[0][0]), (n*n)/P, MPI_INT, rank-1, 0, &(newA[0][0]), (n*n)/P, MPI_INT, rank+1-p, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        else
        {
            MPI_Sendrecv(&(localA[0][0]), (n*n)/P, MPI_INT, rank-1, 0, &(newA[0][0]), (n*n)/P, MPI_INT, rank+1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        
        if(rank ==0)
        {
            commTime += MPI_Wtime()*1000;    
        }
        
        MPI_Barrier(MPI_COMM_WORLD);
        localA = newA;
        localB = newB;

        
        if(rank ==0)
        {
            compTime -= MPI_Wtime()*1000;    
        }
        for(i = 0; i < n/p; i++)
        {
            for(j = 0; j < n/p; j++)
            {
                for(k = 0; k < n/p; k++)
                {
                    //printf("rank{%d} lsum[i][j]=%d  localA[i][j]=%d * localB[i][j]=%d\n", rank, lsum[i][j], localA[i][j], localB[i][j]);
                    lsum[i][j] += localA[i][k] * localB[k][j];
                }
            }
            
        }
        
        if(rank ==0)
        {
            compTime += MPI_Wtime()*1000;    
        }
        
    
        MPI_Barrier(MPI_COMM_WORLD);
    }
    
    if(rank ==0)
    {
        printf("\nNprocs: %d CommTime: %f  CompTime: %f TotalTime: %f\n", nprocs, commTime, compTime, commTime+compTime);
    }

    
    free(localA[0]);
    free(localB[0]);
    free(localA);
    free(localB);
  	MPI_Finalize();
}


