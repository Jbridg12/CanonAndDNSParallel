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
    MPI_Comm NEW_COMM_WORLD;
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
    p = cbrt(P);
    int dims[3];
    int period[3];
    
    dims[0] = p;
    dims[1] = p;
    dims[2] = p;
    
    period[0] = 1;
    period[1] = 1;
    period[2] = 1;
    
    MPI_Cart_create(MPI_COMM_WORLD, 3, dims, period, 1, &NEW_COMM_WORLD);
    MPI_Comm_size(MPI_COMM_WORLD,&nprocs);  
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	
    if(NEW_COMM_WORLD == MPI_COMM_NULL)
        return EXIT_FAILURE;

    int coords[3];
    MPI_Cart_coords(NEW_COMM_WORLD, rank, 3, coords);

    int A[n][n];
    int B[n][n];
    int **localA, **localB, **lsum, **tsum;
    localA = alloc_2d_int(n/p);
    localB = alloc_2d_int(n/p);
    lsum = alloc_2d_int(n/p);
    tsum = alloc_2d_int(n/p);
    int rankA, rankB;
    int coord[3];
    int coordA[3];
    int coordB[3];
    int recvA[(n/p)*(n/p)];
    int recvB[(n/p)*(n/p)];
    int i, j, k, l, h;
    int ranki;
    int bufA[n*n];
    int bufB[n*n];
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
        
        for(i = 0; i <(n/p)*(n/p); i++)
        {
            recvA[i] = bufA[i];
            recvB[i] = bufB[i];
        }
        for(i = 0; i < p; i++)
        {
            for(j = 0; j < p; j++)
            {
                if(i == 0 && j == 0)
                    continue;
    
                int index = ((p * i) + 1) * (n/p);
                coord[0] = i; coord[1] = j; coord[2] = 0;
                MPI_Cart_rank(NEW_COMM_WORLD, coord, &ranki);
                MPI_Send(&(bufA[index]), (n/p)*(n/p), MPI_INT, ranki, 0, NEW_COMM_WORLD);
                MPI_Send(&(bufB[index]), (n/p)*(n/p), MPI_INT, ranki, 0, NEW_COMM_WORLD);
            }
            
        }
    }
    else if(coords[2] == 0)
    {
        MPI_Recv(recvA, (n/p)*(n/p), MPI_INT, 0, 0, NEW_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(recvB, (n/p)*(n/p), MPI_INT, 0, 0, NEW_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    
    if(rank ==0)
    {
        commTime += MPI_Wtime()*1000;
        
    }
    
    


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
        /*
        for(i = 0; i < n/p; i++)
        {
            for(j = 0; j < n/p; j++)
            {
                printf("rank: %d  L[%d][%d]: %d\n", rank, i, j, lsum[i][j]);
            }
            
        }
        */
        printf("\nTotalTime: %f\n", compTime);
        return 0;
        
    }
    
    for(i = 0; i < n/p; i++)
    {
        for(j = 0; j < n/p; j++)
        {
            localA[i][j] = 0;
            localB[i][j] = 0;
            if(coords[2] == 0)
            {
                localA[i][j] = recvA[(n/p) * i + j];
                localB[i][j] = recvB[(n/p) * i + j];
            }
            lsum[i][j] = 0;
        }
        
    }

    int chord[3] = {1,1,0};
    MPI_Cart_rank(NEW_COMM_WORLD, chord, &ranki);
    
    
    if(rank ==0)
    {
        commTime -= MPI_Wtime()*1000;
        
    }
    


    
    if(coords[2] == 0)
    {
        if(coords[1] != 0)
        {
            coordA[0] = coords[0]; coordA[1] = coords[1]; coordA[2] = coords[1];
            MPI_Cart_rank(NEW_COMM_WORLD, coordA, &rankA);
            MPI_Send(&(localA[0][0]), (n/p)*(n/p), MPI_INT, rankA, 0, NEW_COMM_WORLD);
        }
        
        if(coords[0] != 0)
        {
            coordB[0] = coords[0]; coordB[1] = coords[1]; coordB[2] = coords[0];
            MPI_Cart_rank(NEW_COMM_WORLD, coordB, &rankB);
            MPI_Send(&(localB[0][0]), (n/p)*(n/p), MPI_INT, rankB, 0, NEW_COMM_WORLD);
        }
    }
    else
    {
        if(coords[0] == coords[2])
        {
            coordA[0] = coords[0]; coordA[1] = coords[1]; coordA[2] = 0;
        
            MPI_Cart_rank(NEW_COMM_WORLD, coordA, &rankA);
            
            MPI_Recv(&(localA[0][0]), (n/p)*(n/p), MPI_INT, rankA, 0, NEW_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        if(coords[1] == coords[2])
        {
            coordA[0] = coords[0]; coordA[1] = coords[1]; coordA[2] = 0;
        
            MPI_Cart_rank(NEW_COMM_WORLD, coordA, &rankB);
            MPI_Recv(&(localB[0][0]), (n/p)*(n/p), MPI_INT, rankB, 0, NEW_COMM_WORLD, MPI_STATUS_IGNORE);
        }

    }
    
    if(rank == 0)
    {
        commTime += MPI_Wtime()*1000;
        
    }
    
    MPI_Barrier(NEW_COMM_WORLD);



    coordA[0] = coords[0]; coordA[1] = coords[2]; coordA[2] = coords[2];
    coordB[0] = coords[2]; coordB[1] = coords[1]; coordB[2] = coords[2];
        
    MPI_Cart_rank(NEW_COMM_WORLD, coordB, &rankA);
    MPI_Cart_rank(NEW_COMM_WORLD, coordA, &rankB);
        
    if(rank ==0)
    {
        commTime -= MPI_Wtime()*1000;
        
    }
    MPI_Bcast(&(localA[0][0]), (n/p)*(n/p), MPI_INT, rankA, NEW_COMM_WORLD);
    MPI_Bcast(&(localB[0][0]), (n/p)*(n/p), MPI_INT, rankB, NEW_COMM_WORLD);
    
    if(rank == 0)
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
    
    
    MPI_Barrier(NEW_COMM_WORLD);
    
    coordA[0] = coords[0]; coordA[1] = coords[1]; coordA[2] = 0;
    
    MPI_Cart_rank(NEW_COMM_WORLD, coordA, &rankA);
    int lint, tint;
    
    if(rank ==0)
    {
        compTime -= MPI_Wtime()*1000;
        
    }
    MPI_Reduce(&(lsum[0][0]), &(tsum[0][0]), (n/p)*(n/p), MPI_INT, MPI_SUM, rankA, NEW_COMM_WORLD);

    
    if(rank ==0)
    {
        compTime += MPI_Wtime()*1000;
        
    }
    
    MPI_Barrier(NEW_COMM_WORLD);
    
    if(rank ==0)
    {
        printf("\nNprocs: %d CommTime: %f  CompTime: %f TotalTime: %f\n", nprocs, commTime, compTime, commTime+compTime);
    }
    
    free(localA[0]);
    free(localB[0]);
    free(lsum[0]);
    free(tsum[0]);
    free(tsum);
    free(lsum);
    free(localA);
    free(localB);
    MPI_Finalize();
}