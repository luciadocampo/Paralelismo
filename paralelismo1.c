#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>


int main(int argc, char *argv[])
{
    int i, done = 0, n, count, aux, allProcs, proceso; //   Establecemos variables int para los procesos.
    double PI25DT = 3.141592653589793238462643;
    double pi, x, y, z, aprox;

    MPI_Status status;                      //Establecemos un estado para el MPI.
    MPI_Init(&argc,&argv);      //Inicializamos el MPI.
    MPI_Comm_size(MPI_COMM_WORLD, &allProcs); //Obtenemos el número de procesos.
    MPI_Comm_rank(MPI_COMM_WORLD,&proceso); //Obtenemos el número del proceso.

    while (!done)
    {
        if(proceso == 0){    //El proceso 0 pide los datos al usuario que ejecute el programa.
        printf("Enter the number of points: (0 quits) \n");
        scanf("%d",&n);

        for(aux=1; aux<allProcs; aux++){  //Se envian los datos metidos al resto de los procesos.
            MPI_Send(&n,1,MPI_INT,aux,0,MPI_COMM_WORLD);
        }
        }else{
            //En caso de no ser el proceso 0, los demás procesos reciben el mensaje del proceso 0.
            MPI_Recv(&n,1,MPI_INT,0,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
        }

        if (n == 0) break; //Si el buf n es 0, el programa termina.

        count = 0;

        for (i = proceso + 1; i <= n; i+=allProcs) {
            // Get the random numbers between 0 and 1
            x = ((double) rand()) / ((double) RAND_MAX);
            y = ((double) rand()) / ((double) RAND_MAX);

            // Calculate the square root of the squares
            z = sqrt((x*x)+(y*y));

            // Check whether z is within the circle
            if(z <= 1.0)
                count++;
        }
        pi = ((double) count/(double) n)*4.0;

        //A estas alturas del código, todos los procesos tienen su valor local calculado de pi.

        //El resto de procesos envían las aproximaciones parciales al proceso 0.
        if(proceso>0){
            MPI_Send(&pi,1,MPI_DOUBLE,0,0,MPI_COMM_WORLD);
        }else{ //El proceso 0 se encarga de recolectar las aproximaciones enviadas por el resto de procesos y las imprime.
            for(aux = 1; aux<allProcs; aux++){
                MPI_Recv(&aprox,1,MPI_DOUBLE,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
                pi+=aprox;
            }
        }
        printf("pi is approx. %.16f, Error is %.16f\n", pi, fabs(pi - PI25DT));
    }
    //Esperamos a que los procesos terminen y liberamos los recursos reservados con el MPI_Finalize();
    MPI_Finalize();
}