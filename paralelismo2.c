#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

int MPI_FlattreeColectiva(void *sendbuffer, void *recvbuffer, int count, MPI_Datatype datatype, MPI_Op op, int root, MPI_Comm comm){

    //Realizamos el control de errores para los parámetros introducidos en la función
    if(sendbuffer == NULL) {
        return MPI_ERR_BUFFER;
    }else if(count == 0){
        return MPI_ERR_COUNT;
    }else if(datatype != MPI_DOUBLE){
        return MPI_ERR_TYPE;
    }else if(op != MPI_SUM){
        return MPI_ERR_OP;
    }else if(comm != MPI_COMM_WORLD){
        return MPI_ERR_COMM;
    }else{

        int allProcs, proceso, error;
        double *salida = (double *) recvbuffer; //Hacemos un casteo para obtener el valor del buffer que recibimos en la función.
        double *entrada = (double *) sendbuffer; //Hacemos un casteo para obtener el valor del buffer que enviamos con la función.
        double *aux;

        MPI_Status status;                      //Establecemos un estado para el MPI.
        MPI_Comm_size(MPI_COMM_WORLD, &allProcs); //Obtenemos el número de procesos.
        MPI_Comm_rank(MPI_COMM_WORLD,&proceso); //Obtenemos el número del proceso.

        if (proceso == root) { //Se hace esto en vez de poner proceso 0 porque aquí el proceso 0 puede que no sera el proceso raíz.

            aux = malloc(sizeof(double) * count);

            for (int i = 0; i < count; ++i) { //Se cogen los valores del buffer de entrada y se introducen en la salida.
                salida[i] = entrada[i];
            }

            for (int i = 0; i < allProcs; ++i) { //Se recorren todos los procesos excepto el que sea el raíz.

                if (i != root) {
                    error = MPI_Recv(aux, count, datatype, i, 0, comm, &status);
                    //El proceso raíz recibe el mensaje de los otros procesos.

                    if (error != MPI_SUCCESS) {
                        free(aux);
                        return error;
                    }

                    for (int j = 0; j < count; ++j) { //En caso de que haya más de un elemento en el buffer se hace esto.
                        salida[j] += aux[j];
                    }
                }

            }
            free(aux); //Liberamos memoria.
        } else {        //Los procesos que no son el raíz envían el mensaje al proceso raíz.
            error = MPI_Send(sendbuffer, count, datatype, root, 0, comm);
            if (error != MPI_SUCCESS) {
                return error;
            }
        }
    }
    return MPI_SUCCESS;
}

int MPI_BinomialColectiva(void *buffer, int count, MPI_Datatype datatype, int root, MPI_Comm comm){

    //Realizamos el control de errores para los parámetros introducidos en la función
    if(buffer == NULL){
        return MPI_ERR_BUFFER;
    }else if(count == 0){
        return MPI_ERR_COUNT;
    }else if(root != 0) {
        return MPI_ERR_ROOT;
    }else if(comm != MPI_COMM_WORLD){
        return MPI_ERR_COMM;
    }else {

        int allProcs, proceso, error, origen, destino;

        MPI_Status status;                      //Establecemos un estado para el MPI.
        MPI_Comm_size(MPI_COMM_WORLD, &allProcs); //Obtenemos el número de procesos.
        MPI_Comm_rank(MPI_COMM_WORLD,&proceso); //Obtenemos el número del proceso.


        //Ceil lo usamos para redondear al siguiente número entero por encima del que se quede con la operación del logaritmo que puede no ser exacta.
        for (int i = 1; i <= ceil(log2(allProcs)); ++i) { //El tope de las iteraciones al que se llega, es el número más alto que se encuentra para realizar el logartimo de los procesos.
            //Operación de la implementación descrita en el enunciado.
            if (proceso < pow(2, i - 1)) {
                destino = proceso + pow(2, i - 1); //Calculamos el destino del mensaje.
                if (destino < allProcs) {           //Comprobamos que el proceso que va a recibir el mensaje está dentro del rango total de procesos.
                    error = MPI_Send(buffer, count, datatype, destino, 0, comm);
                    if (error != MPI_SUCCESS) {
                        return error;
                    }
                }
            } else {

                if (proceso < pow(2, i)) {          //Comprobamos que el proceso que llega es un posible receptor de la info.
                    origen = proceso - pow(2, i - 1); //Calculamos el proceso origen y recibimos el mensaje de ese proceso.
                    error = MPI_Recv(buffer, count, datatype, origen, 0, comm, &status);
                    if (error != MPI_SUCCESS) {
                        return error;
                    }
                }
            }
        }
    }
    return MPI_SUCCESS;
}


int main(int argc, char *argv[])

//AQUÍ SE COMENTA Y DESCOMENTA SEGÚN LA OPCIÓN COLECTIVA QUE QUIERAS USAR.
{
    int i, done = 0, n, count;
    double PI25DT = 3.141592653589793238462643;
    double pi, x, y, z;

    double aprox;
    int allProcs, proceso;

    MPI_Init(&argc,&argv);      //Inicializamos el MPI.
    MPI_Comm_size(MPI_COMM_WORLD, &allProcs); //Obtenemos el número de procesos.
    MPI_Comm_rank(MPI_COMM_WORLD,&proceso); //Obtenemos el número del proceso.


    while (!done)
    {
        if(proceso == 0){ //El proceso 0 pide los datos al usuario que ejecute el programa.
            printf("Enter the number of points: (0 quits) \n");
            scanf("%d",&n);

        }
        //El proceso 0 envía el valor de n a los demás procesos y estos lo reciben con el Broadcast.
        MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_BinomialColectiva(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);

        if (n == 0) break; //Si el buf n es 0, el programa termina.

        count = 0;

        for (i = proceso + 1; i <= n; i += allProcs) {
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

        //Se realiza la suma de todas las aproximaciones parciales de pi calculadas por los procesos que no son el raíz
         MPI_Reduce(&pi, &aprox, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
         //MPI_FlattreeColectiva(&pi,&aprox,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);

        if(proceso == 0){  //Imprime el proceso 0 el mensaje del cálculo de pi.
            printf("pi is approx. %.16f, Error is %.16f\n", aprox, fabs(aprox - PI25DT));
        }
    }

    //Esperamos a que los procesos terminen y liberamos los recursos reservados con el MPI_Finalize();
    MPI_Finalize();
}