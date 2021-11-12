/*Victor Lelis Soares
2019.1904.038.2*/
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <mpi.h>

/*
int bitExtracted(int number, int k, int p){
    return (((1 << k) - 1) & (number >> (p - 1)));
}
 */
/*int IncOrder(const void *e1, const void *e2);*/
/*void swap(long int*, long int*);
int partition(int*, int, int, int);
void quickSort(int*, int, int);
int validateSort(int *a, int size);
*/
int IncOrder(const void *e1, const void *e2);
long int *quickSort_p(MPI_Comm comm, long int *local_array, long int local_size, int *arrayLength);
void divide (long int *pivot, int *l, int *u, long int *LList, long int *Ulist, long int *array, int size_array);

int main (int argc, char **argv){
	int n;
	int npes;
	int myrank;
	int nlocal;
  /*guarda todos os elementos lidos pelo processo 0*/
	long int *entryelmnts; 
	long int *elemnst;
	long int *sorted;

	double stime, etime;

	MPI_Status status;

	MPI_Init                (&argc, &argv);
	MPI_Comm_size  (MPI_COMM_WORLD, &npes);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	FILE *fp = fopen(argv[1], "r");/*abre o arquivo em modo leitura*/
	fscanf(fp, "%d", &n);/*número de elementos*/

	sorted = (long int*) malloc(n * sizeof(long int));
	nlocal = (n/npes);/*numero de elementos / numero de processos -> elementos locais*/

	/*separa os elementos do arquivo e distribui*/
	if(myrank == 0){
		int i = 0;
		/*aloca espaço para todos os elementos*/
		entryelmnts = (long int *) malloc(n*sizeof(long int));
		
		/*lê todos os elementos e guarda em entryelmnts*/
		while(fscanf(fp, "%ld", (entryelmnts+i)) != EOF){
			/*printf("%d\n", *(entryelmnts+i));*/
			i++;
		}

        double ctime, ftime;
        ctime = MPI_Wtime();
        qsort(entryelmnts, n, sizeof(long int), IncOrder);
        ftime = MPI_Wtime();
        printf("tempo QuickSequencial = %f\n", ftime - ctime);


		printf("nlocal:%d\n", nlocal);
		elemnst = (long int *) malloc(nlocal * sizeof(long int));
		int ini = nlocal;

		int j = 0;
		while(j != nlocal){
			*(elemnst+j) = *(entryelmnts+j);
			j++;
		}

	
		int prox_process = 1;
		
		while(ini < n){
			
			MPI_Send((entryelmnts+ini), nlocal, MPI_LONG, prox_process, 0, MPI_COMM_WORLD);
			prox_process++;
			ini +=nlocal;
		}
		free(entryelmnts);/*já nao preciso mais da entrada separada*/
	}
	else{
		elemnst = (long int *) malloc(nlocal * sizeof(long int));
		MPI_Recv(elemnst, nlocal, MPI_LONG, 0, 0, MPI_COMM_WORLD, &status);
	}
	MPI_Barrier(MPI_COMM_WORLD);
    stime = MPI_Wtime();
	
	/*int j = 0;
	while(j != nlocal){
		printf("%d->%d\n", myrank, *(elemnst+j));
		j++;
	}*/

	sorted = quickSort_p(MPI_COMM_WORLD, elemnst, nlocal, &nlocal);
	/*quickSort(sorted, 0, nlocal-1);*/
    qsort(sorted, nlocal, sizeof(long int), IncOrder);
	etime = MPI_Wtime();
	MPI_Barrier(MPI_COMM_WORLD);
	/*printf("Ordenado!\n");*/
	
    /*if(myrank == 0){
        int k = 0;
    	while(k != nlocal){
    		printf("%d->%ld\n", myrank, *(sorted+k));
    		k++;
    	}
    }
    MPI_Barrier(MPI_COMM_WORLD);
    if(myrank == 1){
        int k = 0;
        while(k != nlocal){
            printf("%d->%ld\n", myrank, *(sorted+k));
            k++;
        }
    }*/
  
  printf("tempo = %f, %d\n", etime - stime, myrank);

  free(sorted);
  MPI_Finalize();

	return 0;
}

/*Perform serial quicksort
void quickSort(int *arr, int left, int right) {
 if (left < right){
    int pivotIndex = (right+left)/2;
    int pivotNewIndex = partition(arr, left, right, pivotIndex);

    quickSort(arr, left, pivotNewIndex - 1);
    quickSort(arr, pivotNewIndex + 1, right);
  }
}*/

/*Quicksort paralelo*/
long int *quickSort_p(MPI_Comm comm, long int *local_array, long int local_size, int *arrayLength){
  int local_rank, other_size, local_nproc;
  int another_paried;
  int l = 0;
  int u = 0;
  long int pivot, *received_array, *new_local, *LList, *Ulist;
  MPI_Request request;
  MPI_Status status;

  /*extract local rank and local nproc*/
  MPI_Comm_rank(comm, &local_rank);
  MPI_Comm_size(comm, &local_nproc);
  
  LList = (long int *)malloc(local_size * sizeof(long int));
  Ulist = (long int *)malloc(local_size * sizeof(long int));

  if(local_nproc > 1){
  
    /*processo zero divulga seu pivô*/
    if(local_rank == 0){
      /*pivo escolhido*/
      pivot = local_array[local_size/2];
    }
    
    /*Broadcast pivo*/
    MPI_Bcast(&pivot, 1, MPI_DOUBLE, 0, comm);
    
    /*divide a lista em elementos maiores que o pivo(u) e menores(l)*/
    divide(&pivot, &l, &u, LList, Ulist, local_array, local_size);
   /* printf("Numero elementos menores que pivo->%d, %d\n",l, local_rank);
    printf("Numero elementos maiores que pivo->%d, %d\n",u, local_rank);
    */
    if(local_rank < local_nproc / 2){
    	/*processador da parte superior que está pareado com o atual*/
    	another_paried = local_rank + (local_nproc / 2);
    	/*envia para os processadores mais altos, seus elementos maiores que o pivô*/
    	MPI_Isend(Ulist, u, MPI_LONG, another_paried, 111, comm, &request);

    	/*verifica o que foi recebido, para alocar o buffer*/
    	MPI_Probe(another_paried, 222, comm, &status);
    	MPI_Get_count(&status, MPI_LONG, &other_size);

    	received_array = (long int *)malloc(other_size*sizeof(long int));

    	/*recebe os menores elementos*/
      MPI_Recv(received_array, other_size, MPI_LONG, another_paried, 222, comm, &status);
      MPI_Wait(&request, &status);

    	/*calculando tamanho para o array mesclado*/
    	local_size = l + other_size;
    	new_local = (long int *) malloc(local_size * sizeof(long int));

      /*	
      for (i = 0; i < other_size; ++i) {
        printf("v->%d, p->%d\n", received_array[i], local_rank);
      }*/
      /*copia os menores elementos para o novo vetor*/
      memcpy(new_local, LList, l*sizeof(long int));
      /*concatena os elementos recebidos com o novo vetor*/
      memcpy(new_local+l, received_array, other_size*sizeof(long int));
      /*
      for (i = 0; i < local_size; ++i) {
      	printf("v->%d, p->%d\n", new_local[i], local_rank);	  
      }*/
    }

    /*processador se encontra na parte superior*/
    if( local_rank >= local_nproc / 2){
    	another_paried = local_rank - (local_nproc / 2);

    	/*envia para os processadores mais baixos, seus elementos menores que o pivô*/
    	MPI_Isend(LList, l, MPI_LONG, another_paried, 222, comm, &request);

    	/*verifica o que foi recebido, para alocar o buffer*/
    	MPI_Probe(another_paried, 111, comm, &status);
    	MPI_Get_count(&status, MPI_LONG, &other_size);

    	received_array = (long int *)malloc(other_size*sizeof(long int));

    	/*recebe os maiores elementos*/
      MPI_Recv(received_array, other_size, MPI_LONG, another_paried, 111, comm, &status);
      MPI_Wait(&request, &status);

    	/*printf("recebi->%d, %d\n", other_size, local_rank);*/

    	/*calculando tamanho para o array mesclado*/
    	local_size = u + other_size;
    	new_local = (long int *) malloc(local_size * sizeof(long int));

    	/*copia os menores elementos para o novo vetor*/
      memcpy(new_local, Ulist, u*sizeof(long int));
      /*concatena os elementos recebidos com o novo vetor*/
      memcpy(new_local+u, received_array, other_size*sizeof(long int));
    	/*
      for (i = 0; i < local_size; ++i) {
      	printf("v->%d, p->%d\n", new_local[i], local_rank);	  
      }
      */
    }

    free(local_array);
    free(received_array);
    free(LList);
    free(Ulist);

    MPI_Comm new_comm;

    int split  = local_rank / (local_nproc / 2);

    MPI_Comm_split(comm, split, 0, &new_comm);
    *arrayLength = local_size;
    return quickSort_p(new_comm, new_local, local_size, arrayLength);
  }
  else{
		*arrayLength = local_size;
    return local_array;  	
  }
}

/*particiona de acordo com o indice do pivo
int partition (int *arr, int left, int right, int pivotIndex){
  int pivotValue = arr[pivotIndex], temp = arr[right];
  int storeIndex = left, i;

  arr[right] = arr[pivotIndex];
  arr[pivotIndex] = temp;

  for(i = left; i < right; i++) {
    if (arr[i] < pivotValue) {
      temp = arr[storeIndex];
      arr[storeIndex] = arr[i];
      arr[i] = temp;
      storeIndex++;
    }
  }

  temp = arr[right];
  arr[right] = arr[storeIndex];
  arr[storeIndex] = temp;

  return storeIndex;
}*/

/* troca dois valores*/
/*void swap (long int *x, long int *y) {
  long int temp = *x;
  *x = *y;
  *y = temp;
}*/

/* Check if array is sorted
int validateSort (int *a, int size) {
  int i;
  for (i = 0; i < size-1; i++) {
    if (a[i] > a[i+1]) {
      printf("%d - %d\n", a[i],a[i+1]);
      return 0;
    }
  }

  return 1;
}*/

/*divide a lista em duas, uma maior que o pivo(u) e outra menor(l)*/
void divide (long int *pivot, int *l, int *u, long int *LList, long int *Ulist, long int *array, int size_array){
    int i = 0;

    for(; i < size_array; i++){

        if( IncOrder(array+i, pivot) <= 0){
          /**(array+i) <= *pivot*/
            *( LList + *(l) ) =  *(array+i);
            ++(*l);
        }
        else{
            *( Ulist + *(u) ) =  *(array+i);
            ++(*u);
        }

    }

}

/*
retorno da função:
  valor <= 0: e1 vem antes de e2 ou os dois são iguais.
  valor  > 0: e1 vem depois de e2
*/
/*int IncOrder(const void *e1, const void *e2){
  
  //printf("recebi-> %d, %d\n", *((int *)e1), *((int *)e2));
  return ( *((int *)e1) - *((int *)e2));
}
*/

/*
Máscara de bits:
  22 primeiros bits -> 11111111 11111111 11111100 00000000 00000000 00000000 00000000 00000000 (64 bits), fffffc0000000000 desloca >> 42
  20 bits do meio   -> 00000000 00000000 00000011 11111111 11111111 11000000 00000000 00000000 (64 bits), 3ffffc00000 desloca >> 22
  22 bits do finais -> 00000000 00000000 00000000 00000000 00000000 00111111 11111111 11111111 (64 bits), 3fffff      desloca >> 0
*/

/*função que compara os bits*/
int IncOrder(const void *e1, const void *e2){
  long int mask_ini = 0xfffffc0000000000, mask_meio = 0x3ffffc00000, mask_fim = 0x3fffff;
  
  /*extrai os primeiros : primeiros 22, 20 do meio e  22 ultimos bits*/
  
  long int ini_e1 = *((long int *)e1) & mask_ini  , 
    meio_e1 = *((long int *)e1) & mask_meio, fim_e1 = *((long int *)e1) & mask_fim;
  
  long int ini_e2 = *((long int *)e2) & mask_ini, 
    meio_e2 = *((long int *)e2) & mask_meio, fim_e2 = *((long int *)e2) & mask_fim;

  if(ini_e1 == ini_e2){/*caso os 22 primeiros bits sejam iguais, comparo os próximos*/
    
    if(meio_e1 == meio_e2){/*compara os 20 bits do meio*/
        if(fim_e1 == fim_e2){/*compara os bits finais*/
            /*caso sejam iguais, nesse caso eu uso eles como um todo,
            mas não faria diferença retornar qualquer um dos outros bits
            como todos são iguais*/
            return ( *((int *)e1) - *((int *)e2));
        }
        else{
            return (int) fim_e1 - (int)fim_e2;
        }
    }
    else{/*caso os bits do meio sejam diferentes*/
        return (int) meio_e1 - (int)meio_e2;
    }
  }
  else{/*caso os 22 primeiros sejam diferentes*/
    return (int) ini_e1 - (int)ini_e2;
  }
}