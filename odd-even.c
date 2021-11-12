/*Victor Lelis Soares
2019.1904.038.2*/
#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#include <time.h>

/*
Máscara de bits:
  22 primeiros bits -> 11111111 11111111 11111100 00000000 00000000 00000000 00000000 00000000 (64 bits), fffffc0000000000 desloca >> 42
  20 bits do meio   -> 00000000 00000000 00000011 11111111 11111111 11000000 00000000 00000000 (64 bits), fffff00000000    desloca >> 22
  22 bits do finais -> 00000000 00000000 00000000 00000000 00000000 00111111 11111111 11111111 (64 bits), 3ffffc00000      desloca >> 0
*/

int IncOrder(const void *e1, const void *e2);
void compareSplit(int nlocal, long int *elmnts, long int *relmnts, long int *wspace, int keepsmall);

int main (int argc, char **argv){
	int n;
	int npes;
	int myrank;
	int nlocal;
	long int *elmnts;/*guarda os elementos locais*/
	long int *relmnts; /*guarda os elementos recebidos*/
	int oddrank;/*rank do processador durante a fase odd*/
	int evenrank;/*rank do processador durante a fase even*/
	long int *wspace; /*para a fase compare e separa*/
	int i;

	MPI_Status status;
	double stime, etime;
	MPI_Init       (&argc, &argv);
	MPI_Comm_size  (MPI_COMM_WORLD, &npes);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);


	FILE *fp = fopen(argv[1], "r");/*abre o arquivo em modo leitura*/
	fscanf(fp, "%d", &n);/*número de elementos*/

	nlocal = (n/npes);/*numero de elementos / numero de processos -> elementos locais*/

	/*alocando espaço*/
	elmnts  =  (long int *) malloc(nlocal * sizeof(long int));
	relmnts =  (long int *) malloc(nlocal * sizeof(long int));
	wspace  =  (long int *) malloc(nlocal * sizeof(long int));

	/*separa arquivos e distribui*/
	if(myrank == 0){
		long int *entryelmnts = NULL;
		int i = 0;
		/*aloca espaço para todos os elementos*/
		entryelmnts = (long int *) malloc(n*sizeof(long int));
		
		/*lê todos os elementos e guarda em entryelmnts*/
		while(fscanf(fp, "%ld", (entryelmnts+i)) != EOF){
			/*printf("%lf\n", *(entryelmnts+i));*/
			i++;
		}

				double ctime, ftime;
        ctime = MPI_Wtime();
        qsort(entryelmnts, n, sizeof(long int), IncOrder);
        ftime = MPI_Wtime();
        printf("tempo QuickSequencial = %f\n", ftime - ctime);

		/*printf("nlocal:%d\n", nlocal);*/
		elmnts = (long int *) malloc(nlocal * sizeof(long int));
		int ini = nlocal;

		int j = 0;
		while(j != nlocal){
			*(elmnts+j) = *(entryelmnts+j);
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
		elmnts = (long int *) malloc(nlocal * sizeof(long int));
		MPI_Recv(elmnts, nlocal, MPI_LONG, 0, 0, MPI_COMM_WORLD, &status);
		
	}

	MPI_Barrier(MPI_COMM_WORLD);
	stime = MPI_Wtime();


	qsort(elmnts, nlocal, sizeof(long int), IncOrder);

	if(myrank % 2 == 0){
		oddrank = myrank - 1;
		evenrank = myrank + 1;
	}
	else{
		oddrank = myrank + 1;
		evenrank = myrank - 1;	
	}

	if(oddrank == -1 || oddrank == npes)
		oddrank = MPI_PROC_NULL;


	if(evenrank == -1 || evenrank == npes)
		evenrank = MPI_PROC_NULL;


	for(i=0; i < npes-1; ++i){
		if(i % 2 == 1)/*fase impar*/
		{
			MPI_Sendrecv(elmnts, nlocal, MPI_LONG, oddrank, 1, relmnts,
				nlocal, MPI_LONG, oddrank, 1, MPI_COMM_WORLD, &status);
		}
		else/*Fase par*/
		{
			MPI_Sendrecv(elmnts, nlocal, MPI_LONG, evenrank, 1, relmnts,
				nlocal, MPI_LONG, evenrank, 1, MPI_COMM_WORLD, &status);	
		}

		/*compare e separa*/
		compareSplit(nlocal, elmnts, relmnts, wspace, myrank < status.MPI_SOURCE);
	}

	etime = MPI_Wtime();
	MPI_Barrier(MPI_COMM_WORLD);

	
/*if(myrank == 0){
		for(i=0; i<nlocal; i++)
			printf("%d->%ld ", myrank, elmnts[i]);
		printf("\n");
	}
	MPI_Barrier(MPI_COMM_WORLD);
	if(myrank == 1){
		for(i=0; i<nlocal; i++)
			printf("%d->%ld ", myrank, elmnts[i]);
		printf("\n");
	}*/
	

	printf("tempo = %f\n", etime - stime);

	free(elmnts);
	free(relmnts);
  free(wspace);
	MPI_Finalize();

	return 0;
}

void compareSplit(int nlocal, long int *elmnts, long int *relmnts, long int *wspace, int keepsmall){
	int i, j, k;

	for(i = 0; i < nlocal; i++){
		wspace[i] = elmnts[i];/*copia elementos*/
	}

	if(keepsmall)/*mantem os menores elementos*/
	{
		/*merge*/
		for(i = j = k = 0; k < nlocal; k++){
			/*if(j == nlocal || (i < nlocal && wspace[i] < relmnts[j]))*/
			if(j == nlocal || (i < nlocal && IncOrder(wspace+i, relmnts+j) < 0))
				elmnts[k] = wspace[i++];
			else
				elmnts[k] = relmnts[j++];
		}

	}
	else{/*mantem os maiores elementos*/
		for(j = i = k = nlocal-1; k>=0; k--){
			/*if(j == 0 || (i >= 0 && wspace[i] >= relmnts[j]))*/
			if(j == 0 || (i >= 0 && IncOrder(wspace+i, relmnts+j) >= 0))
				elmnts[k] = wspace[i--];
			else
				elmnts[k] = relmnts[j--];
		}
	}
}

/*função original para comparação
int IncOrder(const void *e1, const void *e2){
	return ( *((int *)e1) - *((int *)e2) );
}*/

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