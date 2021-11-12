/*Victor Lelis Soares
2019.1904.038.2*/
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <time.h>

int main(int argc, char **argv){
   long int n, max, num, c;
   time_t t1;

   n = atoi(argv[1]);
   max = 999;
   
   srand((unsigned)time(&t1));
   

   printf("%ld\n", n);
   for (c = 1; c <= n; c++){
      num = rand() % max;
      printf("%ld ",num);        
   }

   
   return 0;
}