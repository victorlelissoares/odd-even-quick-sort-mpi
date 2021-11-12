# Compile all
all: quick oddeven generate

#Compile quick
quick: quicksort-mpi.c
	mpicc quicksort-mpi.c -g -o quick -Wall -ansi

#Compile oddeven
oddeven: odd-even.c
	mpicc odd-even.c -g -o oddeven -Wall -ansi

#Compile oddeven
generate: generate_numbers.c
	gcc generate_numbers.c -g -o generate


# Clean up!
clean:
	rm -f *~
	rm -f \#*
	rm -f *.bak *.o core
	rm -f oddeven
	rm -f quick
	rm -f generate
