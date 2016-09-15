kernel: kernel.c
	gcc -Wall kernel.c -o md -fopenmp -lgomp -lm
clean:
	rm md *.exe* out*
