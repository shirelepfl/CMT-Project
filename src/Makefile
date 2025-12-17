all: project run matlab

project: C_code.c
	gcc -Wall C_code.c -o project -lm

run: project
	./project

matlab:
	/usr/local/bin/matlab-2021b -batch "run('Matlab_code.m')"

clean:
	rm -f project

.PHONY: all run matlab clean


