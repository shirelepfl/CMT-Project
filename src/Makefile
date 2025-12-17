all: project run matlab

project: C_code.c
	gcc -Wall C_code.c -o project -lm

run: project
	./project

matlab:
	matlab -batch "run('Matlab_code.m'); exit"

clean:
	rm -f project

.PHONY: all run matlab clean


