spkmeans: spkmeans.o kmeans.o kmeans.h spkmeans.h
	gcc -o spkmeans spkmeans.o kmeans.o -lm
kmeans.o: kmeans.c
	gcc -c -g -ansi -Wall -Wextra -Werror -pedantic-errors kmeans.c
spkmeans.o: spkmeans.c
	gcc -c -g -ansi -Wall -Wextra -Werror -pedantic-errors spkmeans.c

