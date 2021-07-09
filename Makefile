# makefile - leonard
# compile program by typing make
# this thing does wonders, even though i'm still learning
# you should try
# ahh - this particular file only works with unix i guess ¯\_(ツ)_/¯

CC = g++
CFLAGS = -g -Wall
OBJECTS = main.o Fraction.o Matrix.o

main: $(OBJECTS)
	$(CC) $(CFLAGS) -o test.run $^

main.o: main.cpp Matrix.h
	$(CC) $(CFLAGS) -c main.cpp

Matrix.o:  Matrix.h

Fraction.o: Fraction.h

run: main
	./test.run

.PHONY: clean zip test

clean:
	rm -f *.o *.zip *.run

zip: clean main
	zip matrix_manipulation.zip Matrix.h Matrix.cpp main.cpp Fraction.h Fraction.cpp Makefile 
