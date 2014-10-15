main: main.o eQTL.o controller.o
	g++   -fopenmp main.o eQTL.o controller.o -lm -lgsl -lgslcblas -static -o hm_em 
main.o: main.cc 
	g++ -c  main.cc	
eQTL.o: eQTL.cc classdef.h
	g++ -c  eQTL.cc
controller.o: controller.cc classdef.h
	g++  -fopenmp  -c -O3 controller.cc
clean:
	rm *.o hm_em	
