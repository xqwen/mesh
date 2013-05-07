main: main.o qtSSNP.o qtMSNP.o controller.o qtMSNP_ES.o qtMSNP_EE.o qtMSNP_CEFES.o qtMSNP_CEFEE.o paramset.o
	g++ -O3 main.o qtSSNP.o qtMSNP.o qtMSNP_ES.o qtMSNP_EE.o qtMSNP_CEFES.o qtMSNP_CEFEE.o controller.o paramset.o -lm -L /usr/local/lib -lgsl -lgslcblas -o mesh
main.o: main.cc 
	g++ -c -O3 main.cc	
paraset.o: paramset.cc paramset.h
	g++ -c -O3 paramset.cc
controller.o: controller.cc controller.h
	g++ -c -O3 controller.cc
qtSSNP.o: qtSSNP.cc qtSSNP.h
	g++ -c -O3 qtSSNP.cc
qtMSNP.o: qtMSNP.cc qtMSNP.h
	g++ -c -O3 qtMSNP.cc
qtMSNP_ES.o: qtMSNP_ES.h qtMSNP_ES.cc
	g++ -c -O3 qtMSNP_ES.cc
qtMSNP_EE.o: qtMSNP_EE.h qtMSNP_EE.cc
	g++ -c -O3 qtMSNP_EE.cc
qtMSNP_CEFES.o: qtMSNP_CEFES.h qtMSNP_CEFES.cc
	g++ -c -O3 qtMSNP_CEFES.cc
qtMSNP_CEFEE.o: qtMSNP_CEFEE.h qtMSNP_CEFEE.cc
	g++ -c -O3 qtMSNP_CEFEE.cc
clean:
	rm *.o mesh
