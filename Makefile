prep_graph : main.o graph.o
	g++ -O3 -o prep_graph main.o graph.o -lhts
main.o : main.cpp graph.h cmdline.h
	gcc -O3 -c main.cpp
graph.o : graph.cpp graph.h
	gcc -O3 -c graph.cpp
clean :
	rm -f graph.o main.o prep_graph