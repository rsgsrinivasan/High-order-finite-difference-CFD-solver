
SIM: main.o 6th_jacobian.o flow_initialisation.o solver.o roe_average.o viscous_terms.o characteristic_plane.o tecplot_output.o
	mpicc -O3 -lm -std=c99 main.o 6th_jacobian.o flow_initialisation.o solver.o roe_average.o viscous_terms.o characteristic_plane.o tecplot_output.o -o SIM
main.o: main.c function.h
	mpicc -O3    -c -std=c99 main.c
6th_jacobian.o: 6th_jacobian.c function.h
	mpicc -O3    -c -std=c99 6th_jacobian.c
flow_initialisation.o: flow_initialisation.c function.h
	mpicc -O3    -c -std=c99 flow_initialisation.c
solver.o: solver.c function.h
	mpicc -O3    -c -std=c99 solver.c
roe_average.o: roe_average.c function.h
	mpicc -O3    -c -std=c99 roe_average.c
viscous_terms.o: viscous_terms.c function.h
	mpicc -O3    -c -std=c99 viscous_terms.c
characteristic_plane.o: characteristic_plane.c function.h
	mpicc -O3    -c -std=c99 characteristic_plane.c
tecplot_output.o: tecplot_output.c function.h
	mpicc -O3    -c -std=c99 tecplot_output.c

clean:
		rm -rf *.o SIM gnu_plot.txt *.txt metis*
cleanall:
		rm -rf *.o SIM *.plt gnu_plot.txt *.txt metis* restart_file.neu node_neighbour.neu elem_neighbour.neu restart_file_1.neu restart_file_2.neu *.dat
cleanr:
		rm -rf *.o SIM *.plt gnu_plot.txt *.txt metis* restart_file.neu
