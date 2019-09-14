ifeq ($(shell uname),Darwin)
	CC= clang++ 
	CCFLAGS = -Wall -g -std=c++11 -O3 -I /opt/local/include -DArcRemoval -DWarmStart -DPathFollowing #-DWithoutRounding
else
	CC= g++
	CCFLAGS = -Wall -g -std=c++11 -O3 -DArcRemoval -DWarmStart -DPathFollowing #-DWithoutRounding
endif


CCFLAGS_LEMON = -Wall -g -std=c++11 -lemon -I /home/ruben/lemon/include -L /home/ruben/lemon/lib -O3 -DNDEBUG #-DVERBOSE
CCFLAGS_GMP = -Wall -g -std=c++11 -lgmpxx -lgmp -O3 -DNDEBUG #-DVERBOSE
CCFLAGS_delta_wye = -Wall -g -std=c++11 -O3 -DArcRemoval -DWarmStart  #-DNDEBUG  #-DVERBOSE

default: \
	min_cost_flow_read_graph_path_following \
	
min_cost_flow_random_gen_gmp: min_cost_flow_random_gen_gmp.cpp  graph.h potential_reduction.h electrical_flow_problem.h min_cost_flow.h rational_functions.h unbounded_integers.h gmp_helper.h random.h
	$(CC) $(CCFLAGS_GMP) min_cost_flow_random_gen_gmp.cpp -o min_cost_flow_random_gen_gmp

min_cost_flow_random_gen_gmp_with_delta_wye: min_cost_flow_random_gen_gmp_with_delta_wye.cpp  graph.h potential_reduction.h electrical_flow_problem.h min_cost_flow.h rational_functions.h unbounded_integers.h gmp_helper.h random.h
	$(CC) $(CCFLAGS_GMP) min_cost_flow_random_gen_gmp_with_delta_wye.cpp -o min_cost_flow_random_gen_gmp_with_delta_wye

min_cost_flow_read_graph_gmp: min_cost_flow_read_graph_gmp.cpp  graph.h potential_reduction.h electrical_flow_problem.h min_cost_flow.h rational_functions.h unbounded_integers.h gmp_helper.h random.h
	$(CC) $(CCFLAGS_GMP) min_cost_flow_read_graph_gmp.cpp -o min_cost_flow_read_graph_gmp

lemon_capacity_scaling: lemon_capacity_scaling.cc
			$(CC) $(CCFLAGS_LEMON) lemon_capacity_scaling.cc -o lemon_capacity_scaling
lemon_cycle_cancelling: lemon_cycle_cancelling.cc
			$(CC) $(CCFLAGS_LEMON) lemon_cycle_cancelling.cc -o lemon_cycle_cancelling
lemon_network_simplex: lemon_network_simplex.cc
	$(CC) $(CCFLAGS_LEMON) lemon_network_simplex.cc -o lemon_network_simplex

electrical_flow: electrical_flow.cpp graph.h unbounded_integers.h
	$(CC) $(CCFLAGS) electrical_flow.cpp -o electrical_flow

min_cost_flow_random_gen: min_cost_flow_random_gen.cpp  graph.h potential_reduction.h electrical_flow_problem.h min_cost_flow.h unbounded_integers.h random.h
	$(CC) $(CCFLAGS) min_cost_flow_random_gen.cpp -o min_cost_flow_random_gen

min_cost_flow_random_gen_with_delta_wye: min_cost_flow_random_gen_with_delta_wye.cpp  graph.h potential_reduction.h electrical_flow_problem.h min_cost_flow.h unbounded_integers.h  potential_reduction_DS.h electrical_flow_problem_DS.h random.h
	$(CC) $(CCFLAGS_delta_wye) min_cost_flow_random_gen_with_delta_wye.cpp -o min_cost_flow_random_gen_with_delta_wye

min_cost_flow_read_graph: min_cost_flow_read_graph.cpp  graph.h potential_reduction_DS.h electrical_flow_problem_DS.h min_cost_flow.h rational_functions.h unbounded_integers.h random.h
	$(CC) $(CCFLAGS) min_cost_flow_read_graph.cpp -o min_cost_flow_read_graph

min_cost_flow_read_graph_with_delta_wye: min_cost_flow_read_graph_with_delta_wye.cpp  graph.h potential_reduction.h electrical_flow_problem.h min_cost_flow.h rational_functions.h unbounded_integers.h random.h potential_reduction_DS.h electrical_flow_problem_DS.h
	$(CC) $(CCFLAGS_delta_wye) min_cost_flow_read_graph_with_delta_wye.cpp -o min_cost_flow_read_graph_with_delta_wye

min_cost_flow_read_graph_path_following: min_cost_flow_read_graph_path_following.cpp electrical_flow_problem_DS.h  graph.h potential_reduction_path_following.h electrical_flow_problem.h min_cost_flow.h unbounded_integers.h random.h
	$(CC) $(CCFLAGS) min_cost_flow_read_graph_path_following.cpp -o min_cost_flow_read_graph_path_following

min_cost_flow_random_gen_path_following: min_cost_flow_random_gen_path_following.cpp electrical_flow_problem_unrounded.h  graph.h potential_reduction_path_following.h electrical_flow_problem.h min_cost_flow.h unbounded_integers.h random.h
	$(CC) $(CCFLAGS) min_cost_flow_random_gen_path_following.cpp -o min_cost_flow_random_gen_path_following

print_graphlp: print_graphlp.cpp graph.h potential_reduction.h electrical_flow_problem.h unbounded_integers.h
	$(CC) $(CCFLAGS) print_graphlp.cpp -o print_graphlp

min_cost_flow_crossover: min_cost_flow_crossover.cpp  simple_solver.h graph.h potential_reduction.h electrical_flow_problem.h min_cost_flow.h min_cost_crossover.h unbounded_integers.h
	$(CC) $(CCFLAGS) min_cost_flow_crossover.cpp -o min_cost_flow_crossover

min_cost_flow_random_gen_sspvariant: min_cost_flow_random_gen_sspvariant.cpp electrical_flow_problem_unrounded.h  graph.h potential_reduction_path_following.h electrical_flow_problem.h min_cost_flow.h unbounded_integers.h random.h
	$(CC) $(CCFLAGS) min_cost_flow_random_gen_sspvariant.cpp -o min_cost_flow_random_gen_sspvariant

pgm2graph: pgm2graph.cc gridflow.h bucketqueue.h
	$(CC) $(CCFLAGS) pgm2graph.cc -o pgm2graph

recursion: recursion.cc gridflow.h bucketqueue.h
	$(CC) $(CCFLAGS) recursion.cc -o recursion

recursion-variant: recursion.cc gridflow.h bucketqueue.h
	$(CC) $(CCFLAGS) -DSSP_VARIANT recursion.cc -o recursion-variant

recursion2threads: recursion2threads.cc gridflow.h bucketqueue.h
	$(CC) $(CCFLAGS) recursion2threads.cc -o recursion2threads -lpthread

convert_to_dimacs: recursion.cc gridflow.h bucketqueue.h
	$(CC) $(CCFLAGS) -DONLY_CONVERT_TO_DIMACS recursion.cc -o convert_to_dimacs

VALUES = 1 2 4 8 16 32 64 128 256 512 1024
convert_to_dimacs_all: recursion.cc gridflow.h bucketqueue.h
	$(foreach N,$(VALUES),$(CC) $(CCFLAGS) -DONLY_CONVERT_TO_DIMACS -DINSTANCE_CONVERSION_FACTOR=$(N) recursion.cc -o convert_to_dimacs_$(N);)

convert_to_dimacs: recursion.cc gridflow.h bucketqueue.h
	$(CC) $(CCFLAGS) -DONLY_CONVERT_TO_DIMACS recursion.cc -o convert_to_dimacs

convert_to_dimacs: recursion.cc gridflow.h bucketqueue.h
	$(CC) $(CCFLAGS) -DONLY_CONVERT_TO_DIMACS recursion.cc -o convert_to_dimacs

convert_to_dimacs: recursion.cc gridflow.h bucketqueue.h
	$(CC) $(CCFLAGS) -DONLY_CONVERT_TO_DIMACS recursion.cc -o convert_to_dimacs

convert_to_dimacs: recursion.cc gridflow.h bucketqueue.h
	$(CC) $(CCFLAGS) -DONLY_CONVERT_TO_DIMACS recursion.cc -o convert_to_dimacs

output_graph_as_dimacs: min_cost_flow_generic.cpp graph.h random.h min_cost_flow_sspvariant.h min_cost_flow_sspvariant_default.h min_cost_flow_sspvariant_apex_grid.h ppmimage.h unbounded_integers.h bucketqueue.h
	$(CC) $(CCFLAGS) output_graph_as_dimacs.cpp -o output_graph_as_dimacs

exec: electrical_flow
	./electrical_flow


clean:
	rm -f electrical_flow \
	      min_cost_flow_random_gen \
	      min_cost_flow_read_graph \
	      min_cost_flow_read_graph_path_following
