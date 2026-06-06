#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include "structs.h"
#include "mt64.h"
#include "network_functions.h"
#include "quantum_perc_functions.h"



int main (int argc, char **argv)
{



  int pid_id= time(NULL) * getpid();
  init_genrand64((unsigned)pid_id);

  
  clock_t start, end;
  double cpu_time_used;


  struct network *G = (struct network*)malloc(1*sizeof(struct network));
  count_nodes (argv[1], G);

  printf("#Network file: %s\n",argv[1]);
  printf("#Number of nodes: %d\n",G->nr_nodes);
  printf("#Number of edges: %d\n",G->nr_edges);

  allocate_memory_real_network (G);

  read_edges (argv[1], G);


  //
  struct d_search *S = (struct d_search*)malloc(1*sizeof(struct d_search));
  allocate_memory_search (S, G);
  //

  double W = atof(argv[2]);
  printf("#W = %g\n",W);
  G->lambda = W;
  

  int P = atoi(argv[3]);
  printf("#P = %d\n",P);
        

  G->percolation_type = atoi(argv[4]);
  if (G->percolation_type == 1) printf("#Concurrence percolation\n");
  if (G->percolation_type == 2) printf("#Entanglement percolation\n");


  //
  struct entanglement **Q = (struct entanglement**)malloc((P+1)*sizeof(struct entanglement*));
  allocate_memory_entanglement (Q, P);
  //


  G->analysis_mode = atoi(argv[5]);
  printf("#A = %d\n",G->analysis_mode);

  //set A = 0 for standard version, h0
  //set A = 1 for h1 version
  //set A = 2 for h2 standard version
  
  

  
  start = clock();
  ////////////////////////////


  multiple_pair_entanglement_improved (Q, S, G, P);


  
  
  
  //////////////////////////
  end = clock();
  cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
  printf("#estimate entanglement : %g\n",cpu_time_used); fflush(stdout);
  /////////////////////////


  if(G->analysis_mode == 0 || G->analysis_mode == 1 || G->analysis_mode == 2){
    print_results_to_screen_compact_improved (Q, P, G->percolation_type, 0);
  }
  if(G->analysis_mode == 3) print_all_paths_for_equal_equal_weigh_case (Q, G, P);


  //
  deallocate_memory_network (G);
  deallocate_memory_search (S);
  deallocate_memory_entanglement (Q, P);
  //


  return 0;
}
