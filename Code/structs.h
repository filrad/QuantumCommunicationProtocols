struct entanglement{
  int source;
  int target;
  int nr_paths;
  int **paths;
  double *conc_paths;
  double distil_coeff;
  double tmp_distil_coeff;
  int *step_paths;
  int exit_from_loop;
};


struct network{
  int nr_nodes;
  int nr_edges;
  int **adjacency; //adjacency matrix
  int nr_removed_edges;
  int **removed_edges;
  //entangled pairs
  int ent_pairs;
  int percolation_type;
  double lambda;
  int analysis_mode;
};


struct d_search{
  int N;
  int *visited;
  int **dag;
  int *nr_sp;
  int *reset;
  int *vec;
  int *tmp_vec;
  int control;
  int *path;
};

