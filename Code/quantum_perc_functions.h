


int find_best_neighbor (int focal_node, int other_node, struct network *G);
void multiple_pair_entanglement_improved (struct entanglement **Q, struct d_search *S, struct network *G, int P);
void single_pair_entanglement_standard(int source, int target, struct entanglement **Q, struct d_search *S, struct network *G, int p, int step);
void single_pair_entanglement_improved (struct entanglement **Q, struct d_search *S, struct network *G, int p);

void print_all_paths_for_equal_equal_weigh_case (struct entanglement **Q, struct network *G, int P);
void print_results_to_screen_compact_improved (struct entanglement **Q, int P, int percolation_type, int log_var);
void test_individual_pairs_improved (struct entanglement **Q, int p, int percolation_type, int log_var, double *results);
///
void simple_bfs (struct d_search *S, struct network *G, int source, int target);
void bfs_random_path (struct d_search *S, int source, int target);
void reset_variables(struct d_search *S);
void use_path (struct d_search *S, struct network *G, struct entanglement **Q, int p, int step);
///
void allocate_memory_search (struct d_search *S, struct network *G);
void deallocate_memory_search (struct d_search *S);
void allocate_memory_entanglement (struct entanglement **Q, int P);
void deallocate_memory_entanglement (struct entanglement **Q, int P);


