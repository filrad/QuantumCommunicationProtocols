//
void count_nodes (char *filename,struct network *G );
void allocate_memory_real_network (struct network *G);
void deallocate_memory_network (struct network *G);
void read_edges (char *filename, struct network *G);

//
int find_node (int i, int j, int **adjacency);
void delete_edge (struct network *G, int n, int m);
void add_edge (struct network *G, int n, int m);
void restore_edges(struct network *G);
void print_network_structure (struct network *G);

//
double from_schmidt_to_concurrence (double lambda);
double from_concurrence_to_schmidt (double c);



