#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include "mt64.h"
#include "structs.h"
#include "network_functions.h"

void count_nodes (char *filename, struct network *G)
{
  FILE *f;
  int i, j, q;
  int N = 0;
  int E = 0;

  f = fopen(filename,"r");
  while(!feof(f))
    {
      q = fscanf(f,"%d %d",&i,&j);
      if (q<=0) goto exit_file;
      if(i>N) N = i;
      if(j>N) N = j;
      E += 1;
    }
 exit_file:
  fclose(f);
  G->nr_nodes  = N;
  G->nr_edges  = E;
}



void allocate_memory_real_network (struct network *G)
{
  int i;
  int N = G->nr_nodes;
  
  //memory allocation
  G->adjacency = (int **)malloc((N+1)*sizeof(int *));
  G->adjacency[0] = (int *)malloc((N+1)*sizeof(int));
  
  for(i=1;i<=N;i++)
    {
      G->adjacency[0][i] = 0;
      G->adjacency[i] =  (int *)malloc(1*sizeof(int));
    }

  G->nr_removed_edges = 0;
  G->removed_edges = (int **)malloc((G->nr_edges+1)*sizeof(int *));
  for(i=1;i<=G->nr_edges;i++)
    {
      G->removed_edges[i] = (int *)malloc(2*sizeof(int));
      G->removed_edges[i][0] = G->removed_edges[i][1] = 0;
    }

}


void read_edges (char *filename, struct network *G)
{
  FILE *f;
  int i, j, q;


  f = fopen(filename,"r");
  while(!feof(f))
    {
      q = fscanf(f,"%d %d",&i,&j);
      if (q<=0) goto exit_file;
      //add connection
      add_edge (G, i, j);
    }
 exit_file:
  fclose(f);

}

/////////////////////////////

void deallocate_memory_network (struct network *G)
{
  int i;
  for(i=0;i<=G->nr_nodes;i++)
    {
      if(G->adjacency[i] != NULL) free(G->adjacency[i]);
    }
  if(G->adjacency != NULL) free(G->adjacency);


  for(i=1;i<=G->nr_edges;i++)  if(G->removed_edges[i] != NULL) free(G->removed_edges[i]);
  if(G->removed_edges != NULL) free(G->removed_edges);


  if(G != NULL) free(G);
}

////

int find_node (int i, int j, int **adjacency)
{
  int k;
  for(k=1;k<=adjacency[0][i];k++) if(adjacency[i][k]==j) return k;
  return -1;
}


////


void print_network_structure (struct network *G)
{
  int n, i;
  for(n=1;n<=G->nr_nodes;n++)
    {
      printf("-->adj n = %d (%d)",n,G->adjacency[0][n]);
      for(i=1;i<=G->adjacency[0][n];i++) printf("\t%d ",G->adjacency[n][i]);
      printf("\n");
    }
  
}


void add_edge (struct network *G, int n, int m)
{
 //
      G->adjacency[0][n]++;
      G->adjacency[n]= (int *)realloc(G->adjacency[n], (G->adjacency[0][n]+1)*sizeof(int));
      G->adjacency[n][G->adjacency[0][n]] = m;
      G->adjacency[0][m]++;
      G->adjacency[m]= (int *)realloc(G->adjacency[m], (G->adjacency[0][m]+1)*sizeof(int));
      G->adjacency[m][G->adjacency[0][m]] = n;
      //

}

//edges are not relabeled
void delete_edge (struct network *G, int n, int m)
{

  //
  int en, em;
  
  en = find_node (n, m, G->adjacency);
  em = find_node (m, n, G->adjacency);
  
  //remove edge
  G->nr_removed_edges ++;
  G->removed_edges[G->nr_removed_edges][0] = n;
  G->removed_edges[G->nr_removed_edges][1] = m;
  
  //
  G->adjacency[n][en] = G->adjacency[n][G->adjacency[0][n]];
  G->adjacency[0][n]--;
    
  //
  G->adjacency[m][em] = G->adjacency[m][G->adjacency[0][m]];
  G->adjacency[0][m]--;
   
}


void restore_edges(struct network *G)
{

  int e, n, m;
  for(e=1;e<=G->nr_removed_edges;e++)
    {
      n = G->removed_edges[e][0];
      m = G->removed_edges[e][1];
      add_edge (G, n, m);
    }
  
  G->nr_removed_edges = 0;
}






//////////////////////


double from_schmidt_to_concurrence (double lambda)
{
  double c = 2.0 * sqrt(lambda) * sqrt(1.0 - lambda);
  if (c > 1) c = 1;
  if (c < 0) c = 0;
  return c;
}

double from_concurrence_to_schmidt (double c)
{
  double lambda = 0.5 + 0.5 * sqrt(1.0 - c*c);
  return lambda;
}
     










