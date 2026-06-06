#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include "structs.h"
#include "network_functions.h"
#include "quantum_perc_functions.h"
#include "mt64.h"








void multiple_pair_entanglement_improved (struct entanglement **Q, struct d_search *S, struct network *G, int P)
{

  int p;

  for(p=1;p<=P;p++)
    {
      //standard version//////////////////////
      single_pair_entanglement_improved (Q, S, G, p);
      ///////////////////////////////////////
    }
}



void single_pair_entanglement_standard(int source, int target, struct entanglement **Q, struct d_search *S, struct network *G, int p, int step)
{
  //standard version//////////////////////
  S->control = 1;
  Q[p]->tmp_distil_coeff = 1.0;
  Q[p]-> exit_from_loop = -1;
  while(S->control > 0){
    //simple_bfs (S, G, Q[p]->source, Q[p]->target);
    simple_bfs (S, G, source, target);
    //printf("#Variable control = %d\n",S->control);
    if (Q[p]-> exit_from_loop >  0)
      {
	S->control = 0;
	//printf("#*****Exit from loop******\n\n");
      }
    if(S->control > 0)
      {
	//bfs_random_path (S, Q[p]->source, Q[p]->target);
	bfs_random_path (S, source, target);
	use_path (S, G, Q, p, step);
      }
    reset_variables(S);
  }
  ////////////////////////////////////////
}



void single_pair_entanglement_improved (struct entanglement **Q, struct d_search *S, struct network *G, int p)
{
  int best_source, best_target;
  
  int source = (int)(genrand64_real3()*(double)G->nr_nodes) + 1;
  int target = (int)(genrand64_real3()*(double)G->nr_nodes) + 1;

  while(target == source)
    {
      target +=1;
      if (target > G->nr_nodes) target = 1;
    }

  Q[p]->source = source;
  Q[p]->target = target;

  
  if(G->adjacency[0][source] > G->adjacency[0][target])
    {
      Q[p]->source = target;
      Q[p]->target = source;
      target = Q[p]->target;
      source = Q[p]->source;
    }
  
  
  G->nr_removed_edges = 0;


  //////////////////////////


  if (G->analysis_mode == 3)
    {
      G->lambda = 0.99999;
      G->percolation_type = 2;
      //equal weight case version//////////////////////
      single_pair_entanglement_standard(source, target, Q, S, G, p, 1);
      ///////////////////////////////////////
    }

  if (G->analysis_mode == 0)
    {
      //standard version//////////////////////
      single_pair_entanglement_standard(source, target, Q, S, G, p, 1);
      ///////////////////////////////////////
    }
  
  
  if (G->analysis_mode == 1 || G->analysis_mode == 2)
    {
      //improved version version//////////////////////
      best_source = find_best_neighbor (source, target, G);
      best_target = find_best_neighbor (target, source, G);
      //printf("\n\n\n");
      //printf("#source = %d %d\n", source, best_source);
      //printf("#target = %d %d\n", target, best_target);

      //case 0
      if ((best_source == target) || (best_target == source)){
        //printf("#case 0\n");
	single_pair_entanglement_standard(source, target, Q, S, G, p, 1);
      }
      else{
	//case 1
        if ((source != best_source) && (best_source != best_target) && (best_target != target)){
	  //printf("#case 1\n");
	  single_pair_entanglement_standard(source, best_source, Q, S, G, p, 1);
	  single_pair_entanglement_standard(best_target, target, Q, S, G, p, 3);
	  single_pair_entanglement_standard(best_source, best_target, Q, S, G, p, 2);
	}
	//#case 2
        if ((source != best_source) && (best_source != best_target) && (best_target == target)){
	  //printf("#case 2\n");
	  single_pair_entanglement_standard(source, best_source, Q, S, G, p, 1);
	  single_pair_entanglement_standard(best_source, target, Q, S, G, p, 2);
	}
	//#case 3
        if ((source != best_source) && (best_source == best_target) && (best_target != target)){
	  //printf("#case 3\n");
	  single_pair_entanglement_standard(source, best_source, Q, S, G, p, 1);
	  single_pair_entanglement_standard(best_source, target, Q, S, G, p, 2);
	}
	//#case 4
        if ((source != best_source) && (best_source == best_target) && (best_target == target)){
	  //printf("#case 4\n");
	  single_pair_entanglement_standard(source, target, Q, S, G, p, 1);
	}
        //#case 5
        if ((source == best_source) && (best_source != best_target) && (best_target != target)){
	  //printf("#case 5\n");
	  single_pair_entanglement_standard(source, best_target, Q, S, G, p, 1);
	  single_pair_entanglement_standard(best_target, target, Q, S, G, p, 2);
	}
        //#case 6
        if ((source == best_source) && (best_source != best_target) && (best_target == target)){
	  //printf("#case 6\n");
	  single_pair_entanglement_standard(source, target, Q, S, G, p, 1);
	}
        //#case 7
        if ((source == best_source) && (best_source == best_target) && (best_target != target)){
	  //printf("#case 7\n");
	  single_pair_entanglement_standard(source, target, Q, S, G, p, 1);
	}
        //#case 8 #impossible, source and target are the same nodes
        if ((source == best_source) && (best_source == best_target) && (best_target == target)){
	  //printf("#case 8\n");
	  source = source;
	}
      }
      ////////////////////////////////////////////////
    }


  
  //////////////////////////
  
  restore_edges(G);
}



int find_best_neighbor (int focal_node, int other_node, struct network *G)
{
  int j, m, n, i, control;
  int max_degree = G->adjacency[0][focal_node];
  int best_node = focal_node;



  if (G->analysis_mode == 1){
  
    control = 1;
    for(j=1;j<=G->adjacency[0][focal_node];j++)
      {
	m = G->adjacency[focal_node][j];
	if(G->adjacency[0][m] > max_degree && control > 0)
	  {
	    max_degree = G->adjacency[0][m];
	  best_node = m;
	  }
      if (m == other_node)
	{
	  control = -1;
	  max_degree = G->adjacency[0][m];
	  best_node = m;
	}
      }
  }


  if (G->analysis_mode == 2){
  
    control = 1;
    for(j=1;j<=G->adjacency[0][focal_node];j++)
      {
	m = G->adjacency[focal_node][j];
	if(G->adjacency[0][m] > max_degree && control > 0)
	  {
	    max_degree = G->adjacency[0][m];
	    best_node = m;
	  }
	if (m == other_node)
	  {
	    control = -1;
	    max_degree = G->adjacency[0][m];
	    best_node = m;
	  }


	for(i=1;i<=G->adjacency[0][m];i++)
	  {
	    n = G->adjacency[m][i];
	    if(G->adjacency[0][n] > max_degree && control > 0)
	      {
		max_degree = G->adjacency[0][n];
		best_node = n;
	      }
	    if (n == other_node)
	      {
		control = -1;
		max_degree = G->adjacency[0][n];
		best_node = n;
	      }
	  }
	
	
      }
  }

    
  return best_node;


}




void print_results_to_screen_compact_improved (struct entanglement **Q, int P, int percolation_type, int log_var)
{



  
  double avl1, avc1, ave, total;
  double totalP =(double)P;
  
  avl1 = 0.0;
  avc1 = 0.0;
  total = ave = 0.0;

  
  double res[2], tmp;
  int p;

 
  for(p=1;p<=P;p++)
    {
      
      test_individual_pairs_improved (Q, p, percolation_type, log_var, res);
      if(res[1] > -1)
	{
	  total += 1;
	  ave += res[1];
	}
	  
      if (percolation_type ==1)
	{
	  avl1 += 2.0*(1.0 - res[0]);
	  tmp = from_schmidt_to_concurrence (res[0]);
	  avc1 += tmp;
	}
      
      if (percolation_type ==2)
	{
	  avl1 += res[0];
	  tmp = 1.0 - 0.5*res[0];
	  tmp = from_schmidt_to_concurrence (tmp);
	  avc1 += tmp;
	}


      
    }
  

  avl1 /= totalP;
  avc1 /= totalP;
  ave /= totalP;
  total /= totalP;
  
  
  printf("%g %g %g %g\n",avl1,avc1,ave,total);
  
}





void test_individual_pairs_improved (struct entanglement **Q, int p, int percolation_type, int log_var, double *results)
{

  int i, j, control = 1;
  int max_steps;
  double coefficient_steps[4];
  double final_coeff;
  int intermediate_steps[4][2];

  for(i=0;i<4;i++) for(j=0;j<2;j++) intermediate_steps[i][j] = -1;

  results[0] = -1;
  results[1] = 0;


  if(log_var ==1) printf("\n--- \n\n\n");

  if(Q[p]->nr_paths ==0) control = -1;

  else{

     max_steps = -1;
     coefficient_steps[0] = coefficient_steps[1] = coefficient_steps[2] = coefficient_steps[3] = 1.0; 

      for(i=1;i<=Q[p]->nr_paths;i++){
	results[1] +=  Q[p]->nr_paths-1;
	
	if (log_var == 1) printf("path nr %d :", i);
	if (log_var == 1) for(j=1;j<Q[p]->paths[i][0];j++) printf("%d -> ",Q[p]->paths[i][j]);
	if (log_var == 1) printf("%d : %g %d\n",Q[p]->paths[i][Q[p]->paths[i][0]], Q[p]->conc_paths[i], Q[p]->step_paths[i]);
	
	
	if (percolation_type == 1) coefficient_steps[Q[p]->step_paths[i]] *= from_concurrence_to_schmidt (Q[p]->conc_paths[i]);
	if (percolation_type == 2) coefficient_steps[Q[p]->step_paths[i]] *= 1.0 - Q[p]->conc_paths[i];
	intermediate_steps[Q[p]->step_paths[i]][1] = Q[p]->paths[i][1];
	intermediate_steps[Q[p]->step_paths[i]][0] = Q[p]->paths[i][Q[p]->paths[i][0]];
	
	if (Q[p]->step_paths[i] > max_steps) max_steps = Q[p]->step_paths[i];
	
      }
     
  

      if (log_var == 1) printf("\nsource %d \t target %d\n",Q[p]->source, Q[p]->target);
      
      if (log_var == 1) for(i=1;i<=max_steps;i++) printf("%d -> %d \t",intermediate_steps[i][0], intermediate_steps[i][1]);
	if (log_var == 1) printf("\n");
	if(intermediate_steps[1][0]!=Q[p]->source || intermediate_steps[max_steps][1]!=Q[p]->target) control = -1;
	for(i=1;i<max_steps;i++) if(intermediate_steps[i][1] != intermediate_steps[i+1][0]) control = -1;
  
  }
  


  if(control < 0 )
    {
      if (log_var == 1) printf("No path found!\n");
      if (percolation_type == 1) results[0] = 1.0; 
      if (percolation_type == 2) results[0] = 0.0; 
      results[1] = -1;
      return;
    }

  if(control >0){
    
    final_coeff = 1.0;
    for(i=1;i<=max_steps;i++)
      {
	if (percolation_type == 1) if (coefficient_steps[i] < 0.5) coefficient_steps[i] = 0.5;
	if (percolation_type == 2) coefficient_steps[i] = 1.0 - coefficient_steps[i];
	if (log_var == 1) printf("%d coeff = %g\n", i, coefficient_steps[i]);
	if (percolation_type == 1) final_coeff *= from_schmidt_to_concurrence(coefficient_steps[i]);
	if (percolation_type == 2)  final_coeff *= coefficient_steps[i];
      }
    if (log_var == 1)
      {
	if (percolation_type == 1) printf("final %g\n", from_concurrence_to_schmidt(final_coeff));
	if (percolation_type == 2) printf("final %g\n", final_coeff);
      }
    
    if (percolation_type == 1)
      {
	results[0] = from_concurrence_to_schmidt(final_coeff); //return from_concurrence_to_schmidt(final_coeff);
	if (results[0] < 0.5) results[0] = 0.5;
      }
    if (percolation_type == 2) results[0] = final_coeff; //return final_coeff;


    Q[p]->distil_coeff = results[0];
  
  return;

    
  }


   return;
 
  
  
}





void print_all_paths_for_equal_equal_weigh_case (struct entanglement **Q, struct network *G, int P)
{
  int p, n, m;
  for(p=1;p<=P;p++)
    {
      n = Q[p]->source;
      m = Q[p]->target;
      printf("%d %d %d",G->adjacency[0][n], G->adjacency[0][m], Q[p]->nr_paths);
      for(n=1;n<=Q[p]->nr_paths;n++) printf(" %d",Q[p]->paths[n][0]);
      printf("\n");
    }
}



///////////////////////////



void simple_bfs (struct d_search *S, struct network *G, int source, int target)
{

  int i, j, m, n, d;

  S->vec[0] = 1;
  S->vec[S->vec[0]] = source;
  S->visited[source] = 0; 
  S->dag[0][source] = 0;
  d = 0;
  S->nr_sp[source] = 1;


  S->reset[0] = 1;
  S->reset[1] = source;
  S->control = 0;


  while(S->vec[0]>0 && S->control ==0)
    {

      S->tmp_vec[0] = 0;
      d += 1;

      for(i=1;i<=S->vec[0];i++)
	{
	  n = S->vec[i];
	  for(j=1;j<=G->adjacency[0][n];j++)
	    {
	      m = G->adjacency[n][j];
	      if (m == target) S->control = 2;
	      if(S->visited[m]<0){
		S->tmp_vec[0] +=1;
		S->tmp_vec[S->tmp_vec[0]] = m;
		S->visited[m] = S->visited[n] + 1;
		S->reset[0] += 1;
		S->reset[S->reset[0]] = m;
	      }
	      
	      if (S->visited[m] == d)
		{
		  S->dag[0][m] +=1;
		  S->dag[m][S->dag[0][m]] = n;
		  S->nr_sp[m] += S->nr_sp[n];
		}
	      
	      
	    }
	  
	}
      for(i=0;i<=S->tmp_vec[0];i++) S->vec[i] = S->tmp_vec[i]; 
    }

  /*
  printf("#distance %d\n",d);
  for(i=1;i<=G->nr_nodes;i++)
    {
      printf("%d: ",i);
      for(j=1;j<=S->dag[0][i];j++) printf("%d ",S->dag[i][j]);
      printf("\n");
    }
  */
}



void bfs_random_path (struct d_search *S, int source, int target)
{

   

  int i, n , q, T, found, control;
  
  control = 0;
  S->path[0] = 1;
  S->path[1] = target;

  //printf("s = %d t = %d\n",source,target);
  //printf("--> %d\t",S->path[S->path[0]] );
  
   while(control==0)
    {
      n = S->path[S->path[0]];
      //printf("A) %d %d\n",n,S->dag[0][n]);

      if(S->dag[0][n]>0)
	{
	  //find next node on the path
	  T = 0;
	  for(i=1;i<=S->dag[0][n];i++) T += S->nr_sp[S->dag[n][i]];
	  
	  q = (int)(genrand64_real3()*(double)T) + 1;
	  if (q>T) q=1;

	  //printf("B) T = %d\t q = %d\n",T,q);


	  i = T = found = 0;
	  while(found==0)
	    {	      
	      //printf("C) T = %d\t q = %d\t\t%d %d\n",T,q,i,S->dag[0][n]);
	      i++;
	      T += S->nr_sp[S->dag[n][i]];
	      if(q<=T)
		{
		  //printf("->i = %d\n",i);
		  found = S->dag[n][i];
		}
	    }
	  //printf("f = %d %d\n\n",i,found);
	  /////////////////

	  S->path[0]++;
	  S->path[S->path[0]] = found;
	  //printf("--> %d\t",S->path[S->path[0]]);
	  
	 if(S->path[S->path[0]] == source) control = 1;

	}
      else{
	S->path[0]=0;
	//printf("\tN\n\n");
	return;
      }
    }
   //printf("\tO\n\n");
   return;
}




void use_path (struct d_search *S, struct network *G, struct entanglement **Q, int p, int step)
{
  
  Q[p]->nr_paths ++;
  Q[p]->paths = (int **)realloc(Q[p]->paths, (Q[p]->nr_paths+1)*sizeof(int *));
  Q[p]->paths[Q[p]->nr_paths] = (int *)malloc((S->path[0]+1)*sizeof(int));
  Q[p]->paths[Q[p]->nr_paths][0] = 0;
  Q[p]->conc_paths = (double *)realloc(Q[p]->conc_paths, (Q[p]->nr_paths+1)*sizeof(double));
  Q[p]->step_paths = (int *)realloc(Q[p]->step_paths, (Q[p]->nr_paths+1)*sizeof(int));


  int i, n, m;

  for (i=1;i<=S->path[0];i++)
    {
      n = S->path[i];
      m = Q[p]->paths[Q[p]->nr_paths][0] + 1;
      Q[p]->paths[Q[p]->nr_paths][0] += 1;
      Q[p]->paths[Q[p]->nr_paths][m] = n;

      //delete the edge
      if (i<S->path[0]) delete_edge (G, S->path[i], S->path[i+1]);
     

      //if (i<S->path[0]) printf("%d -> ",n);
    }
  //printf("%d\n", S->path[S->path[0]]);


  double tmp;
  //
  //
  Q[p]->step_paths[Q[p]->nr_paths] = step;
  if(G->percolation_type == 1){
    tmp = (S->path[0]-1)*log(from_schmidt_to_concurrence (G->lambda));
    Q[p]->conc_paths[Q[p]->nr_paths] = exp(tmp);
    Q[p]->tmp_distil_coeff *= from_concurrence_to_schmidt (exp(tmp));
    if(Q[p]->tmp_distil_coeff < 0.5)
      {
	Q[p]->exit_from_loop = 1;
      Q[p]->tmp_distil_coeff = 0.5;
      }
  }
  if(G->percolation_type == 2){
    tmp = (S->path[0]-1)*log(2*(1-G->lambda));
    Q[p]->conc_paths[Q[p]->nr_paths] = exp(tmp);
    Q[p]->tmp_distil_coeff = 1.0;
  }
  //
  //
  
}
  



void reset_variables(struct d_search *S)
{

  int i, n;

  for(i=1;i<=S->path[0];i++) S->path[0] = 0;

  
  for(i=1;i<=S->reset[0];i++)
    {
      n = S->reset[i];
      S->visited[n] = -1;
      S->dag[0][n] = 0;
      S->nr_sp[n] = 0;
    }

}



////////////////////////////



////////////////////////////

void allocate_memory_search (struct d_search *S, struct network *G)
{
  int i;
  //printf("#Number of nodes %d\n", G->nr_nodes); fflush(stdout);
  int N = G->nr_nodes;
  S->N = N;
  S->visited = (int *)malloc((N+1)*sizeof(int));
  S->reset = (int *)malloc((N+1)*sizeof(int));
  S->vec = (int *)malloc((N+1)*sizeof(int));
  S->path = (int *)malloc((N+1)*sizeof(int));
  S->tmp_vec = (int *)malloc((N+1)*sizeof(int));
  for(i=0;i<=N;i++) S->visited[i] = S->reset[i] = S->vec[i] = S->tmp_vec[i] =  0;

  ////
  S->dag = (int **)malloc((N+1)*sizeof(int *));
  S->dag[0] = (int *)malloc((N+1)*sizeof(int));
  S->nr_sp = (int *)malloc((N+1)*sizeof(int));
 
  
  for(i=1;i<=N;i++)
    {
      S->dag[i] = (int *)malloc((G->adjacency[0][i]+1)*sizeof(int));
      S->dag[0][i] = 0;
      S->visited[i] = -1;
      S->nr_sp[i] = 0;
    }
  //////
  S->reset[0] = 0;

}



void deallocate_memory_search (struct d_search *S)
{
  int i;
  free(S->visited);
  for(i=0;i<=S->N;i++) free(S->dag[i]);
  free(S->dag);
  free(S->path);
  free(S->nr_sp);
  free(S->reset);
  free(S->vec);
  free(S->tmp_vec);
  free(S);
}



//
void allocate_memory_entanglement (struct entanglement **Q, int P)
{
  int p;
  for(p=1;p<=P;p++)
    {
      Q[p] = (struct entanglement*)malloc(1*sizeof(struct entanglement));
      Q[p]->source = -1;
      Q[p]->target = -1;
      Q[p]->nr_paths = 0;
      Q[p]->paths = (int **)malloc(1*sizeof(int *));
      Q[p]->paths[0] = (int *)malloc(1*sizeof(int));
      Q[p]->paths[0][0] = 0;
      Q[p]->conc_paths = (double *)malloc(1*sizeof(double));
      Q[p]->conc_paths[0] = 0.0;
      Q[p]->step_paths = (int *)malloc(1*sizeof(int));
      Q[p]->distil_coeff = 1.0;
      Q[p]->tmp_distil_coeff = 1.0;
      Q[p]->exit_from_loop = -1;
    }
}



void deallocate_memory_entanglement (struct entanglement **Q, int P)
{
  int p, i;
  for(p=1;p<=P;p++)
    {
      for(i=0;i<=Q[p]->nr_paths;i++)  if(Q[p]->paths[i] != NULL) free(Q[p]->paths[i]);
      if(Q[p]->paths != NULL) free(Q[p]->paths);
      if(Q[p]->conc_paths != NULL) free(Q[p]->conc_paths);
      if(Q[p]->step_paths != NULL) free(Q[p]->step_paths);
      if(Q[p] != NULL) free(Q[p]);
    }
  if(Q != NULL) free(Q);
}








