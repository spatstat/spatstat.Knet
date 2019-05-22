/*

   netKcode.h

   Function definitions for linear K function
   -- definitions that depend on options --

   This code is #included multiple times in *.c files
   with different values of the macro 'INHOM'
*/

#ifndef NETWORKDEF_H_
#include "networkdef.h"
#endif

#undef Crash
#ifndef INHOM
   /* Homogeneous K-function
      Data type 'Crash' is defined as 'MultipleCrash'
      (includes multiplicity but not lambda)
      Function names are not modified
    */
#define Crash MultipleCrash
#define FN(BASE) BASE
#else
   /* Inhomogeneous K-function
      Data type 'Crash' is defined as 'WeightedCrash'
      (includes multiplicity and lambda)
      Function names are prefixed by "I_"
    */
#define Crash WeightedCrash
#define FN(BASE) I_ ## BASE
#endif

#define OK(BASE) BASE

/**************  Prototypes for functions in this file *****************/

int FN(graph_building_with_crash)(Graph *g,
				  int *no_of_vertices,
				  int *no_of_edges,
				  int *no_of_crashes,
				  int *crash_seg,
				  double *crash_tp,
				  int *crash_freq,
#ifdef INHOM
				  double *crash_lambda,
#endif
				  int *vert_id1,
				  int *vert_id2,
				  double *edge_length);

int FN(k_function_v1)(Graph *graph,
		      PathVertex *start,
		      double *distance,
		      double *tme_up_lt_vec,
		      int *m_val_vec,
		      int *vec_len,
		      double *inv_mvals,
		      int *source_cr_freq,
#ifdef INHOM
		      double *source_cr_lambda,
#endif
		      double *MAX_Distance,
		      double *MIN_Distance,
		      int *no_of_distance);

int FN(sum_of_inv_mvals_on_part_of_edge_v1)(PathVertex *adj_vert,
					    double *pth_vert_d,
					    double *dist_ratio,
					    double *tme_up_lt_vec,
					    int *m_val_vec,
					    int *vec_len,
					    double *inv_mvals,
					    int *source_cr_freq,
#ifdef INHOM
					    double *source_cr_lambda,
#endif
					    double *MAX_Distance,
					    double *MIN_Distance,
					    int *no_of_distance);

int FN(sum_of_inv_mvals_for_all_pts_on_edge_v1)(PathVertex *adj_vert,
						double *pth_vert_d,
						double *tme_up_lt_vec,
						int *m_val_vec,
						int *vec_len,
						double *inv_mvals,
						int *source_cr_freq,
#ifdef INHOM
						double *source_cr_lambda,
#endif
						double *MAX_Distance,
						double *MIN_Distance,
						int *no_of_distance);

void FN(createGraphNet)(int *no_of_vertices,
		    int *no_of_edges,
		    int *no_of_crashes,
		    int *crash_seg,
		    double *crash_tp,
		    int *crash_freq,
#ifdef INHOM
		    double *crash_lambda,
#endif
		    int *vert_id1,
		    int *vert_id2,
		    double *edge_length,
		    double *MAX_Distance,
		    double *MIN_Distance,
		    int *no_of_distance,
		    double *max_r,
		    int *verboseIter,
		    double *K_r);

int FN(restricted_shortest_v2)(Graph *graph,
			       const PathVertex *start,
			       List *paths,
			       Graph *exhaustiveTree,
			       double Rmax);

int FN(insert_edge_in_subnet)(Graph *g,PathVertex *p,PathVertex *q);

int FN(crash_compare)(Crash *crash1,Crash *crash2);

int FN(ord_list_ins_next)(List *ord_list,void *data);

int FN(ins_crsh_pthvrtx_list)(PathVertex *pv,PathVertex *pv_rev,Crash *cr_lst,
                          int srch_id,int *no_of_crashes);

int FN(crash_point_in_graph_as_vertex)(Graph *graph,Crash *crashPt,PathVertex *adj_vert1, PathVertex *adj_vert2,int *no_of_vertices);

int FN(copy_crash_list_v2)(List *copy_from,List *copy_to);

int FN(copy_crash_list_rev_v2)(List *copy_from,List *copy_to);

void FN(store_edge_before_deleting)(PathVertex *pth_vertex,
				    PathVertex *adj_vertex,
				    PathVertex *adj_vert1,
				    PathVertex *adj_vert2);

int FN(delete_single_crash_from_crashlist)(List *crlist, Crash *cr);


int FN(extended_sh_path_tree_restricted)(Graph *gLink,
					 Graph *gTree,
					 List *Vlist,
					 List *newPath,
					 double R,
					 int *no_of_vertices,
					 int *no_of_edges);

int FN(break_crash_list_into_two_lists_v2)(List *cr_list1,
					   List *cr_list2,
					   List *crash_list,
					   double dist1,
					   double edgeLength,
					   int e_id1,
					   int e_id2);

int FN(break_crash_list_into_one_list)(List *cr_list,
				       List *crash_list,
				       double dist,
				       int e_id);

int FN(break_crash_list_before_max_dist)(List *cr_list,
					 List *crash_list,
					 double dist,
					 double edgeLength,
					 int e_id);

int FN(spTree_restricted_v1)(Graph *graph,
			     Graph *graphTree,
			     List *P,
			     PathVertex *start);

/******** Declarations of functions in this file **********************/

/********************  ins_crsh_pthvrtx_list  ***************************/

int FN(ins_crsh_pthvrtx_list)(PathVertex *pv,
			      PathVertex *pv_rev,
			      Crash *cr_lst,
			      int srch_id,
			      int *no_of_crashes){
int     id = srch_id,
        retval;

Crash   cr = cr_lst[id];

int Ncr = *no_of_crashes;

while(pv->edgeID == cr.edgeId){
  /** Allocate space for crash structure to be included
      in the crashList of the PathVertex **/
  Crash *crash,*crash_rev;
  crash = (Crash *)Calloc(1,Crash);
  /*if(crash==NULL){Rprintf("malloc is not working !!");return -1;}*/

  crash->tp        = cr.tp;
  crash->edgeId    = cr.edgeId;
  crash->frequency = cr.frequency;
#ifdef INHOM
  crash->lambda    = cr.lambda;
#endif

  retval = FN(ord_list_ins_next)(&pv->crashList,(void *)crash);
  if(retval != 0){
    Rprintf("ord_list_ins_next did not work in ins_crsh_pthvrtx_list!\n");
    return -1;
  }

  crash_rev = (Crash *)Calloc(1,Crash);
  /*if(crash_rev==NULL){Rprintf("malloc is not working !!");return -1;}*/

  crash_rev->tp        = (1 - cr.tp);
  crash_rev->edgeId    = cr.edgeId;
  crash_rev->frequency = cr.frequency;
#ifdef INHOM
  crash_rev->lambda    = cr.lambda;
#endif

  retval = FN(ord_list_ins_next)(&pv_rev->crashList,(void *)crash_rev);
  if(retval != 0){
    Rprintf("ord_list_ins_next did not work in ins_crsh_pthvrtx_list!\n");
    return -1;
  }

  ++id;
  if(id < Ncr){
    cr = cr_lst[id];
  } else {
    break;
  }
 }
 return id;
}

/*****************  graph_building_with_crash  ***************************/

int FN(graph_building_with_crash)(Graph *g,
			      int *no_of_vertices,
			      int *no_of_edges,
                              int *no_of_crashes,
			      int *crash_seg,
			      double *crash_tp,
                              int *crash_freq,
#ifdef INHOM
			      double *crash_lambda,
#endif
			      int *vert_id1,
			      int *vert_id2,
                              double *edge_length){
  /***********************************************************************
   *                       Include Vertices in the Graph                 *
   * Since we're allocating memory here for every vertex, we do not need *
   * to allocate any memory while creating the adjacency lists for these *
   * vertices.                                                           *
   ***********************************************************************/
  int Nv = *no_of_vertices;
  int Ne = *no_of_edges;
  int Ncr = *no_of_crashes;

  int i;
  for(i=0;i<Nv;++i){
    /********************************************************
     *      Memory is allocated for storing vertex ids      *
     *******************************************************/
    PathVertex *pv_ptr;

    pv_ptr = (PathVertex *) Calloc(1,PathVertex);
    /*if(pv_ptr == NULL){Rprintf("malloc is not working !!");return -1;}*/

    pv_ptr->data = Calloc(1,int);

    *((int *)(pv_ptr->data)) = (i+1);

    list_init(&pv_ptr->crashList,NULL);

    int retval_vertex = OK(graph_ins_vertex)(g,pv_ptr); /** graph_ins_vertex() **/

    if(retval_vertex != 0){
      Rprintf("Vertex insertion has failed!\n");
      return -1;
    }
  }
  /**********************************************************************
   *                Build a list of crashes                             *
   *   Assume that the crash data is already sorted                     *
   *   in nondecreasing order of edge id.                               *
   **********************************************************************/

  Crash	*crash_points;

  /** Allocate space in the heap section of the memory
      for storing crash events **/
  crash_points = (Crash *) Calloc((*no_of_crashes),Crash);
  /*if(crash_points == NULL){Rprintf("malloc is not working !!");return -1;}*/

  /** Get "edgeId" ,"tp" and "frequency" values for every crash
      and store them in "crash_points" **/
  int crash_count;
  for(crash_count = 0; crash_count < Ncr; ++crash_count){
    crash_points[crash_count].edgeId = crash_seg[crash_count];
    crash_points[crash_count].tp = crash_tp[crash_count];
    crash_points[crash_count].frequency = crash_freq[crash_count];
#ifdef INHOM
    crash_points[crash_count].lambda = crash_lambda[crash_count];
#endif
  }

  /*****************************************************************
   *            INCLUDE EDGES IN THE GRAPH                         *
   *****************************************************************/
  int search_id = 0,
    edge_count;

  for(edge_count=0; edge_count < Ne; ++edge_count){

    int *vert1 = (int *)Calloc(1,int);/*if(vert1==NULL){Rprintf("malloc is not working !!");return -1;}*/
    int *vert2 = (int *)Calloc(1,int);/*if(vert2==NULL){Rprintf("malloc is not working !!");return -1;}*/
    double edgeWeight;

    *vert1 = vert_id1[edge_count];
    *vert2 = vert_id2[edge_count];
    edgeWeight = edge_length[edge_count];

    PathVertex      pv_ptr1;
    pv_ptr1.data = (void *)vert1;

    PathVertex *pv_ptr2 = (PathVertex *) Calloc(1,PathVertex);
    /*if (pv_ptr2 == NULL){Rprintf("malloc is not working !!");return -1;}*/

    pv_ptr2->data = (void *)(vert2);
    pv_ptr2->weight = (edgeWeight);
    pv_ptr2->edgeID = (edge_count+1);

    OK(list_init)(&pv_ptr2->crashList,&destroy_crash);

    PathVertex  pv_ptr1_rev;
    pv_ptr1_rev.data = (void *)vert2;

    PathVertex *pv_ptr2_rev = (PathVertex *) Calloc(1,PathVertex);
    /*if(pv_ptr2_rev == NULL){Rprintf("malloc is not working !!");return -1;}*/

    pv_ptr2_rev->data = (void *)(vert1);
    pv_ptr2_rev->weight = (edgeWeight);
    pv_ptr2_rev->edgeID = (edge_count+1);

    OK(list_init)(&pv_ptr2_rev->crashList,&destroy_crash);

    if(search_id < Ncr){
      search_id = FN(ins_crsh_pthvrtx_list)(pv_ptr2,
					pv_ptr2_rev,
					crash_points,
					search_id,
					no_of_crashes);
    }

    if(OK(graph_ins_edge)(g, &pv_ptr1, pv_ptr2) !=0){
      Rprintf("Edge insertion has failed!\n");
      return -1;
    }

    if(OK(graph_ins_edge)(g, &pv_ptr1_rev, pv_ptr2_rev) !=0){
      Rprintf("Edge insertion has failed!\n");
      return -1;
    }

  }

  Free(crash_points);

  return 0;
}

/********************* copy_crash_list_v2 *********************************/
int FN(copy_crash_list_v2)(List *copy_from,
		       List *copy_to){
  /** Number of crashes in the list to be copied **/
  int crash_no = copy_from->size;
  if(crash_no==0) return 0;   /** If the list is empty then return **/
  /** Run a for loop to go over each list element
      and copy that to the newly created list **/
  ListElmt *element;
  for(element=list_head(copy_from);element!=NULL;element=list_next(element)){
    Crash   *cr;
    cr=(Crash *)(list_data(element));
    Crash   *crash; /** Allocate memory for the crash **/
    crash = (Crash *) Calloc(1,Crash);
    /*if(crash == NULL){Rprintf("malloc is not working !!");return -1;}*/
    crash->edgeId = cr->edgeId;
    crash->tp     = cr->tp;
    crash->frequency = cr->frequency;
#ifdef INHOM
    crash->lambda = cr->lambda;
#endif
    FN(ord_list_ins_next)(copy_to,(void *)(crash));
  }
  return 0;
}

/**************** copy_crash_list_rev_v2 *********************************/

int FN(copy_crash_list_rev_v2)(List *copy_from,
			   List *copy_to){
  /** Number of crashes in the list to be copied **/
  int crash_no = copy_from->size;
  if(crash_no==0) return 0;  /** If the list is empty then return **/

  /** Run a for loop to go over each list element
      and copy that to the newly created list **/
  ListElmt *element;
  for(element=list_head(copy_from);element!=NULL;element=list_next(element)){
    Crash   *cr;
    cr=(Crash *)(list_data(element));

    Crash   *crash; /** Allocate memory for the crash **/
    crash = (Crash *)Calloc(1,Crash); /*if(crash == NULL){Rprintf("malloc is not working !!");return -1;}*/
    crash->edgeId = cr->edgeId;
    crash->tp     = (1-cr->tp);
    crash->frequency = cr->frequency;
#ifdef INHOM
    crash->lambda = cr->lambda;
#endif
    FN(ord_list_ins_next)(copy_to,(void *)(crash));
  }
  return 0;
}

/********** break_crash_list_into_two_lists_rev_v2 **************************/

int FN(break_crash_list_into_two_lists_rev_v2)(List *cr_list1,
					   List *cr_list2,
					   List *crash_list,
					   double dist1,
					   double edgeLength,
                                           int e_id1,
					   int e_id2){
  /*
   * This function is used in the computation of extended shortest path tree.
   * The edges missing in the shortest path tree are needed to be broken
   * into two new linking edges. Consequently, the crashList on that edge
   * is divided into two separate crash lists which eventually become the
   * crashLists of newly constructed linking edges.
   *
   *  cr_list1 and cr_list2 are lists which are defined outside this function
   *  (memory is allocated in the heap) and initialized.
   *   The references of these lists are passed to this function so that
   *   we can include the crashes on the two linking edges into these lists.
   *
   *  Suppose v1---->v2 is broken into two separate linking edges:
   *      (1) v1---->b1 and (2) v2---->b2
   *  Then cr_list1 and cr_list2 will be crash lists for the two linking edges
   *
   * The input crash_list is the crashList of the edge v1---->v2.
   * dist1 is the distance between v1 and b1.
   * edgeLength is the edge length of the edge v1---->v2.
   ................................................................... */

  int crash_no = crash_list->size;
  /** Number of crashes in the inputted crash list **/

  if(crash_no==0) {return 0;}
  /** If the list is empty then return two empty lists **/

  /** Define variables in order to copy the data
      from the given list to the newly generated list **/
  ListElmt *element;
  Crash *cr;

  double line_prop = (dist1/edgeLength);
  double dist2 = edgeLength - dist1;
  double dist_upto_crash;

  /** Run a for loop to go over each list element
      and allocate crashes to appropriate lists **/
  for(element=list_head(crash_list);element!=NULL;element=list_next(element)){
    cr=(Crash *)(list_data(element));
    if((cr->tp) <= line_prop){
      Crash *crash1 = (Crash *)Calloc(1,Crash);
      /*if (crash1==NULL){Rprintf("malloc is not working !!");return -1;}*/
      crash1->edgeId = e_id1;
      dist_upto_crash = (cr->tp)*edgeLength;
      crash1->tp = (1 - ((dist_upto_crash)/dist1));
      crash1->frequency = cr->frequency;
#ifdef INHOM
      crash1->lambda = cr->lambda;
#endif
      if(FN(ord_list_ins_next)(cr_list1,(void *)(crash1)) != 0) return -1;
    }
    else{
      Crash *crash2 = (Crash *) Calloc(1,Crash);
      /*if (crash2==NULL){Rprintf("malloc is not working !!");return -1;}*/
      crash2->edgeId = e_id2;
      dist_upto_crash = (cr->tp)*edgeLength;
      crash2->tp = (1 - (dist2-(dist_upto_crash-dist1))/dist2);
      crash2->frequency = cr->frequency;
#ifdef INHOM
      crash2->lambda = cr->lambda;
#endif
      if(FN(ord_list_ins_next)(cr_list2,(void *)(crash2)) != 0) return -1;
    }
  }
  return 0;
}

/******************* break_crash_list_into_one_list ****************/
int FN(break_crash_list_into_one_list)(List *cr_list,
				   List *crash_list,
				   double dist,
				   int e_id)
{
  int crash_no = crash_list->size;

  if(crash_no==0) return 0;

  ListElmt *element;

  Crash *cr;

  if(dist == 0){
    /** This means dist1 = 0 i.e. there is no linking edge from pth_vertex.
	The entire crash list will be copied (in reverse order)
	inside the linking edge adj_vertex---->brk_pt **/

    for(element=list_head(crash_list);
	element!=NULL;
	element=list_next(element)){

      cr = (Crash *)(list_data(element));
      Crash *crash = (Crash *)Calloc(1,Crash);
      /*if (crash==NULL){Rprintf("malloc is not working !!");return -1;}*/
      crash->edgeId = e_id;
      crash->tp = (1 - (cr->tp));
      crash->frequency = cr->frequency;
#ifdef INHOM
      crash->lambda = cr->lambda;
#endif
      if(FN(ord_list_ins_next)(cr_list,(void *)(crash)) != 0) return -1;
    }
  } else{
    /** This means dist2 = 0 i.e. there is no linking edge from adj_vertex.
	The entire crash list will be copied
	inside the linking edge pth_vertex---->brk_pt **/
    for(element=list_head(crash_list);
	element!=NULL;
	element=list_next(element)){
      cr = (Crash *)(list_data(element));
      Crash *crash = (Crash *)Calloc(1,Crash);
      /*if (crash==NULL){Rprintf("malloc is not working !!");return -1;}*/
      crash->edgeId = e_id;
      crash->tp = (cr->tp);
      crash->frequency = cr->frequency;
#ifdef INHOM
      crash->lambda = cr->lambda;
#endif
      if(FN(ord_list_ins_next)(cr_list,(void *)(crash)) != 0) return -1;
    }
  }
  return 0;
}

/**************** break_crash_list_into_two_lists_v2 *******************/

int FN(break_crash_list_into_two_lists_v2)(List *cr_list1,
				       List *cr_list2,
				       List *crash_list,
				       double dist1,
				       double edgeLength,
				       int e_id1,
				       int e_id2){
  /** This function is used in the computation of extended shortest path tree.
      The edges missing in the shortest path tree are needed to be broken
      into two new linking edges. Consequently, the crashList on that edge
      is divided into two separate crash lists which eventually become the
      crashLists of newly constructed linking edges.

      cr_list1 and cr_list2 are lists which are defined outside this function
      (memory is allocated in the heap) and initialized.
      The references of these lists are passed to this function so that
      we can include the crashes on the two linking edges into these lists.

      Suppose v1---->v2 is broken into two separate linking edges:
      (1) v1---->b1 and (2) v2---->b2

      Then cr_list1 and cr_list2 will be crash lists for the two linking edges.

      The input crash_list is the crashList of the edge v1---->v2.
      dist1 is the distance between v1 and b1.
      edgeLength is the edge length of the edge v1---->v2.
  **/

  int crash_no = crash_list->size;
  /** Number of crashes in the inputted crash list **/

  if(crash_no==0) {return 0;}
  /** If the list is empty then return two empty lists **/

  /** Define variables in order to copy the data
      from the given list to the newly generated list **/
  ListElmt *element;
  Crash *cr;

  double line_prop = (dist1/edgeLength);
  double dist2 = edgeLength - dist1;
  double dist_upto_crash;

  /** Run a for loop to go over each list element
      and allocate crashes to appropriate lists **/
  for(element=list_head(crash_list);element!=NULL;element=list_next(element)){
    cr=(Crash *)(list_data(element));
    if((cr->tp) <= line_prop){
      Crash *crash1 = (Crash *)Calloc(1,Crash);
      /*if(crash1==NULL){Rprintf("malloc is not working !!");return -1;}*/
      crash1->edgeId = e_id1;
      dist_upto_crash = (cr->tp)*edgeLength;
      crash1->tp = (dist_upto_crash)/dist1;
      crash1->frequency = cr->frequency;
#ifdef INHOM
      crash1->lambda = cr->lambda;
#endif
      if(FN(ord_list_ins_next)(cr_list1,(void *)(crash1)) != 0) return -1;
      /** list_ins_next uses malloc to allocate memory for ListElmt. **/
    } else {
      Crash *crash2 = (Crash *)Calloc(1,Crash);
      /*if(crash2==NULL){Rprintf("malloc is not working !!");return -1;}*/
      crash2->edgeId = e_id2;
      dist_upto_crash = (cr->tp)*edgeLength;
      crash2->tp = (dist2-(dist_upto_crash-dist1))/dist2;
      crash2->frequency = cr->frequency;
#ifdef INHOM
      crash2->lambda = cr->lambda;
#endif
      if(FN(ord_list_ins_next)(cr_list2,(void *)(crash2)) != 0) return -1;
    }
  }

  return 0;
}

/************** break_crash_list_before_max_dist *************************/

int FN(break_crash_list_before_max_dist)(List *cr_list,
				     List *crash_list,
				     double dist,
				     double edgeLength,
				     int e_id){

  int crash_no = crash_list->size;
  /** Number of crashes in the inputted crash list **/

  if(crash_no==0) {return 0;}
  /** If the list is empty then return two empty lists **/

  /** Define variables in order to copy the data
      from the given list to the newly generated list **/
  ListElmt *element;
  Crash *cr;

  double line_prop = (dist/edgeLength);

  double dist_upto_crash;

  /** Run a for loop to go over each list element
      and allocate crashes to appropriate lists **/
  for(element=list_head(crash_list);
      element!=NULL;
      element=list_next(element)){
    cr=(Crash *)(list_data(element));
    if((cr->tp) <= line_prop){
      Crash *crash1 = (Crash *)Calloc(1,Crash);
      /*if (crash1==NULL){Rprintf("malloc is not working !!");return -1;}*/
      crash1->edgeId = e_id;
      dist_upto_crash = (cr->tp)*edgeLength;
      crash1->tp = ((dist_upto_crash)/dist);
      crash1->frequency = cr->frequency;
#ifdef INHOM
      crash1->lambda = cr->lambda;
#endif
      if(FN(ord_list_ins_next)(cr_list,(void *)(crash1)) != 0) return -1;
    }
  }
  return 0;
}

/********************* k_function_v1 ************************************/

int FN(k_function_v1)(Graph *graph,
		  PathVertex *start,
		  double *distance,
		  double *tme_up_lt_vec,
                  int *m_val_vec,
		  int *vec_len,
		  double *inv_mvals,
		  int *source_cr_freq,
#ifdef INHOM
		  double *source_cr_lambda,
#endif
                  double *MAX_Distance,
		  double *MIN_Distance,
		  int *no_of_distance){
  /** This function calculates the number of crashes within "r" distance
      from starting node **/

  ListElmt *element,
           *member;

  AdjList  *adjlist;

  PathVertex *pth_vertex,
             *adj_vertex;

  double     rem_dist;

  int        retval;

  for(element=list_head(&graph->adjlists);
      element!=NULL;
      element=list_next(element)){

    adjlist = (AdjList *)list_data(element);
    pth_vertex = (PathVertex *)(adjlist->vertex);

    if(graph->match(start,pth_vertex))
      break;
  }

  if(element == NULL){
    Rprintf("Graph is empty!\n");
    return -1;
  }

  /** Base case in recursion **/
  member = list_head(&adjlist->adjacent);
  if (member == NULL) return 0;

  double Distance = *distance;

  while(member != NULL){

    adj_vertex = (PathVertex *)list_data(member);
    /** Important step **/

    if(adj_vertex->weight > Distance){
      /** Add up (1/mvalue) for all crash points in this edge **/

      double dist_ratio = ((Distance)/(adj_vertex->weight));
      retval = FN(sum_of_inv_mvals_on_part_of_edge_v1)(adj_vertex,
						   &pth_vertex->d,
						   &dist_ratio,
						   tme_up_lt_vec,
						   m_val_vec,
						   vec_len,
						   inv_mvals,
						   source_cr_freq,
#ifdef INHOM
						   source_cr_lambda,
#endif
						   MAX_Distance,
						   MIN_Distance,
						   no_of_distance);
      if(retval != 0){
            Rprintf("sum_of_inv_mvals_on_part_of_edge_v1 did not work inside k_function!\n");
            return -1;
      }
    } else {
      rem_dist = (Distance) - (adj_vertex->weight);

      retval = FN(sum_of_inv_mvals_for_all_pts_on_edge_v1)(adj_vertex,
						       &pth_vertex->d,
						       tme_up_lt_vec,
						       m_val_vec,
						       vec_len,
						       inv_mvals,
						       source_cr_freq,
#ifdef INHOM
						       source_cr_lambda,
#endif
						       MAX_Distance,
						       MIN_Distance,
						       no_of_distance);
      if(retval != 0){
	Rprintf("sum_of_inv_mvals_on_part_of_edge_v1 did not work inside k_function!\n");
	return -1;
      }

      /** Depth-first search implementation **/
      if(FN(k_function_v1)(graph,
		       adj_vertex,
		       &rem_dist,
		       tme_up_lt_vec,
		       m_val_vec,
		       vec_len,
		       inv_mvals,
		       source_cr_freq,
#ifdef INHOM
		       source_cr_lambda,
#endif
		       MAX_Distance,
		       MIN_Distance,
		       no_of_distance) != 0) return -1;
    }

    member = list_next(member);
}
return 0;
}


/**************** sum_of_inv_mvals_on_part_of_edge_v1 ******************/

int FN(sum_of_inv_mvals_on_part_of_edge_v1)(PathVertex *adj_vert,
					double *pth_vert_d,
					double *dist_ratio,
					double *tme_up_lt_vec,
					int *m_val_vec,
					int *vec_len,
                                        double *inv_mvals,
					int *source_cr_freq,
#ifdef INHOM
					double *source_cr_lambda,
#endif
					double *MAX_Distance,
                                        double *MIN_Distance,
					int *no_of_distance){

  if(list_size(&adj_vert->crashList) == 0) return 0;

  ListElmt *cr_element = list_head(&(adj_vert->crashList));

  int vec_index = 0,
      vec_id,
      len = *vec_len;
  /** total number of time points considered for the K-function computation **/

  double  pth_d = *pth_vert_d,
         adj_wt = adj_vert->weight;
#ifdef INHOM
  double lambda1 = *source_cr_lambda;
#endif

  int freq1 = *source_cr_freq;
  double Dratio = *dist_ratio;
  while(cr_element != NULL){
    Crash       *crash = (Crash *)(list_data(cr_element));
    int          freq = crash->frequency;
#ifdef INHOM
    double       lambda2 = crash->lambda;
#endif
/** Sh-Pth-Dist from "start" node to the crash point **/
    if(crash->tp <= Dratio){

      double      dist_frm_start = (pth_d) + ((adj_wt)*(crash->tp));

      for(vec_id = vec_index; vec_id < len ; ++vec_id){

	if (dist_frm_start <= tme_up_lt_vec[vec_index]){

	  int mv = m_val_vec[vec_index];
#ifndef INHOM
	  double inv_m_val= freq * freq1/(((double)mv));
#else
	  double inv_m_val= freq * freq1/(((double)mv) * lambda1 * lambda2);
#endif
	  OK(allot_inv_mvals_in_dist_array)(&dist_frm_start,
					&inv_m_val,
					inv_mvals,
					MAX_Distance,
					MIN_Distance,
					no_of_distance);
	  break;
	}

	if (((vec_id +1) < len) &&
	    (tme_up_lt_vec[vec_id] < dist_frm_start) &&
	    (dist_frm_start <= tme_up_lt_vec[vec_id+1])){

	  int mv = m_val_vec[vec_id+1];
#ifndef INHOM
	  double inv_m_val= freq * freq1/(((double)mv));
#else
	  double inv_m_val= freq * freq1/(((double)mv) * lambda1 * lambda2);
#endif
	  OK(allot_inv_mvals_in_dist_array)(&dist_frm_start,
					&inv_m_val,
					inv_mvals,
					MAX_Distance,
					MIN_Distance,
					no_of_distance);
	  break;
	}

      }
    } else {
      break;
    }

    vec_index = vec_id;
    /** Since Crash points are ordered by sh-pth-distances, we can start
	our search from the most recently searched distance **/
    cr_element = list_next(cr_element);
  }
  return 0;
}


/****************** sum_of_inv_mvals_for_all_pts_on_edge_v1 ************/

int FN(sum_of_inv_mvals_for_all_pts_on_edge_v1)(PathVertex *adj_vert,
					    double *pth_vert_d,
					    double *tme_up_lt_vec,
					    int *m_val_vec,
                                            int *vec_len,
					    double *inv_mvals,
					    int *source_cr_freq,
#ifdef INHOM
					    double *source_cr_lambda,
#endif
                                            double *MAX_Distance,
					    double *MIN_Distance,
					    int *no_of_distance){

  if(list_size(&adj_vert->crashList) == 0) return 0;

  ListElmt *cr_element = list_head(&(adj_vert->crashList));

  int vec_index = 0,
      vec_id,
      len = *vec_len;

  int freq1 = *source_cr_freq;

  double  pth_d = *pth_vert_d,
          adj_wt = adj_vert->weight;
#ifdef INHOM
  double lambda1 = *source_cr_lambda;
#endif

while(cr_element != NULL){

  Crash       *crash = list_data(cr_element);
  int         freq = crash->frequency;
#ifdef INHOM
  double      lambda2 = crash->lambda;
#endif

/** Sh-Pth-Dist from "start" node to the crash point **/
  double      dist_frm_start = (pth_d) + ((adj_wt)*(crash->tp));

  for(vec_id = vec_index; vec_id < len ; ++vec_id){

    if (dist_frm_start <= tme_up_lt_vec[vec_index]){

      double inv_m_val;
      int mv = m_val_vec[vec_index];
      if(mv > 0) {
#ifndef INHOM
	inv_m_val = freq * freq1 /( ((double) mv) );
#else
	inv_m_val = freq * freq1 /( ((double) mv) * lambda1 *lambda2);
#endif
      } else {
            Rprintf("m-value must be positive!\n");
            return -1;
      }

      OK(allot_inv_mvals_in_dist_array)(&dist_frm_start,
				    &inv_m_val,
				    inv_mvals,
				    MAX_Distance,
				    MIN_Distance,
				    no_of_distance);
      break;
    }

    if (((vec_id +1) < len) &&
	(tme_up_lt_vec[vec_id] < dist_frm_start) &&
	(dist_frm_start <= tme_up_lt_vec[vec_id+1])){

      double inv_m_val;
      int mv = m_val_vec[vec_id+1];
      if(mv > 0) {
#ifndef INHOM
	inv_m_val = freq * freq1 /( ((double) mv));
#else
	inv_m_val = freq * freq1 /( ((double) mv) * lambda1 * lambda2);
#endif
      } else {
            Rprintf("m-value must be positive!\n");
            return -1;
      }

      OK(allot_inv_mvals_in_dist_array)(&dist_frm_start,
				    &inv_m_val,
				    inv_mvals,
				    MAX_Distance,
				    MIN_Distance,
				    no_of_distance);
      break;
    }
  }

  vec_index = vec_id;
  /** Since Crash points are ordered by sh-pth-distances, we can start
      our search from the most recently searched distance **/
  cr_element = list_next(cr_element);
 }
 return 0;
}

/**************************************************************************
 *                                                                        *
 *  This is the main function which will be called within R session.      *
 *                                                                        *
 **************************************************************************/
void FN(createGraphNet)(int *no_of_vertices,
		    int *no_of_edges,
		    int *no_of_crashes,
		    int *crash_seg,
		    double *crash_tp,
		    int *crash_freq,
#ifdef INHOM
		    double *crash_lambda,
#endif
		    int *vert_id1,
		    int *vert_id2,
		    double *edge_length,
		    double *MAX_Distance,
		    double *MIN_Distance,
		    int *no_of_distance,
		    double *max_r,
		    int *verboseIter,
		    double *K_r){

  /************ Compute |L|: Total length of the network ********************/
int Nv = *no_of_vertices;
int Ne = *no_of_edges;
int Ncr = *no_of_crashes;
int Ndis = *no_of_distance;
double Max_R = *max_r;
double t = *MAX_Distance;
#ifndef INHOM
  double  L=0;
  int edge_count;
  for(edge_count = 0; edge_count < Ne; ++edge_count)
    L += edge_length[edge_count];
#endif
 int vIter = *verboseIter;
  /********** Create the distance vector. The K-function will be computed
	      for distances in this vector **/

  double   *dist_vec;
  dist_vec = (double *) Calloc(Ndis,double);
  /*if (dist_vec == NULL) {Rprintf("Malloc is not working!\n");return;}*/
  double dist_interval = ((*MAX_Distance) - (*MIN_Distance))/(Ndis - 1);
  int dist_id;
  for(dist_id = 0; dist_id < Ndis; ++dist_id){
    dist_vec[dist_id] = *MIN_Distance + (dist_interval*dist_id);
  }

  double *store_inv_mvals;
  store_inv_mvals = (double *)Calloc(Ndis,double);
  /*if (store_inv_mvals == NULL) {Rprintf("Malloc is not working!\n");return;}*/

  /**************************************************************************/

  int (*match_ptr)(const void*,const void*);
  match_ptr = match_graph;  /** match_graph() **/

  void (*destroy_ptr)(void *data);
  destroy_ptr = pathVertex_destroy;  	/** Graph_destroy() **/

/***********************************************************************
  Step-0: Store all crash points in an array so that it can be accessed
  for creating extended shortest path tree for every crash point.
*************************************************************************/
  Crash	*crash_points;
#ifdef INHOM
  double sum_of_inv_lambdas = 0;
#else
  int no_of_total_crashes = 0;
#endif

  /** Allocate space in the heap section of the memory
      for storing crash events **/
  crash_points = (Crash *)Calloc(Ncr,Crash);
  /*if(crash_points == NULL) {Rprintf("Malloc is not working!\n");return;}*/
  /** Store "edgeId" and "tp" values for every crash in "crash_points" **/
  int crash_count;
  for(crash_count = 0; crash_count < Ncr; ++crash_count){
    crash_points[crash_count].tp = crash_tp[crash_count];
    crash_points[crash_count].edgeId = crash_seg[crash_count];
    crash_points[crash_count].frequency = crash_freq[crash_count];
#ifdef INHOM
    crash_points[crash_count].lambda = crash_lambda[crash_count];
    sum_of_inv_lambdas +=
      crash_freq[crash_count]/crash_lambda[crash_count];
#else
    no_of_total_crashes += crash_freq[crash_count];
#endif
  }

  /**********************************************************************
   *  STEP-1: BUILD THE GRAPH USING EDGE LIST AND CRASH LIST            *
   **********************************************************************/
  Graph *gBig=(Graph*)Calloc(1,Graph);
  /*if(gBig==NULL) {Rprintf("Malloc is not working!\n");return;}*/

  OK(graph_init)(gBig,match_ptr,destroy_ptr);

  if(FN(graph_building_with_crash)(gBig,
			       no_of_vertices,
			       no_of_edges,
			       no_of_crashes,
			       crash_seg,
			       crash_tp,
			       crash_freq,
#ifdef INHOM
			       crash_lambda,
#endif
			       vert_id1,
			       vert_id2,
			       edge_length) != 0)
    {Rprintf("graph_building_with_crash_v2 Failed!\n");return;}


  //Rprintf("Step-1: Building gBig is finished!\n");

  int iteration = 0;

  while(iteration < Ncr){

    int retval;

    /**********************************************************************
     *  STEP-2: INSERT CRASH POINT AS VERTEX IN THE GRAPH                 *
     **********************************************************************/

    Crash *cr = (Crash *)Calloc(1,Crash);
    /*if(cr==NULL) {Rprintf("Malloc is not working!\n");return;}*/

    cr->tp = crash_points[iteration].tp;
    cr->edgeId = crash_points[iteration].edgeId;
    cr->frequency = crash_points[iteration].frequency;
#ifdef INHOM
    cr->lambda = crash_points[iteration].lambda;
#endif

    /**
	"crash_point_in_graph_as_vertex"
	insert the crash event in the graph as a new vertex
	with vertex id: (graph->vcount + 1).

	The function includes two edges from the crash point
	to the closest vertices in the graph. It also deletes the edge
	that hold the crash event.

	Hence, in order to restore gBig to its initial form,
	we need to delete the newly formed edges and insert the edges
	that are deleted.
    **/

    PathVertex *Vertex1 = (PathVertex *)Calloc(1,PathVertex);
    /*if(Vertex1==NULL) {Rprintf("Malloc is not working!\n");return;}*/
    PathVertex *Vertex2 = (PathVertex *)Calloc(1,PathVertex);
    /*if(Vertex2==NULL) {Rprintf("Malloc is not working!\n");return;}*/

    /** Vertex1 and Vertex2 are references ot the edges
	those are deleted from gBig **/
    retval =  FN(crash_point_in_graph_as_vertex)(gBig,
					     cr,
					     Vertex1,
					     Vertex2,
					     no_of_vertices);
         /** No memory leak in this step **/

    if(retval != 0 ){
      Rprintf("crash_point_in_graph_as_vertex Failed!\n");
      return;
    }

    Free(cr);

    //Rprintf("Step-2: crash_point_in_graph_as_vertex is successful!\n");
    /*******************************************************************
     *  STEP-3: CREATE a LOCALIZED NETWORK ROOTED at the CRASH EVENT.  *
     *******************************************************************/

    /**
	"g" shall hold the localized network rooted at the crash event.
	"g" will be passed inside "restricted_shortest_v2" as an argument,
	and once the function run is completed, "g" should consist of
	all edges emanating from the reachable nodes.
    **/

    Graph *g=(Graph*)Calloc(1,Graph);
    /*if(g==NULL) {Rprintf("Malloc is not working!\n");return;}*/

    OK(graph_init)(g,match_ptr,destroy_ptr);

    PathVertex *start=(PathVertex *)Calloc(1,PathVertex);
    /*if(start==NULL) {Rprintf("Malloc is not working!\n");return;}*/

    start->data = Calloc(1,int);
    *((int *)start->data) = (Nv + 1);

    OK(list_init)(&start->crashList,NULL);

    List *Path = (List *)Calloc(1,List);
    /*if (Path==NULL) {Rprintf("Malloc is not working!\n");return;}*/

    retval = FN(restricted_shortest_v2)(gBig, start, Path, g, Max_R);
    if (retval != 0) {
      Rprintf("Computing shortest path is not successful\n!");
      return;
    }


    //Rprintf("Step-3: Computing shortest path is successful!\n");

    /***********************************************************************
     * STEP-4: BUILD GRAPH CONSISTING of the EDGES in SHORTEST-PATH TREE. *
     ***********************************************************************/
    Graph *gT=(Graph*)Calloc(1,Graph);
    /*if(gT==NULL) {
      Rprintf("Malloc is not working!\n");
      return;
    }*/
    OK(graph_init)(gT,match_ptr,destroy_ptr);

    retval = FN(spTree_restricted_v1)(gBig,gT,Path,start);
    if(retval!=0){
      Rprintf("constructing shortest path tree is not successful!\n");
      return;
    } /** No memory leak **/

    //Rprintf("Step-4: Shortest path tree construction is successful!\n");

    /***********************************************************************
     *  STEP-5: SELECT ALL EDGES THAT ARE NOT INCLUDED IN THE
     * SHORTEST-PATH TREE, "g" CONTAINS THAT INFORMATION.
     ******************************************************************/
    retval = OK(graphEdgeDifference_v2)(gT,g);
    if(retval!=0){
      Rprintf("computing graph edge difference is not successful!\n");
      return;
    }

    //Rprintf("Step-5: graphEdgeDifference_v2 is successful!\n");

    /**********************************************************************
     * STEP-6: DELETE THE SECOND RECORD OF THE EDGES THOSE ARE RECORDED TWICE.*
     *************************************************************************/
    retval = OK(deleteSameEdge)(g);
    if(retval!=0){
      Rprintf("deleting same edge is not successful!\n");
      return;
    }

    //Rprintf("Step-6: deleteSameEdge is successful!\n");
    /*******************************************************************
     *  STEP-7: BUILD THE EXTENDED SHORTEST PATH TREE FROM THE "START" NODE. *
     ************************************************************************/
    List *newPath = (List *)Calloc(1,List);
    /*if (newPath==NULL) {
      Rprintf("Malloc is not working!\n");
      return;
    }*/

    retval = FN(extended_sh_path_tree_restricted)(g,
					      gT,
					      Path,
					      newPath,
					      Max_R,
					      no_of_vertices,
					      no_of_edges);
    if(retval!=0){
      Rprintf("building extended shortest path tree is not successful with return value:%d!\n",retval);
      return;
    }

    //Rprintf("Step-7:  extended_sh_path_tree_restricted is successful!\n");

    /**
	"newPath" is the list that contains pointers to the nodes
	which are reachable from the "start" node. The PathVertex element
	corresponding to each node in the list contains the shortest-path
	distance from "start".
    **/

    /*********************************************************************
     * STEP-8: LIST the DISTANCES FROM START NODE
     * AND SORT THE NODES IN NON-DECREASING ORDER OF d-VALUES.
     ********************************************************************/
    List *sortedList;
    sortedList = (List *)Calloc(1,List);
    /*if(sortedList == NULL) {
      Rprintf("Malloc is not working!\n");
      return;
    }*/
    OK(list_init)(sortedList,NULL);
    retval = OK(sorted_list)(newPath,sortedList);
    if(retval!=0){
      Rprintf("sorted_list is not successful!\n");
      return;
    }

    //Rprintf("Step-8:  sorted_list is successful!\n");
    /*******************************************************************
     * STEP-9: INCLUDE THE DEGREE INFORMATION FOR EVERY NODE IN PATHVERTEX. *
     ************************************************************************/
    OK(include_vertex_degree_info_inside_pthvert)(gT, sortedList);

    //Rprintf("Step-9:  include_vertex_degree_info_inside_pthvert is successful!\n");
    /********************************************************************
     *  STEP-10: Compute m-value for computing k-function
     **********************************************************************/
    int sortedList_size = list_size(sortedList);

    double *time_array = (double *)Calloc(sortedList_size,double);
    /*if(time_array == NULL) {Rprintf("Malloc is not working!\n");return;}*/

    int *degree_array = (int *)Calloc(sortedList_size,int);
    /*if(degree_array == NULL) {Rprintf("Malloc is not working!\n");return;}*/

    /** time_array is introduced to store the shortest path distances
	in a non-increasing order from "start" node to other vertices
	in the extended shortest path tree  **/
    /** degree_array stores the degree informations
	of the corresponding nodes **/
    OK(create_distance_and_degree_array)(sortedList,time_array,degree_array);

    double *tme_uppr_lmt = (double *)Calloc(sortedList_size,double);
    /*if(tme_uppr_lmt == NULL) {Rprintf("Malloc is not working!\n");return;}*/

    int *m_val = (int *)Calloc(sortedList_size,int);
    /*if(m_val == NULL) {Rprintf("Malloc is not working!\n");return;}*/

    int array_size=0;

    OK(create_m_val_array)(time_array,
		       degree_array,
		       tme_uppr_lmt,
		       m_val,
		       sortedList_size,
		       &array_size);

    tme_uppr_lmt = (double *)Realloc(tme_uppr_lmt,array_size,double);
    m_val = (int *)Realloc(m_val,array_size,int);

    //Rprintf("Step-10:  create_m_val_array is successful!\n");
    /******************************************************************
     *  STEP-11: Compute K-function for a distance vector
     ****************************************************************/
    /*double t = *MAX_Distance;*/ /* t is defined outside the while loop */
    int cr_freq_wt = crash_points[iteration].frequency;
#ifdef INHOM
    double cr_lambda = crash_points[iteration].lambda;
#endif

    retval = FN(k_function_v1)(gT,
			   start,
			   &t,
			   tme_uppr_lmt,
			   m_val,
			   &sortedList_size,
			   store_inv_mvals,
			   &cr_freq_wt,
#ifdef INHOM
			   &cr_lambda,
#endif
			   MAX_Distance,
			   MIN_Distance,
			   no_of_distance);

    if (retval != 0){Rprintf("k_function did not work.\n"); return;}

    //Rprintf("Step-11:  k_function_v1 is successful!\n");
    /**************************************************************
     *               FREE EVERYTHING
     **************************************************************/
    OK(delete_vertex_from_graph)(gBig,start);
    /** Delete the crash vertex from gBig **/

    //Rprintf("Step-12:  delete START vertex from graph is successful!\n");
    if(OK(graph_ins_edge)(gBig,Vertex1,Vertex2)!=0){
      Rprintf("graph_ins_edge did not work!!\n");
      return;
    }

    if(OK(graph_ins_edge)(gBig,Vertex2,Vertex1)!=0){
      Rprintf("graph_ins_edge did not work!!\n");
      return;
    }

    OK(pathVertex_destroy)(start);

    //Rprintf("Step-13:  PathVertex START is freed successfully!\n");

    OK(graph_destroy)(g);
    Free(g);

    //Rprintf("Step-14:  The Graph g is freed successfully!\n");

    OK(path_destroy)(Path);
    Free(Path);

    OK(path_destroy)(newPath);
    Free(newPath);

    OK(path_destroy)(sortedList);
    Free(sortedList);

    OK(graph_destroy)(gT);
    Free(gT);

    //Rprintf("Step-15:  The Graph gT is freed successfully!\n");

    Free(time_array);Free(degree_array);Free(tme_uppr_lmt);Free(m_val);
    if(vIter == 1){
      int iter1;
      iter1 = (iteration + 1);
      Rprintf("iteration no: %d\n",iter1);
    }
    ++iteration;
  }

  /*** Store the K(r) values for all r in an array **/
  double sum_of_kvals_for_all_points = 0,
         K_Function_Value,
         denom;

#ifdef INHOM
  denom = sum_of_inv_lambdas;
#else
  denom = no_of_total_crashes * (no_of_total_crashes - 1)/L;
#endif

  for (dist_id = 0; dist_id < Ndis; ++dist_id){
    sum_of_kvals_for_all_points += store_inv_mvals[dist_id];
    K_Function_Value = sum_of_kvals_for_all_points/denom;
    K_r[dist_id] = K_Function_Value;
  }

  Free(crash_points);
  Free(store_inv_mvals);
  Free(dist_vec);
  OK(graph_destroy)(gBig); Free(gBig);
  return;
}


/*****************************  crash_compare  ***************************/

int FN(crash_compare)(Crash *crash1,Crash *crash2){
if((crash2->tp)>(crash1->tp))
    return 1;
else
    return 0;
}


/*****************************  ord_list_ins_next  ***************************/

int FN(ord_list_ins_next)(List *ord_list,void *data){
ListElmt *element,*prev_element;
Crash *crash;
int retval;
if(list_size(ord_list)==0){
  retval=OK(list_ins_next)(ord_list,NULL,data);
  if(retval!=0){
    Rprintf("list_ins_next did not work inside ord_list_ins_next.\n");
    return -1;
  }
}
else{
    prev_element=NULL;
    element=list_head(ord_list);
    crash = (Crash *)(element->data);

    while(element!=NULL){
        crash = (Crash *)element->data;
        if(FN(crash_compare)((Crash *)data,crash)==1){
            break;
        }
        else{
            prev_element=element;
            element=element->next;
        }
    }
    retval = OK(list_ins_next)(ord_list,prev_element,data);
    if(retval !=0){
      Rprintf("list_ins_next did not work inside ord_list_ins_next.\n");
      return -1;
    }
 }
 return 0;
}


/************* delete_single_crash_from_crashlist ***********************/

int FN(delete_single_crash_from_crashlist)(List *crlist, Crash *cr){

  ListElmt   *element, *prev_element;

  Crash      *crash, *data2;

  int         retval;

  if(list_size(crlist) == 0){
    Rprintf("Crash List can not be empty!\n");
    return -1;
  }

  prev_element = NULL;

  for(element = list_head(crlist);
      element != NULL;
      element = list_next(element)){

    crash = (Crash *)(element->data);

    if((crash->tp == cr->tp) && (crash->edgeId == cr->edgeId)){

      break;
    }
    prev_element = element;
  }
  if (element == NULL) {
    Rprintf("No crash match has been observed!\n");
    return -1;
  }

  retval = OK(list_rem_next)(crlist,prev_element,(void **)&data2);

  if(retval == 0){Free(data2);} else{
    Rprintf("list_rem_next did not work!\n"); return -1;
  }
  return 0;
}

/******************** store_edge_before_deleting ************************************/

void FN(store_edge_before_deleting)(PathVertex *pth_vertex,
				    PathVertex *adj_vertex,
				    PathVertex *adj_vert1,
				    PathVertex *adj_vert2){

 /** pth_vertex---->adj_vertex is the edge that hold the crash event under consideration. **/
 /** Before deleting the edges pth_vertex---->adj_vertex and adj_vertex---->pth_vertex from
     the "graph" inside "crash_point_in_graph_as_vertex_v1", we save these edge informations
     inside "adj_vert1" and "adj_vert2".
     The space for "adj_vert1" and "adj_vert2" had already been allocated in the main function. **/

    /** pth_vertex information will be stored in adj_vert1 and adj_vertex information will be stored in adj_vert2 **/

    adj_vert2->data = Calloc(1,int); /*if(adj_vert2->data==NULL){Rprintf("malloc did not work!!\n");return;}*/
    *((int *)(adj_vert2->data)) = *((int *)(adj_vertex->data));
    adj_vert2->edgeID = adj_vertex->edgeID;
    adj_vert2->weight = adj_vertex->weight;
    OK(list_init)(&(adj_vert2->crashList),&destroy_crash);
    if(FN(copy_crash_list_v2)(&(adj_vertex->crashList),
			      &(adj_vert2->crashList)) != 0){
      Rprintf("copy_crash_list did not work!!\n");return;
    }

    /** Insert the pth_vertex information in adj_vert1 **/
    adj_vert1->data = Calloc(1,int); /*if(adj_vert1->data==NULL){Rprintf("malloc did not work!!\n");return;}*/
    *((int *)(adj_vert1->data)) = *((int *)(pth_vertex->data));
    adj_vert1->edgeID = adj_vertex->edgeID;
    adj_vert1->weight = adj_vertex->weight;
    OK(list_init)(&(adj_vert1->crashList),&destroy_crash);
    if(FN(copy_crash_list_rev_v2)(&(adj_vertex->crashList),
				  &(adj_vert1->crashList)) != 0){
      Rprintf("copy_crash_list_rev did not work!!\n");return;
    }
    return;
}

/***************************** crash_point_in_graph_as_vertex ************************************/

int FN(crash_point_in_graph_as_vertex)(Graph *graph,Crash *crashPt,PathVertex *adj_vert1,
                                      PathVertex *adj_vert2,int *no_of_vertices){

/*************************************************************************************************
    This function includes a single crash point in the graph as a new vertex.                    *
    We first identify the edge that contains the crash point in its crash list, e.g. suppose     *
    v1------(crash)-->v2. Then we include two new edges in the graph, (crash)------>v1 and       *
    (crash)-->v2 destroying two existing edges v1--->v2 and v2--->v1. In the process we split    *
    the crash list that exists between v1 and v2 and distribute the crashes to the newly created *
    edges.This function will be called before the computation of the extended shortest path tree *
    from the given crash point.                                                                  *
*************************************************************************************************/

          /*********************************************************************
          * The memory allocated to the edges that are removed in this function*
          * are returned back to the heap at the end of the program using the  *
          * "pathVertex_destroy" function.                                     *
          *********************************************************************/
ListElmt     *element, *member;

AdjList      *adjlist;

PathVertex   *pth_vertex,*adj_vertex,
             *pth_node,*adj_node,
             *crashNode;

double       dist1,dist2;

int          ret_val;

int          element_key = 1, edgeCount = (graph->ecount)/2;
/** Let us find the edge that has the crash in its crash list **/

    element = list_head(&graph->adjlists);

    if(element==NULL){Rprintf("The Graph is empty!\n");return -1;}

    while(element!=NULL){

        adjlist = (AdjList *)list_data(element);
        pth_vertex = (PathVertex *)(adjlist->vertex);

        /** Check if there is any outgoing edge from the pth_vertex. **/
        member = list_head(&adjlist->adjacent);
        if(member == NULL){
        /** If There is no outgoing edge, move to the next pth_vertex in the graph. **/
                element = list_next(element);
        }
        else{
        /** Visit all outgoing edges listed in the adjacency-list of the pth_vertex until
            find the match with the crash edgeid or reach the end of the adjacency-list. **/
                while(member != NULL){

                    adj_vertex = list_data(member);
                /** Check if crash has occurred on this edge **/
                    if(adj_vertex->edgeID == crashPt->edgeId){
                        /** We found the match.So get out of this "while" loop. **/
                        element_key = 0;
                        break;          /**** BREAK BREAK BREAK IS HERE ****/
                    }
                    else{
                        /** Crash is not on this edge, so check the next edge on the adjacency list. **/
                        member = list_next(member);
                    }

                }
        /** If crash was found on any edge emanating from pth_vertex, element_key=0. Then we need to set
            "element=NULL". Otherwise, we set "element=list_next(element)". **/
                if(element_key == 0){
                    element = NULL;
                }
                else{
                    element = list_next(element);
                    element_key = 1;
                }

        }

    }

    /** Return if there is no match **/
    if(member == NULL){Rprintf("There is no match found between crash segment id and edge id in the graph!\n"); return -1;}

/** If the above piece of code runs properly, then crashPt belongs to the edge pth_vertex--->adj_vertex **/

    /** Before proceeding towards entering new edges and deleting the old ones, we shall save the information of the edges which are going to be deleted **/
    FN(store_edge_before_deleting)(pth_vertex,adj_vertex,adj_vert1,adj_vert2);

    int Nv = *no_of_vertices;
    /** Start inserting crash as a node in the graph **/
    crashNode = (PathVertex *)Calloc(1,PathVertex);/*if(crashNode==NULL){Rprintf("malloc is not working !!");return -1;}*/
    crashNode->data = Calloc(1,int);/* if(crashNode->data==NULL){Rprintf("malloc is not working !!");return -1;}*/
    *((int *)(crashNode->data)) = (Nv + 1);
    OK(list_init)(&(crashNode->crashList),NULL);
    ret_val = OK(graph_ins_vertex)(graph,crashNode); if(ret_val!= 0){Rprintf("graph_ins_vertex did not work!\n"); return -1;}

    /** Prepare path vertices that will enter in the adjacency list of the crashNode **/
        /** Firstly, (crashPt)------>pth_node   **/
    dist1 = (crashPt->tp)*(adj_vertex->weight);
    pth_node = (PathVertex *)Calloc(1,PathVertex); /*if(pth_node==NULL){Rprintf("malloc is not working !!");return -1;}*/
    pth_node->data = Calloc(1,int); /* if(pth_node->data==NULL){Rprintf("malloc is not working !!");return -1;} */
    *((int *)pth_node->data) = *((int *)pth_vertex->data);
    pth_node->weight = dist1;
    pth_node->edgeID = (edgeCount +1);
    OK(list_init)(&(pth_node->crashList),&destroy_crash);
        /** Secondly, (crashPt)------>adj_node   **/
    dist2 = (adj_vertex->weight - dist1);
    adj_node = (PathVertex *)Calloc(1,PathVertex); /* if(adj_node==NULL){Rprintf("malloc is not working !!");return -1;} */
    adj_node->data = Calloc(1,int); /* if(adj_node->data==NULL){Rprintf("malloc is not working !!");return -1;} */
    *((int *)adj_node->data) = *((int *)adj_vertex->data);
    adj_node->weight = dist2;
    adj_node->edgeID = (edgeCount +2);
    OK(list_init)(&(adj_node->crashList),&destroy_crash);


    /** Delete the crash point from the crashList of adj_vertex **/
    ret_val = FN(delete_single_crash_from_crashlist)(&adj_vertex->crashList,crashPt);
    if(ret_val != 0){Rprintf("Deleting the crash from the crashList of adj_vertex has not worked!\n");return -1;}

    /** Now break the crashList in adj_vertex into two separate lists for newly created edges **/
    ret_val =
      FN(break_crash_list_into_two_lists_rev_v2)(&(pth_node->crashList),
						 &(adj_node->crashList),
						 &adj_vertex->crashList,
			 ((crashPt->tp)*(adj_vertex->weight)),
						 (adj_vertex->weight),
						 (edgeCount +1),
						 (edgeCount +2));

    /** Now insert the edges to the graph **/
    ret_val = OK(graph_ins_edge)(graph,crashNode,pth_node);
    if(ret_val!=0){Rprintf("Edge insertion did not work!\n");return -1;}
    ret_val = OK(graph_ins_edge)(graph,crashNode,adj_node);
    if(ret_val!=0){Rprintf("Edge insertion did not work!\n");return -1;}

    /** Delete the edge pth_vertex--->adj_vertex from the graph **/
    /** Before deleting the egde from "graph", store the adj_vertex information in "adj_vert1" and "adj_vert2" **/

    void *data1,*vd_ptr, **data2,
         *data1_x,*vd_ptr_x, **data2_x;

    data1 = pth_vertex;
    vd_ptr = adj_vertex;
    data2 = &vd_ptr;
    ret_val = OK(graph_rem_edge)(graph,data1,data2);
    if(ret_val!=0){Rprintf("Edge removal did not work!\n");return -1;}

    /** Now, data2 shall hold the address of the pointer to the PathVertex structure that has just been removed! **/
    /** We need to destroy data2 **/
    /** Delete the edge adj_vertex--->pth_vertex from the graph **/

    data1_x = adj_vertex; /** Problem is here! Once adj_vertex is destroyed we should not try to access it **/
    vd_ptr_x = pth_vertex;
    data2_x = &vd_ptr_x;
    ret_val = OK(graph_rem_edge)(graph,data1_x,data2_x);
    if(ret_val!=0){Rprintf("Edge removal did not work!\n");return -1;}

    OK(pathVertex_destroy)(*data2);
    OK(pathVertex_destroy)(*data2_x);

    return 0;

}

/*********** extended_sh_path_tree_restricted ****************************/

int FN(extended_sh_path_tree_restricted)(Graph *gLink,
				     Graph *gTree,
				     List *Vlist,
                                     List *newPath,
				     double R,
				     int *no_of_vertices,
                                     int *no_of_edges)
{
  ListElmt *element,*member,*el;

  AdjList  *adjlist;

  PathVertex *pth_vertex,*adj_vertex,
             *brk_pt1,*brk_pt2,*brk_pt,
             *path_vert,*parent_vert;

  double     dist1,dist2,link_dist1,link_dist2,link_dist,
             total_link_dist,dist;

  int         ret_val,retval,
              eCount = (*no_of_edges + 2),
  /** eCount is the total no. of edges we started with in the original graph **/
              vCount = (*no_of_vertices + 1),
              case_indicator=0,
              l=0,m=0;

  /**
      gLink contains all edges that are not present in the
      shortest path tree gTree.

      gLink also contains all shortest path distances from the start node
      to all other vertices in the graph

      We shall go through every edge in gLink using a for loop
      in order to include the breaking points and
      newly constructed linking edges to gTree
  **/
  for(element = list_head(&gLink->adjlists);
      element != NULL;
      element = list_next(element)){
    adjlist = (AdjList *)list_data(element);
    pth_vertex = adjlist->vertex;
    member = list_head(&adjlist->adjacent);
    /**
	If there is no edge present in the adjacency list of pth_vertex,
	move to the next vertex in the graph
    **/
    if(member==NULL){continue;} else {
      dist1 = pth_vertex->d;
      /** dist1 is the shortest path distance from the start node
	  to pth_vertex **/
      /** Visit every vertices in the adjaceny list of the pth_vertex
	  to include the breaking points and linking edges **/
      while(member != NULL){
	adj_vertex = list_data(member);
	/** Get "adj_vertex" vertex of the edge pth_vertex----->adj_vertex **/
	for(el=list_head(Vlist);el!=NULL;el=list_next(el)){
	  parent_vert = (PathVertex *)list_data(el);
	  if(*((int *)parent_vert->data) ==  *((int *)pth_vertex->data)){
	    break;
	  }
	}
	if (el == NULL){
	  Rprintf("The pth_vertex does not exist in Vlist, i.e. pth_vertex is not reachable!!\n");
	  return -1;
	}

	for(el=list_head(Vlist);el!=NULL;el=list_next(el)){
	  path_vert = (PathVertex *)list_data(el);
	  if(*((int *)path_vert->data) ==  *((int *)adj_vertex->data)){
	    case_indicator = 1;
	    dist2 = path_vert->d;
	    break;
	  }
	}
	if(el == NULL){case_indicator = 2;}
	/** el==NULL implies that there's no match between the nodes
	    in the Vlist and adj_vertex.**/

	/**************************************************************
	 * (case_indicator = 1) => The edge pth_vertex---->adj_vertex
	 * needs to be broken into two linking edges:
	 * pth_vertex--->brk_pt1 and adj_vertex--->brk_pt2.
	 * Both breakpoints need to be inserted in gTree.

	 * (case_indicator = 2) => The node adj_vertex is not present in gTree.
	 * Therefore, first the node needs be included.
	 * Then only the edge pth_vertex--->adj_vertex will be
	 * included in gTree.
	 ************************************************************/

	if(case_indicator == 1){
	  total_link_dist = adj_vertex->weight;
	  /** Length of the edge pth_vertex-------->adj_vertex **/
	  dist = ((dist2-dist1)+total_link_dist)/2.0;
     /** pth_vert(d-value) + dist = adj_vert(d-value) + (edge_length - dist) **/
	  link_dist1 = dist; link_dist2 = (total_link_dist - dist);
	  /** Include the breaking points as new vertices
	      into gTree graph structure **/
	  if(link_dist1 != 0){
	    PathVertex *node1 = (PathVertex *)Calloc(1,PathVertex);
	    /*if(node1==NULL){Rprintf("malloc is not working !!");return -1;}*/
	    /**
		Create the first break point and add the edge
		from pth_vertex to brk_pt1 (pth_vertex--->brk_pt1)
	    **/
	    brk_pt1 = (PathVertex *)Calloc(1,PathVertex);
	    /*if(brk_pt1 == NULL){Rprintf("malloc is not working !!");return -1;}*/
	    ++l; ++m;
	    /** graph_ins_breakpt modifies gTree by entering a new vertex
		"node1". It adds the edgeID and weight information
		in "brk_pt1" and also initializes its crashList.
		Consequently, it prepares "brk_pt1" for the adjacency list
		of path_vertex **/
	    node1->d = pth_vertex->d + link_dist1;
	    node1->parent = parent_vert;

	    ret_val = OK(graph_ins_breakpt)(gTree,
					node1,
					brk_pt1,
					(vCount + l),
					(eCount + m),
					link_dist1);
	    if (ret_val != 0){
	      Rprintf("graph_ins_breakpt did not work!\n");
	      return -1;
	    }
	  }
	  if(link_dist2 != 0){
	    PathVertex *node2 = (PathVertex *)Calloc(1,PathVertex);
	    /*if(node2==NULL){Rprintf("malloc is not working !!");return -1;}*/

	    brk_pt2 = (PathVertex *)Calloc(1,PathVertex);
	    /*if(brk_pt2 == NULL){Rprintf("malloc is not working !!");return -1;}*/
	    ++l;++m;

	    node2->d = dist2 + link_dist2;
	    node2->parent = path_vert;
	    /** Using arguments similar to above we have the following **/
	    ret_val = OK(graph_ins_breakpt)(gTree,
					node2,
					brk_pt2,
					(vCount + l),
					(eCount + m),
					link_dist2);
	    if (ret_val != 0){
	      Rprintf("graph_ins_breakpt did not work!\n");
	      return -1;
	    }
	  }
	  /** Divide the crash list that exist between pth_vertex
	      and adj_vertex. Subsequently, construct two new crash lists
	      for the linking edges **/
	  if((link_dist1 != 0) && (link_dist2 != 0)){
	    ret_val =
	      FN(break_crash_list_into_two_lists_v2)(&(brk_pt1->crashList),
						 &(brk_pt2->crashList),
						 &(adj_vertex->crashList),
						 link_dist1,
						 total_link_dist,
						 brk_pt1->edgeID,
						 brk_pt2->edgeID);
	    if(ret_val == 0){
	      if(OK(graph_ins_edge)(gTree,pth_vertex,brk_pt1)!=0) return -1;
	      if(OK(graph_ins_edge)(gTree,adj_vertex,brk_pt2)!=0) return -1;

	      member = list_next(member);
	      /** Move to the next edge, if the insertion of linking edges
		  is successful! **/
	    } else {
	      Rprintf("Crash List distribution did not work!\n");
	      return -1;
	    }
	  } else {
	    /** Either dist1=0 or dist2=0. Both cannot be zero together **/

	    if(link_dist1 == 0){
	      /** There is no linking edge emanating from pth_vertex.
		  The entire crashList is going to be copied
		  (tp values in reverse order) inside adj_vertex----->brk_pt2.
	      **/
	      ret_val =
		FN(break_crash_list_into_one_list)(&(brk_pt2->crashList),
					       &(adj_vertex->crashList),
					       link_dist1,
					       brk_pt2->edgeID);
	      if(ret_val == 0){
		if(OK(graph_ins_edge)(gTree,adj_vertex,brk_pt2)!=0) return -1;
		/** Insert linking edge adj_vertex----->brk_pt2 **/
		member = list_next(member);
	      } else {
		Rprintf("Crash List distribution did not work!\n");
		return -1;
	      }
	    } else {
	      /** link_dist2 = 0. This means that there is no linking edge
		  emanating from adj_vertex. The entire crashList
		  going to be copied inside pth_vertex----->brk_pt1
	      **/
	      ret_val =
		FN(break_crash_list_into_one_list)(&(brk_pt1->crashList),
					       &(adj_vertex->crashList),
					       link_dist1,
					       brk_pt1->edgeID);
	      if(ret_val == 0){
		if(OK(graph_ins_edge)(gTree,pth_vertex,brk_pt1)!=0) return -1;
		member = list_next(member);
	      } else {
		Rprintf("Crash List distribution did not work!\n");
		return -1;
	      }
	    }
	  }
	  /** Before moving to the next edge (adj_vertex=list_data(member)),
	      we should destroy the edge in pth_vertex----->adj_vertex
	      in gLink
	  **/
	  OK(delete_edge_from_graph)(gLink,pth_vertex,adj_vertex);
	}
	/** case_indicator == 2 **/
	if (case_indicator == 2){
	  total_link_dist = adj_vertex->weight;
	  if(total_link_dist <= 0){
	    Rprintf("Invalid edge weight!\n");
	    return -1;
	  }
	  if ((dist1 + total_link_dist) < R) {
	    Rprintf("The node should have been included in the tree!!\n");
	    return -1;
	  }
	  link_dist = R - dist1;
	  /** Need to include a breaking point in the tree **/
	  PathVertex *adj_node = (PathVertex *)Calloc(1,PathVertex);
	  /*if(adj_node==NULL){Rprintf("malloc is not working !!");return -1;}*/
	  brk_pt = (PathVertex *)Calloc(1,PathVertex);
	  /*if(brk_pt==NULL){Rprintf("malloc is not working !!");return -1;}*/

	  adj_node->d = R;
	  adj_node->parent = parent_vert;

	  ++l; ++m;

	  retval = OK(graph_ins_breakpt)(gTree,
				     adj_node,
				     brk_pt,
				     (vCount+l),
				     (eCount+m),
				     link_dist);
	  if (retval != 0){
	    Rprintf("graph_ins_breakpt did not work!\n");
	    return -1;
	  }
	  retval =
	    FN(break_crash_list_before_max_dist)(&(brk_pt->crashList),
					     &(adj_vertex->crashList),
					     link_dist,
					     total_link_dist,
					     brk_pt->edgeID);
	  if (retval !=0){
	    Rprintf("break_crash_list_before_max_dist did not work!\n");
	    return -1;
	  }

	  retval = OK(graph_ins_edge)(gTree,pth_vertex,brk_pt);
	  if(retval != 0){
	    Rprintf("graph_ins_edge did not work!!\n");
	    return -1;
	  }
	  member = list_next(member);

	  OK(delete_edge_from_graph)(gLink,pth_vertex,adj_vertex);

	}

      }
    }
  }

  OK(list_init)(newPath, NULL);

  for (element = list_head(&graph_adjlists(gTree));
       element != NULL;
       element = list_next(element))
    {
      pth_vertex = ((AdjList *)list_data(element))->vertex;

      if (OK(list_ins_next)(newPath, list_tail(newPath), pth_vertex) != 0)
	{
	  OK(list_destroy)(newPath);
	  return -1;
	}
    }
  return 0;
}



/******************** restricted_shortest_v2 *****************************/

int FN(restricted_shortest_v2)(Graph *graph, const PathVertex *start, List *paths,
                           Graph *exhaustiveTree, double Rmax)
{
  AdjList                     *adjlist;

  PathVertex                  *pth_vertex,
                              *adj_vertex;

  ListElmt                    *element,
                              *member;

  double                      minimum;

  int                         found, i;

  /***********************************************************************
   *                                                                      *
   * Initialize all of the vertices in the graph.                         *
   *                                                                      *
   ***********************************************************************/

  found = 0;

  for (element = list_head(&graph_adjlists(graph));
       element != NULL;
       element = list_next(element))
    {
      pth_vertex = ((AdjList *)list_data(element))->vertex;

      if (graph->match(pth_vertex, start))
        {
	  /***************************************************************
	   * Initialize the start vertex.                                 *
	   ****************************************************************/
	  pth_vertex->color = white;
	  pth_vertex->d = 0;
	  pth_vertex->parent = NULL;
	  found = 1;
        } else {
	  /**************************************************************
	   *                                                             *
	   * Initialize vertices other than the start vertex.            *
	   *                                                             *
	   ***************************************************************/
	  pth_vertex->color = white;
	  pth_vertex->d = DBL_MAX;
	  pth_vertex->parent = NULL;
        }
    }
  /***********************************************************************
   *                                                                      *
   * Return if the start vertex was not found.                            *
   *                                                                      *
   ************************************************************************/
  if (!found)
    {
      Rprintf("The start vertex was not found!");
      return -1;
    }
    /*******************************************************************
    *
    * Use Dijkstra's algorithm to compute shortest paths
    * from the start vertex.
    *
    ********************************************************************/

  i = 0;
  int Nv = graph_vcount(graph);
  while (i < Nv)
    {
      /*********************************************************************
       * Select the white vertex with the smallest shortest-path estimate. *
       ********************************************************************/
      minimum = DBL_MAX;

      for (element = list_head(&graph_adjlists(graph));
	   element != NULL;
	   element = list_next(element))
        {
	  pth_vertex = ((AdjList *)list_data(element))->vertex;
	  if (pth_vertex->color == white && pth_vertex->d < minimum)
            {
	      minimum = pth_vertex->d;
	      adjlist = list_data(element);
            }
	}
      /******************************************************************
       *                                                                *
       * Color the selected vertex black.                               *
       *                                                                *
       *****************************************************************/
      /** Before coloring the vertex, check if its d-value estimate
	    is within Rmax **/
      if(((PathVertex *)(adjlist->vertex))->d > Rmax){break;}

      ((PathVertex *)(adjlist->vertex))->color = black;

      /*****************************************************************
       *                                                               *
       * Traverse each vertex adjacent to the selected vertex.         *
       *                                                               *
       ****************************************************************/

      for (member = list_head(&adjlist->adjacent);
	   member != NULL;
	   member = list_next(member))
        {
	  adj_vertex = list_data(member);

	  /****************************************************************
	   *                                                              *
	   * Find the adjacent vertex in the
	   * list of adjacency-list structures.
	   ****************************************************************/
	  for (element = list_head(&graph_adjlists(graph));
	       element != NULL;
	       element = list_next(element))
            {
	      pth_vertex = ((AdjList *)list_data(element))->vertex;

	      if (graph->match(pth_vertex, adj_vertex))
                {
		  /********************************************************
		   *
		   * Relax the adjacent vertex in the list of
		   * adjacency-list structures.       *
		   *********************************************************/

		  relax(adjlist->vertex, pth_vertex, adj_vertex->weight);
		  /** insert the edge in exhaustiveTree **/
		  if(FN(insert_edge_in_subnet)(exhaustiveTree,
					   ((PathVertex *)(adjlist->vertex)),
					   adj_vertex)!=0){
		    Rprintf("insert_edge_in_subnet did not work in restricted_shortest_v1!\n");
		    return -1;
		  }
                }
            }
	}
      /*****************************************************************
       *                                                               *
       * Prepare to select the next vertex.                            *
       *                                                               *
       *****************************************************************/
      i++;
    }

  /*******************************************************************
   *
   * Load the vertices with their path information into a list.
   *                                                                  *
   *******************************************************************/

  list_init(paths, NULL);
  for (element = list_head(&graph_adjlists(graph));
       element != NULL;
       element = list_next(element))
    {
      /**********************************************************************
       *                                                                    *
       * Load each black vertex from the list of adjacency-list structures. *
       **********************************************************************/

      pth_vertex = ((AdjList *)list_data(element))->vertex;
      if (pth_vertex->color == black)
        {

	  if (OK(list_ins_next)(paths, list_tail(paths), pth_vertex) != 0)
            {
	      OK(list_destroy)(paths);
	      return -1;
            }
        }
    }
  return 0;
}


/***************************** insert_edge_in_subnet ************************************/

int FN(insert_edge_in_subnet)(Graph *g,PathVertex *p,PathVertex *q){

/** This function is called inside "restricted_shortest_v1" function. We intend to build a
    localized network around the crash event by including all edges that are reachable within
    a given max distance. The graph would consist of all edges emanating from nodes that
    are reachable from the crash event subject to the max distance constraint.      **/

/** The information regarding which edge to include would be passed by the PathVertex elements p and q **/

/** Check if "p" and "q" are present in "g" as vertices **/
ListElmt        *element,*member;

AdjList         *adjlist, *adjlist1;

PathVertex      *pth_vertex,*adj_vertex;

int             indicator1 = 0,indicator2 = 0,indicator12 = 0;

for(element = list_head(&g->adjlists);element != NULL;element = list_next(element)){

    adjlist = (AdjList *)(list_data(element));
    pth_vertex = (PathVertex *)(adjlist->vertex);

    if(g->match(p,pth_vertex)){

        adjlist1 = adjlist;
        ((PathVertex *)(adjlist1->vertex))->d = p->d;
        ((PathVertex *)(adjlist1->vertex))->parent = p->parent;
        indicator1 = 11;
    }
    if(g->match(q,pth_vertex)){

        indicator2 = 11;
    }

}

/** In subnetwork "g", the appropriate d-value estimate would not be possible to store **/
if(indicator1 != 11){

    /** Insert "p" as a vertex in "g" **/
    PathVertex *pth_vert1 = (PathVertex *)Calloc(1,PathVertex); /*if(pth_vert1==NULL){Rprintf("malloc is not working !!");return -1;}*/

    pth_vert1->data = Calloc(1,int); /*if(pth_vert1->data==NULL){Rprintf("malloc is not working !!");return -1;}*/

    *((int *)(pth_vert1->data)) = *((int *)(p->data));

    OK(list_init)(&pth_vert1->crashList,&destroy_crash);

    if(OK(graph_ins_vertex)(g,pth_vert1) != 0){
      Rprintf("graph_ins_vertex did not work inside insert_edge_in_subnet!\n");
      return -1;
    }
 } else {

    /** Check if the adjacecny list of "p" contain "q". Then we do not need insert the edge **/
  for(member = list_head(&adjlist1->adjacent);
      member != NULL;
      member=list_next(member)){

    adj_vertex = (PathVertex *)(list_data(member));
        if(g->match(adj_vertex,q)){indicator12 = 11;break;}
    }
 }
/** If indicator12 == 11, then no need to do anything. Just return **/
 if(indicator12 == 11){
   return 0;
 } else {
   /** Check whether we need to insert "q" as a vertex or not **/
   if(indicator2 != 11){
     /** Insert "q" as a vertex in "g" **/
     PathVertex *pth_vert2 = (PathVertex *)Calloc(1,PathVertex);
     /*if(pth_vert2==NULL){Rprintf("malloc is not working !!");return -1;}*/
     pth_vert2->data = Calloc(1,int);
     /*if(pth_vert2->data==NULL){Rprintf("malloc is not working !!");return -1;}*/
     *((int *)(pth_vert2->data)) = *((int *)(q->data));

     list_init(&pth_vert2->crashList,&destroy_crash);
     if(OK(graph_ins_vertex)(g,pth_vert2) != 0){
       Rprintf("graph_ins_vertex did not work inside insert_edge_in_subnet!\n");
       return -1;
     }
   }
   /** Insert "q" in the adjacency list of "p" **/
   PathVertex *adj_node = (PathVertex *)Calloc(1,PathVertex);
   /*if(adj_node==NULL){Rprintf("malloc is not working !!");return -1;}*/
   adj_node->data = Calloc(1,int);
   /*if(adj_node->data==NULL){Rprintf("malloc is not working !!");return -1;}*/

   *((int *)(adj_node->data)) = *((int *)(q->data));
   adj_node->edgeID = q->edgeID;
   adj_node->weight = q->weight;
   list_init(&adj_node->crashList,&destroy_crash);
   if(FN(copy_crash_list_v2)(&q->crashList,
			     &adj_node->crashList)!=0){
     Rprintf("copy_crash_list_v2 did not work inside insert_edge_in_subnet!\n");
     return -1;
   }
   if(OK(graph_ins_edge)(g,p,adj_node)!=0){
     Rprintf("graph_ins_edge did not work inside insert_edge_in_subnet!\n");
     return -1;
   }
 }
 return 0;
}

/************************ spTree_restricted_v1 ****************************/

int FN(spTree_restricted_v1)(Graph *graph,
			     Graph *graphTree,
			     List *P,
			     PathVertex *start)
{
    AdjList              *adjlist;

    ListElmt             *element1,
                         *element2,
                         *el,
                         *member;

    PathVertex           *Node,
                         *prevNode,
                         *pth_vertex1,
                         *pth_vertex2,
                         *lst_vertex,
                         *adj_vertex;

    int                   retval;

/** Allocate space for all the vertices to be included in "graphTree". **/
/** The number of vertices is same as the list size of Path. **/
/** Path contains all the nodes that are reachable from the start node **/

  for(el=list_head(P);el!=NULL;el=list_next(el)){
    PathVertex *pv_ptr;
    pv_ptr = (PathVertex *)Calloc(1,PathVertex);
    /*if (pv_ptr == NULL){Rprintf("malloc is not working !!");return -1;}*/
    /** It will constitute the vertex list of graphTree **/
    pv_ptr->data = Calloc(1,int);/*if(pv_ptr->data==NULL){Rprintf("malloc is not working !!");return -1;}*/

    lst_vertex = (PathVertex *)list_data(el);

    *((int *)pv_ptr->data) = *((int *)lst_vertex->data);

    pv_ptr->d = lst_vertex->d;
    pv_ptr->parent = lst_vertex->parent;

    OK(list_init)(&pv_ptr->crashList,NULL);

    retval = OK(graph_ins_vertex)(graphTree,pv_ptr);
    if(retval != 0){
      Rprintf("Vertex insertion has failed!\n");
      return -1;
    }
  }

/** Since we have already allocated space for the vertices, we only need to allocate space for their adjacent vertices. **/
/*********************************************************************************************************
 * In order to construct the shortest-path tree, which is the union of the shortest paths from every node*
 * of the network to the "start" vertex, we visit every node of the network using the following for loop.*
 * Subsequently, we enter all edges in graphTree that constitute the shortest path from "start" vertex   *
 * to the network node.                                                                                  *
 ********************************************************************************************************/

  for (element1 = list_head(&graph->adjlists); element1 != NULL; element1 =
	 list_next(element1)){

    pth_vertex1 = ((AdjList *)list_data(element1))->vertex;
    Node = pth_vertex1;

/** Include the shortest path from "Node" to the "start" vertex in "graphTree". **/
/** Until we reach the "start" node from the selected vertex of the network via shortest path route, run the following while loop **/
    while((graphTree->match(Node, start)) != 1){

      if(Node->color != black){break;}

      prevNode = Node->parent;

/** Include the edge prevNode------->Node in graphTree. **/
/** Allocate space for "Node" to be included in the adjacency list of "prevNode". **/
      PathVertex  *adj_node;

      adj_node = (PathVertex *)Calloc(1,PathVertex);
      /*if (adj_node == NULL){Rprintf("malloc is not working !!");return -1;}*/

      adj_node->data = Calloc(1,int);
      /*if (adj_node->data == NULL){Rprintf("malloc is not working !!");return -1;}*/

      *((int *)(adj_node->data)) = *((int *)(Node->data));

      /** Include "weight" and "crashList" inside prevNode------->Node.
	  Search the adjacency list of "prevNode" in graph. **/
      for (element2 = list_head(&graph->adjlists); element2 != NULL;
	   element2 = list_next(element2)){

	adjlist = (AdjList *)list_data(element2);
	pth_vertex2 = ((AdjList *)list_data(element2))->vertex;

	if (graphTree->match(prevNode,pth_vertex2)){
	  break;
	}
      }

      for (member = list_head(&adjlist->adjacent); member != NULL;
	   member = list_next(member)){

	adj_vertex = list_data(member);

	if(graphTree->match(adj_vertex, Node)){
	  break;
	}
      }

      adj_node->weight = adj_vertex->weight;
      adj_node->edgeID = adj_vertex->edgeID;

      OK(list_init)(&(adj_node->crashList),&destroy_crash);
      retval = FN(copy_crash_list_v2)(&(adj_vertex->crashList),
				      &(adj_node->crashList));
      if(retval != 0){
	Rprintf("copy_crash_list_v2 did not work in spTree_restricted!\n");
	return -1;
      }

      retval = OK(graph_ins_edge)(graphTree, pth_vertex2, adj_node);

      if (retval == 0){
	Node = prevNode;
      } else {
        /** Do not require to insert the edges that are already inserted in a previous shortest-path inclusion. **/
        /** If "retval != 0",  the edge has already been inserted, therefore move to the next edge in the shortest-path. **/
        /** Break from the while loop, because all edges in the current shoretest-path that trace back to the "start" vertex has already been inserted. **/
        /** Before breaking out from the loop, free the memory allocated to adj_node,adj_vert and the crashList allotted to adj_node **/
	pathVertex_destroy(adj_node);
	break;
      }
    }
  }
  return 0;
}


