/*********************************************************************
*     netbase.c                                                      *
*     Basic functions for computing in networks                      *
*********************************************************************/

#include <R.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>

#ifndef NETWORKDEF_H_
#include "networkdef.h"
#endif

#include "netbase.h"

int malloc_not_working()
{
    Rprintf("malloc is not working !!");
    return -1;
}

/*****************************  relax  ***************************/

void relax(PathVertex *u, PathVertex *v, double weight)
{
    /***************************************************************
    *                                                              *
    * Relax edge between two vertices u and v.                     *
    *                                                              *
    ****************************************************************/
    if(v->d > u->d + weight)
    {
        v->d = u->d + weight;
        v->parent = u;
    }
    return;
}

/***************************** pth_vert_compare  ***************************/

int pth_vert_compare(PathVertex *pv1,PathVertex *pv2){
if((pv2->d)>(pv1->d))
    return 1;
else
    return 0;
}

/*****************************  list_init  ***************************/
void list_init(List *list, void (*destroy)(void *data))
{

    /*********************************************************************
    *                                                                    *
    * Initialize the list.                                               *
    *                                                                    *
    *********************************************************************/

    list->size = 0;
    list->destroy = destroy;
    list->head = NULL;
    list->tail = NULL;

    return;
}

/*****************************  list_destroy  ***************************/

void list_destroy(List *list)
{
    void             *data;

    /*********************************************************************
    *                                                                    *
    * Remove each element.                                               *
    *                                                                    *
    *********************************************************************/

    while (list_size(list) > 0)
    {
        if (list_rem_next(list, NULL, (void **)&data) == 0 && list->destroy
                != NULL)
        {
            /*************************************************************
            *                                                            *
            * Call a user-defined function                               *
	    * to free dynamically allocated data.                        *
            *                                                            *
            **************************************************************/
            list->destroy(data);
        }
    }

    /*************************************************************************
    *                                                                        *
    * No operations are allowed now, but clear the structure as a precaution.*
    *                                                                        *
    *************************************************************************/

    memset(list, 0, sizeof(List));
    return;
}

/*****************************  match_graph  ***************************/

int match_graph(const void *pvPtr1, const void *pvPtr2)
{
    if (*((int *)(((PathVertex *)pvPtr1)->data)) ==
	*((int *)(((PathVertex *)pvPtr2)->data)))
      return 1;
    return 0;
}

/*****************************  destroy_crash  ***************************/

void destroy_crash(void *cr){
    Free(cr);
}

/*****************************  path_destroy  ***************************/

void path_destroy(List *list){
  if (list_size(list) == 0) {
    memset(list,0,sizeof(List));
    return;
  }
  ListElmt    *element,
              *next_element;
  element = list_head(list);
  while(element != NULL){
    next_element = list_next(element);
    Free(element);
    element = next_element;
  }
  memset(list,0,sizeof(List));
}

/**************************  pathVertex_destroy  ***************************/

void pathVertex_destroy(void *pth_data){
  /** Delete the crashList in the path vertex structure **/

  if((&(((PathVertex *)pth_data)->crashList))->destroy != NULL){
        list_destroy(&(((PathVertex *)pth_data)->crashList));
        //Rprintf("Crash List destroyed successfully!\n");
    }

  /** Now free the PathVertex structure **/
  int *id =  ((PathVertex *)pth_data)->data;
  Free(id);
  Free(pth_data);
  return ;
}

/*****************************  graph_init  ***************************/

void graph_init(Graph *graph, int (*match)(const void *key1, const void *key2),
                void (*destroy)(void *data))
{

    /*********************************************************************
    *                                                                    *
    * Initialize the graph.                                              *
    *                                                                    *
    *********************************************************************/

    graph->vcount = 0;
    graph->ecount = 0;
    graph->match = match;
    graph->destroy = destroy;

    /*********************************************************************
    *                                                                    *
    * Initialize the list of adjacency-list structures.                  *
    *                                                                    *
    *********************************************************************/

    list_init(&graph->adjlists, NULL);

    return;
}

/*****************************  set_init  ***************************/

void set_init(Set *set, int (*match)(const void *key1, const void *key2),
              void (*destroy)(void *data))
{

    /*********************************************************************
    *                                                                    *
    * Initialize the set.                                                *
    *                                                                    *
    *********************************************************************/

    list_init(set, destroy);
    set->match = match;

    return;
}

/*****************************  set_remove  ***************************/

int set_remove(Set *set, void **data)
{

    ListElmt              *member,
                          *prev;

    /*********************************************************************
    *                                                                    *
    * Find the member to remove.                                         *
    *                                                                    *
    *********************************************************************/

    prev = NULL;

    for (member = list_head(set); member != NULL; member = list_next(member))
    {
        if (set->match(*data, list_data(member)))
            break;

        prev = member;

    }

    /*********************************************************************
    *                                                                    *
    * Return if the member was not found.                                *
    *                                                                    *
    *********************************************************************/

    if ( member == NULL)
    {
        return -1;
    }


    /*********************************************************************
    *                                                                    *
    * Remove the member.                                                 *
    *                                                                    *
    *********************************************************************/

    return list_rem_next(set, prev, data);

}

/*****************************  list_ins_next  ***************************/

int list_ins_next(List *list, ListElmt *element, const void *data)
{

    ListElmt                 *new_element;

    /*************************************************************************
    *                                                                        *
    * Allocate storage for the element.                                      *
    *                                                                        *
    *************************************************************************/
    new_element = (ListElmt *)Calloc(1,ListElmt);

    /*if ((new_element = (ListElmt *)Calloc(1,ListElmt)) == NULL)*/
    /*{*/

         /*Rprintf("malloc is not working !!");*/
         /*return -1;*/

    /*}*/
    /*************************************************************************
    *                                                                        *
    * Insert the element into the list.                                      *
    *                                                                        *
    *************************************************************************/

    new_element->data = (void *)data;

    if (element == NULL)
    {

        /*************************************************************************
        *                                                                        *
        * Handle insertion at the head of the list.                              *
        *                                                                        *
        *************************************************************************/
        if (list_size(list) == 0)
            list->tail = new_element;

        new_element->next = list->head;
        list->head = new_element;

    }
    else
    {

        /*************************************************************************
        *                                                                        *
        * Handle insertion somewhere other than at the head.                     *
        *                                                                        *
        *************************************************************************/
        if (element->next == NULL)
            list->tail = new_element;

        new_element->next = element->next;
        element->next = new_element;

    }

    /*************************************************************************
    *                                                                        *
    * Adjust the size of the list to account for the inserted element.       *
    *                                                                        *
    *************************************************************************/

    list->size++ ;

    return 0;

}

/*****************************  graph_ins_vertex  ***************************/

int graph_ins_vertex(Graph *graph, const void *data)
{

    ListElmt            *element;

    AdjList             *adjlist;

    int                  retval;

    /****************************************************************************
    *                                                                           *
    * Do not allow the insertion of duplicate vertices.                         *
    *                                                                           *
    ****************************************************************************/

    for (element = list_head(&(graph->adjlists)); element != NULL; element =
                list_next(element))
    {

        if (graph->match(data, ((AdjList *)list_data(element))->vertex))  /* Compare the data members of PathVertex structures */

            return 1;
    }

    /****************************************************************************
    *                                                                           *
    * Insert the vertex                                                         *
    *                                                                           *
    ****************************************************************************/

    if ((adjlist = (AdjList *)Calloc(1,AdjList)) == NULL)
        return -1;


    (adjlist->vertex) = (void *)data;

    set_init(&(adjlist->adjacent), graph->match, graph->destroy);

    if ((retval = list_ins_next(&graph->adjlists, list_tail(&graph->adjlists),
                                adjlist)) != 0)
    {

        return retval;
    }

    /****************************************************************************
    *                                                                           *
    * Adjust the vertex count to account for the inserted vertex.               *
    *                                                                           *
    ****************************************************************************/

    graph->vcount++;

    return 0;

}

/*****************************  set_insert  ***************************/

int set_insert(Set *set, const void *data)
{

    /*********************************************************************
    *                                                                    *
    * Do not allow the insertion of duplicates.                          *
    *                                                                    *
    *********************************************************************/

    if (set_is_member(set, data))
    {

        return 1;

    }

    /*********************************************************************
    *                                                                    *
    * Insert the data.                                                   *
    *                                                                    *
    *********************************************************************/

    return list_ins_next(set, list_tail(set), data);

}

/*****************************  graph_ins_edge  ***************************/

int graph_ins_edge(Graph *graph, const void *data1, const void *data2)
{

    ListElmt                *element;

    int                     retval;


    /****************************************************************************
    *                                                                           *
    * Do not allow insertion of an edge without both its vertices in the graph. *
    *                                                                           *
    ****************************************************************************/

    for (element = list_head(&graph->adjlists); element != NULL; element =
                list_next(element))
    {

        if (graph->match(data2, ((AdjList *)list_data(element))->vertex))
            break;
    }

    if (element == NULL) {return -1;}

    for (element = list_head(&graph->adjlists); element != NULL; element =
                list_next(element))
    {

        if (graph->match(data1, ((AdjList *)list_data(element))->vertex))
            break;
    }

    if (element == NULL){return -1;}

    /****************************************************************************
    *                                                                           *
    * Insert the second vertex into the adjacency list of the first vertex.     *
    *                                                                           *
    ****************************************************************************/

    if ((retval = set_insert(&((AdjList *)list_data(element))->adjacent, data2))
            != 0)
    {

        return retval;

    }


    /****************************************************************************
    *                                                                           *
    * Adjust the edge count to account for the inserted edge.                   *
    *                                                                           *
    ****************************************************************************/

    graph->ecount++;

    return 0;

}

/*****************************  set_is_member  ***************************/

int set_is_member(const Set *set, const void *data)
{

    ListElmt                 *member;

    /*********************************************************************
    *                                                                    *
    * Determine if the data is a member of the set.                      *
    *                                                                    *
    *********************************************************************/

    for (member = list_head(set); member != NULL; member = list_next(member))
    {

        if (set->match(data, list_data(member)))
            return 1;

    }

    return 0;

}

/*****************************  graph_rem_vertex  ***************************/

int graph_rem_vertex(Graph *graph, void **data)
{

    ListElmt                 *element,
                             *temp,
                             *prev;

    AdjList                  *adjlist;

    int                      found;


    /****************************************************************************
    *                                                                           *
    * Traverse each adjacency list and the vertices it contains.                *
    *                                                                           *
    ****************************************************************************/

    prev = NULL;
    found = 0;

    for (element = list_head(&graph->adjlists); element != NULL; element =
                list_next(element))
    {

        /****************************************************************************
        *                                                                           *
        * Do not allow removal of the vertex if it is in an adjacency list.         *
        *                                                                           *
        ****************************************************************************/

        if (set_is_member(&((AdjList *)list_data(element))->adjacent, *data))
            return -1;

        /****************************************************************************
        *                                                                           *
        * Keep a pointer to the vertex to be removed.                               *
        *                                                                           *
        ****************************************************************************/

        if (graph->match(*data, ((AdjList *)list_data(element))->vertex))
        {

            temp = element;
            found = 1;
        }

        /****************************************************************************
        *                                                                           *
        * Keep a pointer to the vertex before the vertex to be removed.             *
        *                                                                           *
        ****************************************************************************/

        if (!found)
            prev = element;

    }

    /****************************************************************************
    *                                                                           *
    * Return if the vertex was not found.                                       *
    *                                                                           *
    ****************************************************************************/

    if (!found)
        return -1;

    /****************************************************************************
    *                                                                           *
    * Do not allow removal of the vertex if its adjacency list is not empty.    *
    *                                                                           *
    ****************************************************************************/

    if (set_size(&((AdjList *)list_data(temp))->adjacent) > 0)
        return -1;

    /****************************************************************************
    *                                                                           *
    * Remove the vertex.                                                        *
    *                                                                           *
    ****************************************************************************/

    if (list_rem_next(&graph->adjlists, prev, (void **)&adjlist) != 0)
        return -1;

    /****************************************************************************
    *                                                                           *
    * Free the storage allocated by the abstract datatype.                      *
    *                                                                           *
    ****************************************************************************/

    *data = adjlist->vertex;
    Free(adjlist);

    /****************************************************************************
    *                                                                           *
    * Adjust the vertex count to account for the removed vertex.                *
    *                                                                           *
    ****************************************************************************/

    graph->vcount--;

    return 0;

}

/*****************************  graph_rem_edge  ***************************/

int graph_rem_edge(Graph *graph, void *data1, void **data2)
{

    ListElmt         *element;

    for (element = list_head(&graph->adjlists); element != NULL; element =
                list_next(element))
    {

        if (graph->match(data1, ((AdjList *)list_data(element))->vertex))
            break;
    }

    if (element == NULL){

        return -1;
    }

    /****************************************************************************
    *                                                                           *
    * Remove the second vertex from the adjacency list of the first vertex.     *
    *                                                                           *
    ****************************************************************************/

    if (set_remove(&((AdjList *)list_data(element))->adjacent, data2) != 0)
    {

            return -1;
    }



    /****************************************************************************
    *                                                                           *
    * Adjust the edge count to account for the removed edge.                    *
    *                                                                           *
    ****************************************************************************/

    graph->ecount--;

    return 0;

}

/*****************************  graph_destroy ***************************/

void graph_destroy(Graph *graph)
{

    AdjList         *adjlist;

    /****************************************************************************
    *                                                                           *
    * Remove each adjacency-list structure and destroy its adjacency list.      *
    *                                                                           *
    ****************************************************************************/

    while (list_size(&graph->adjlists) > 0)
    {

        if (list_rem_next(&graph->adjlists, NULL, (void **)&adjlist) == 0)
        {

            set_destroy(&adjlist->adjacent);

            if (graph->destroy != NULL){

                graph->destroy(adjlist->vertex);

            }
            Free(adjlist);

        }
    }

    /****************************************************************************
    *                                                                           *
    * Destroy the list of adjacency-list structures, which is now empty.        *
    *                                                                           *
    ****************************************************************************/

    list_destroy(&graph->adjlists);

    /****************************************************************************
    *                                                                           *
    * No operations are allowed now, but clear the structure as a precaution    *
    *                                                                           *
    ****************************************************************************/

    memset(graph, 0, sizeof(Graph));

    return;

}



/**************************** list_rem_next ***************************************/

int list_rem_next(List *list, ListElmt *element, void **data)
{

    ListElmt              *old_element;


    /*************************************************************************
    *                                                                        *
    * Do not allow removal from an empty list.                               *
    *                                                                        *
    *************************************************************************/

    if(list_size(list) == 0)
        return -1;

    /*************************************************************************
    *                                                                        *
    * Remove the element from the list.                                      *
    *                                                                        *
    *************************************************************************/

    if ( element == NULL)
    {

        /*************************************************************************
        *                                                                        *
        * Handle removal from the head of the list.                               *
        *                                                                        *
        *************************************************************************/
        *data = list->head->data;
        old_element = list->head;
        list->head = list->head->next;

        if (list_size(list) == 1)
            list->tail = NULL;

    }

    else
    {

        /*************************************************************************
        *                                                                        *
        * Handle removal from somewhere other than the head.                     *
        *                                                                        *
        *************************************************************************/
        if (element->next == NULL)
            return -1;

        *data = element->next->data;
        old_element = element->next;
        element->next = element->next->next;

        if (element->next == NULL)
            list->tail = element;

    }

    /*************************************************************************
    *                                                                        *
    * Free the storage allocated by the abstract datatype.                   *
    *                                                                        *
    *************************************************************************/

    Free(old_element);

    /*************************************************************************
    *                                                                        *
    * Adjust the size of the list to account for the removed element.        *
    *                                                                        *
    *************************************************************************/

    list->size--;

    return 0;

}




/***************************** graphEdgeDifference_v2 ************************************/

int graphEdgeDifference_v2(Graph *small_graph, Graph *large_graph)
{
/*********************************************************************************
** This function removes all edges from the road network (large_graph) that are **
** not present in the shortest path tree (small_graph). This in turn converts   **
** the large graph into a graph containing all edges that are absent in the     **
** shortest-path tree.                                                          **
*********************************************************************************/
    ListElmt    *element_S,*member_S;

    AdjList     *adjlist_S;

    PathVertex  *pth_vertex_S,*adj_vertex_S;

    int         retval=0;
/************************************************************************************
* The small_graph is a directed graph. Further, each edge appears only once in the  *
* shortest-path tree (small_graph). However, we need to remove any chosen edge twice*
* from the large_graph because it is an undirected graph and each edge is recorded  *
* twice within it. This works perfectly for every edge in the small_graph except the*
* edges emanating from the crash point (added as a vertex in the graph). The reason *
* for this is that we have included edges from the crash point as directed edges in *
* the large_graph. This complicacy is handled by the "delete_edge_from_graph" funct *
* -ion at the end of this function.                                                 *
************************************************************************************/

/** Start selecting edges in the shortest-path tree (small_graph) one by one.**/
element_S = list_head(&small_graph->adjlists);

    if(element_S == NULL) /** Checking if the small_graph is empty **/
    {
        Rprintf("Small graph is empty.\n");
        retval = -1;
        return retval;
    }

    while(element_S != NULL) /** Loop through every vertex of the small_graph **/
    {
        adjlist_S = ((AdjList *)list_data(element_S));
        pth_vertex_S = (PathVertex *)(adjlist_S->vertex);

            member_S = list_head(&adjlist_S->adjacent);
            if(member_S == NULL) /** Checking if there is any edge out from pth_vertex_S **/
            {

                element_S = list_next(element_S); /** Go to the next vertex in the graph **/
                continue;
            }
            else
            {
                while(member_S != NULL)
                {
                    adj_vertex_S = list_data(member_S);

                    /** Delete the edge pth_vertex------>adj_vertex and free up the memory if the edge removal is successful! **/

                    delete_edge_from_graph(large_graph,pth_vertex_S,adj_vertex_S);

                    /** Delete the same edge recorded as adj_vertex------>pth_vertex and free up the memory if the edge removal is successful! **/
                    delete_edge_from_graph(large_graph,adj_vertex_S,pth_vertex_S);

                    member_S = list_next(member_S);
                }
            }

        element_S = list_next(element_S);
    }


return retval;
}

/***************************** delete_edge_from_graph ************************************/

void delete_edge_from_graph(Graph *graph,PathVertex *pth_vert,PathVertex *adj_vert){

/***************************************************************************************
**       This function deletes the edge (pth_vert------>adj_vert) from the graph.     **
***************************************************************************************/
void *data1, *vd_ptr, **data2;

int retval;

data1 = pth_vert;
vd_ptr = adj_vert;
data2 = &vd_ptr;

retval = graph_rem_edge(graph,data1,data2);

if (retval == 0){
    pathVertex_destroy(*data2);
    return;
}

return;
}

/***************************** deleteSameEdge ************************************/

int deleteSameEdge(Graph *graph){
ListElmt *element,*member;
AdjList *adjlist;
PathVertex *pth_vertex,*adj_vertex;

    element = list_head(&graph->adjlists);
    if(element==NULL){Rprintf("Graph is empty for deleting same edges.\n");return -1;}

    while(element != NULL){
        adjlist = (AdjList *)list_data(element);
        pth_vertex = (PathVertex *)(adjlist->vertex);

        /** Check whether there is any edge emanating from the pth_vertex **/
        member = list_head(&adjlist->adjacent);

        /** If there are no edges, move to the next vertex in the graph **/
        if(member != NULL){
            while(member != NULL){
                adj_vertex = list_data(member);
                delete_edge_from_graph(graph,adj_vertex,pth_vertex);
                member = list_next(member);
            }
        }
        element = list_next(element);
    }

return 0;
}


/***************************** graph_ins_breakpt ************************************/

int graph_ins_breakpt(Graph *graph,PathVertex *node,PathVertex *brk_pt,
                      int vert_id,int edge_id,double link_dist){
/** This function will be called inside the function "extended_sh_path_tree" if link_dist != 0 **/
int retval;

node->data = Calloc(1,int); /*if(node->data == NULL){Rprintf("malloc is not working !!");return -1;}*/ /** allocate space to store vertex id **/

*((int *)node->data) = vert_id;

list_init(&node->crashList,NULL);

retval = graph_ins_vertex(graph,node);

/** Create the break point and add the edge informations in order to create (pth_vertex--->brk_pt1) **/
brk_pt->data = Calloc(1,int); /*if(brk_pt->data == NULL){Rprintf("malloc is not working !!");return -1;}*/

*((int *)brk_pt->data) = vert_id;

brk_pt->weight = link_dist;

brk_pt->edgeID = edge_id;

/** Initialize the crash lists associated with the linking edges **/
list_init(&(brk_pt->crashList), &destroy_crash);

if(retval == 0) return retval;

return 0;
}

/***************************** sorted_list ************************************/

int sorted_list(List *list, List *ord_list){
  ListElmt   *element;
  PathVertex *pth_ptr;
  int retval;

  for(element=list_head(list);element!=NULL;element=list_next(element)){
    pth_ptr = (PathVertex *)list_data(element);
    retval = ord_list_ins_next_pthVert(ord_list, (void *)pth_ptr);
    if(retval != 0){
      Rprintf("ord_list_ins_next did not work!\n");
      return -1;
    }
  }
  return 0;
}

/******************** ord_list_ins_next_pthVert ****************************/

int ord_list_ins_next_pthVert(List *ord_list,void *data){
  ListElmt *element,*prev_element;
  PathVertex *pth_ptr;
  int retval;
  if(list_size(ord_list)==0){
    retval=list_ins_next(ord_list,NULL,data);
    if(retval!=0){
      Rprintf("list_ins_next did not work inside ord_list_ins_next_pthVert.\n");
      return -1;}
  } else {
    prev_element=NULL;
    element=list_head(ord_list);
    pth_ptr = (PathVertex *)(element->data);

    while(element!=NULL){
      pth_ptr = (PathVertex *)element->data;
      if(pth_vert_compare((PathVertex *)data,pth_ptr)==1){
	break;
      } else {
	prev_element=element;
	element=element->next;
      }
    }
    retval=list_ins_next(ord_list,prev_element,data);
    if(retval !=0){
      Rprintf("list_ins_next did not work inside ord_list_ins_next.\n");
      return -1;
    }
  }
  return 0;
}

/*********************** vertex_degree **********************************/

int vertex_degree(Graph *graph, PathVertex *pth_ptr){

  ListElmt    *element;
  AdjList     *adjlist;
  PathVertex  *pth_vertex;
  int         count=0,retval;

  for(element=list_head(&graph->adjlists);
      element!=NULL;
      element=list_next(element)){

    adjlist = (AdjList *)list_data(element);
    pth_vertex = adjlist->vertex;

    if(graph->match(pth_vertex,pth_ptr))
      break;
  }

  if(element==NULL) return -1;

  count = list_size(&adjlist->adjacent);

  if(pth_ptr->d == 0){
    retval = (count);
  } else{
    retval = (count+1);
  }
  return retval;
}

/******* create_distance_and_degree_array ********************************/

void create_distance_and_degree_array(List *srtd_lst,
				      double *tme_val,
				      int *degree_val){
  /** Space for tme_val and m_val had already been allocated.
      The array size in both occasions are equal to the size of the srtd lst.
      This function simply fill these two arrays with corresponding
      time values and degree values.
  **/

  ListElmt   *element;
  PathVertex *pth_vertex;

  int     j, lst_sze = list_size(srtd_lst);

  element = list_head(srtd_lst);

  for(j=0;j<lst_sze;++j){
    pth_vertex = (PathVertex *)list_data(element);
    tme_val[j] = pth_vertex->d;
    degree_val[j] = pth_vertex->degree;
    element = list_next(element);
  }
  return;
}

/************ include_vertex_degree_info_inside_pthvert *********************/

void include_vertex_degree_info_inside_pthvert(Graph *graph, List *paths){

  ListElmt    *element;
  PathVertex  *pth_vertex;
  int retval;

  for(element=list_head(paths);element!=NULL;element=list_next(element)){
    pth_vertex = (PathVertex *)list_data(element);
    retval = vertex_degree(graph,pth_vertex);
    pth_vertex->degree = retval;
  }
  return;
}

/*************** create_m_val_array ************************************/

void create_m_val_array(double *tme_val,
			int *degree_val,
			double *tme_uppr_lmt,
			int *m_val,
                        int lst_size,
			int *array_size){
  int j = 0,
      k = 0,
      m_count = 0,
      delta_not = degree_val[0]; /** degree of the "start" node.
				     In cases where no crashes had occurred
				     on the node, delta_not = 2 **/

    /** tme_val[0] will always be 0.
	Therefore, tme_uppr_lmt[0] = tme_val[1] and m_val[0] = degree_val[0] **/
  ++j;
  tme_uppr_lmt[k] = tme_val[j];
  m_val[k] = delta_not;

  --lst_size;

  while(j < lst_size){
    if(tme_val[j] < tme_val[j+1]){
      /** For every distance that is greater than the previous entry
	  in the sorted distance array, we have a new entry in the
	  time upper limit and corresponding m-value.
      **/
      ++k;
      tme_uppr_lmt[k] = tme_val[j+1];
      m_count += (degree_val[j] - 2);
      m_val[k] = m_count + delta_not;
    } else {
      m_count += (degree_val[j]-2);
    }
    ++j;
  }

  *array_size = (k+1);
  return;
}


/***************** allot_inv_mvals_in_dist_array ****************************/
void allot_inv_mvals_in_dist_array(double *dist_ptr,
				   double *inv_mv_ptr,
				   double *inv_mval_vec,
				   double *MAX_Distance,
				   double *MIN_Distance,
				   int *no_of_distance){
  int index = 0;
  double epsilon = 0.00000001;
  double distance = *dist_ptr;
  double dist = (distance - *MIN_Distance - epsilon);
  double interval = (*MAX_Distance - *MIN_Distance)/(*no_of_distance - 1);
  if(dist >= 0){
    index = (int)(dist/interval) + 1;
  }
  inv_mval_vec[index] = inv_mval_vec[index] + *(inv_mv_ptr);
  return;
}

/******************* delete_vertex_from_graph *****************************/

void delete_vertex_from_graph(Graph *g,PathVertex *pth_vert){
  ListElmt        *element,*member;
  AdjList         *adjlist;
  PathVertex      *pth_vertex,*adj_vertex;
  int retval;

  for(element = list_head(&g->adjlists);
      element != NULL;
      element=list_next(element)){
    adjlist = (AdjList *)(list_data(element));
    pth_vertex = (PathVertex *)(adjlist->vertex);
    if(g->match(pth_vert,pth_vertex)) break;
  }

  /** First delete adjacent edges of pth_vert if there are any **/
  void *data1; data1 = pth_vertex;

  member=list_head(&adjlist->adjacent);

  while(member!=NULL){
    adj_vertex = (PathVertex *)(list_data(member));
    void *vd_ptr,**data2;
    vd_ptr = adj_vertex;
    data2 = &vd_ptr;
    member = list_next(member);
    retval = graph_rem_edge(g,data1,data2);
    if(retval == 0){
      pathVertex_destroy(*data2);
    } else {
      Rprintf("graph_rem_edge did not work!\n");
      return;
    }
    //delete_edge_from_graph(g,pth_vertex,adj_vertex);
    //pathVertex_destroy(adj_vertex);
    //Rprintf("Adjacent edge deletion is successful!\n");
  }

  /** Now delete the vertex **/

  void    **data22;

  //data1 = pth_vertex;
  data22 = &data1;

  if(graph_rem_vertex(g,data22)!=0){
    Rprintf("graph_rem_vertex did not work\n");
    return;
  } else {
    pathVertex_destroy(*data22);
    //Rprintf("Node deletion is successful!\n");
  }
  return;
}

