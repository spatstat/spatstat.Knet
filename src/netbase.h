/*********************************************************************
*     netbase.h                                                      *
*     Declarations for functions defined in netbase.c                *
*********************************************************************/

#ifndef NETWORKDEF_H_
#include "networkdef.h"
#endif

int malloc_not_working(void);

void relax(PathVertex *u, PathVertex *v, double weight);

int pth_vert_compare(PathVertex *pv1,PathVertex *pv2);

void list_init(List *list, void (*destroy)(void *data));

void list_destroy(List *list);

int match_graph(const void *pvPtr1, const void *pvPtr2);

void destroy_crash(void *cr);

void    path_destroy(List *list);

void pathVertex_destroy(void *pth_data);

void graph_init(Graph *graph, int (*match)(const void *key1, const void *key2),
                void (*destroy)(void *data));

void set_init(Set *set, int (*match)(const void *key1, const void *key2), void (*destroy)(void *data));

int set_remove(Set *set, void **data);

int list_ins_next(List *list, ListElmt *element, const void *data);

int graph_ins_vertex(Graph *graph, const void *data);

int set_insert(Set *set, const void *data);


int graph_ins_edge(Graph *graph, const void *data1, const void *data2);


int set_is_member(const Set *set, const void *data);

int graph_rem_vertex(Graph *graph, void **data);

int graph_rem_edge(Graph *graph, void *data1, void **data2);

void graph_destroy(Graph *graph);

int list_rem_next(List *list, ListElmt *element, void **data);


void delete_edge_from_graph(Graph *graph,PathVertex *pth_vert,PathVertex *adj_vert);

int graphEdgeDifference_v2(Graph *small_graph, Graph *large_graph);

int deleteSameEdge(Graph *graph);

int graph_ins_breakpt(Graph *graph,PathVertex *node,PathVertex *brk_pt, int vert_id,int edge_id,double link_dist);

int sorted_list(List *list, List *ord_list);

int ord_list_ins_next_pthVert(List *ord_list,void *data);

int vertex_degree(Graph *graph, PathVertex *pth_ptr);

void create_distance_and_degree_array(List *srtd_lst,double *tme_val,int *degree_val);

void include_vertex_degree_info_inside_pthvert(Graph *graph, List *paths);

void create_m_val_array(double *tme_val,int *degree_val,double *tme_uppr_lmt,int *m_val, int lst_size,int *array_size);

void allot_inv_mvals_in_dist_array(double *dist_ptr,double *inv_mv_ptr,
				   double *inv_mval_vec,
				   double *MAX_Distance,double *MIN_Distance,
				   int *no_of_distance);

void delete_vertex_from_graph(Graph *g,PathVertex *pth_vert);

