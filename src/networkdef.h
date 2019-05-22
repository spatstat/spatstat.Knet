/*
     networkdef.h

     Data type definitions for linear networks

*/

/*********************************************************************
*                                                                    *
*  Define a structure for linked list elements.                      *
*                                                                    *
*********************************************************************/

typedef struct ListElmt_
{
    void                 *data;
    struct  ListElmt_    *next;
} ListElmt;


/*********************************************************************
*                                                                    *
*  Define a structure for linked lists.                              *
*                                                                    *
*********************************************************************/

typedef struct List_
{

    int                 size;

    int                 (*match)(const void *key1, const void *key2);
    void                (*destroy)(void *data);

    ListElmt            *head;
    ListElmt            *tail;

} List;

#define list_size(list) ((list)->size)

#define list_head(list) ((list)->head)

#define list_tail(list) ((list)->tail)

#define list_is_head(list, element) ((element) == (list)->head ? 1 : 0)

#define list_is_tail(element) ((element)->next ==NULL ? 1 : 0)

#define list_data(element) ((element)->data)

#define list_next(element) ((element)->next)

/*********************************************************************
*                                                                    *
*  Implement sets as linked lists.                                   *
*                                                                    *
*********************************************************************/

typedef List Set;

#define set_destroy list_destroy

#define set_size(set) ((set)->size)

/*********************************************************************
*                                                                    *
*  Define a structure for adjacency lists.                           *
*                                                                    *
*********************************************************************/

typedef struct AdjList_
{

    void         *vertex;
    Set          adjacent;

} AdjList;


/*********************************************************************
*                                                                    *
*  Define a structure for graphs.                                    *
*                                                                    *
*********************************************************************/

typedef struct Graph_
{

    int          vcount;
    int          ecount;
    int          (*match)(const void *key1, const void *key2);
    void         (*destroy)(void *data);

    List         adjlists;
} Graph;

#define graph_adjlists(graph) ((graph)->adjlists)

#define graph_vcount(graph) ((graph)->vcount)

#define graph_ecount(graph) ((graph)->ecount)


/*********************************************************************
*                                                                    *
*  Define colors for vertices in graphs.                             *
*                                                                    *
*********************************************************************/

typedef enum VertexColor_ {white, gray, black} VertexColor;


/*********************************************************************
*                                                                    *
*  Define a structure for the point events on the network            *
*                                                                    *
*********************************************************************/
typedef struct SimpleCrash_
{
	double	tp;
	int	edgeId;
} SimpleCrash;

typedef struct MultipleCrash_
{
	double	tp;
	int	edgeId;
	int     frequency;
} MultipleCrash;

typedef struct WeightedCrash_
{
	double	tp;
	int	edgeId;
	int     frequency;
	double  lambda;
} WeightedCrash;

/*********************************************************************
*                                                                    *
* Define a structure for vertices in shortest-path problems.         *
*                                                                    *
*********************************************************************/
typedef struct PathVertex_
{
  void		*data;
  double	weight;
  VertexColor	color;
  double	d;

  int		edgeID;
  int           degree;
  List		crashList;

  struct PathVertex_   *parent;

} PathVertex;

#define NETWORKDEF_H_
