#include <Python.h>
#include <zlib.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <assert.h>

#include "kseq.h"

KSEQ_INIT(gzFile, gzread)

#define GRAPHEXT ".gfasta"
#define GRAPHEXTGZ ".gfasta.gz"

#if LARGEIDX == 1
	#include <divsufsort64.h>
	typedef int64_t idx_t;
	typedef uint64_t uidx_t;
	#define PYIDXTYPE "L"
	#define PYSONAME "GSA_64"
	#define PYIDXFORMAT "%lld"
	#define PYUIDXFORMAT "%llu"
#else
	#include <divsufsort.h>
	typedef int32_t idx_t;
	typedef uint32_t uidx_t;
	#define PYIDXTYPE "i"
	#define PYSONAME "GSA"
	#define PYIDXFORMAT "%d"
	#define PYUIDXFORMAT "%u"
#endif

static PyObject *GSAError;

typedef struct {
    PyObject_HEAD
    unsigned int vi;
    unsigned int ei;
    PyObject* vertices; 	//Dict of vertices
    PyObject* edges; 		//set of graph_Edge type objects TODO: check if this can be removed!?
	PyObject* origins;		//set containing filenames of the input samples
} Graph;

typedef struct {
    PyObject_HEAD
    unsigned int id;
	int input_origin; 		//field for faster comparison, corresponds to the index of the graph to be aligned
    idx_t contig_start;
    idx_t contig_end;
	idx_t saoffset;
    idx_t rcsaoffset;
	int indexed;				//whether the primary (0) the reverse complement (1) or both (2) of the sequences are indexed
	PyObject *edges_from; 		//set of graph_Edge type objects pointing from this node
    PyObject *edges_to; 		//set of graph_Edge type objects pointing to this node
	PyObject *attributes; 		//additional attributes (only for flexibility). Essential attributes should be in c struct, (TODO) should be removed at some point...
	PyObject *originid_set; 	//set of originids pointing to graph->samples
	PyObject *contigid_set; 	//set of contigids pointing to graph->contigs
	PyObject *coord_origin; 	//id pointing to the file which was used for the coordinates used
	PyObject *coord_contig; 	//id pointing to the contig that was used for the coordinates used
} Vertex;

typedef struct {
    PyObject_HEAD
    Vertex* source;
    Vertex* target;
    uint8_t orientation; 	//same=0, inny=1, outty=2
} Edge;

typedef struct {
    PyObject_HEAD
    idx_t *SA;				//suffix array, has to be signed for libdivsufsort...
    idx_t *SAi; 			//inverse of SA
    unsigned int * LCP;		//LCP array
    unsigned char * T;		//the text indexed in the suffix array
    Vertex **TIiv; 			//array of Vertex pointers for every suffix in inverse SA order
    idx_t n; 				//number of characters in the suffix array
    idx_t sep; 			//indicates the position of the sentinel seperating the two input sets
    idx_t orgn; 			//size of the original suffix array, used by indices extracted from the initial index
    uint8_t rcindex; 		//size of the original suffix array, used by indices extracted from the initial index
    Graph *graph;			//graph object associated with the index (for now only root index has a graph associated with it)
    PyObject * root; 		//reference to GSA type object that points to the root index (for access to TIiv and SAi in sub-indexes)
} GSA;

static Vertex * add_vertex_c(Graph *self, int input_origin, int indexed, char * coord_origin, char * coord_contig, idx_t contig_start, idx_t contig_end, idx_t saoffset, idx_t rcsaoffset);
static Edge * add_edge_c(Graph *self, Vertex* source, Vertex* target, uint8_t orientation);
static int buildGSA(GSA *self, char *input1, char* input2);
static PyObject* split(GSA *self, PyObject *args);
static GSA* extract(GSA *self, PyObject *args);
static PyObject* relabel(GSA *self, PyObject *args);
static PyObject* update(GSA* self, PyObject* args, PyObject* kwds);
static PyObject* updateVertex(GSA* self, PyObject* args, PyObject* kwds);
static PyObject * Vertex_new(PyTypeObject *type, PyObject *args, PyObject *kwds);
static idx_t loadFile(GSA *self, char * input1, int inputid, idx_t f, int useRC);


char comp[128]; //reverse complement translation table

static void rc(char* T, int l) {
    char tmp;
    int i;
    for (i=0; i<((l+1)/2); i++)
    {
		if (i==l-i-1){ 
			T[i]=comp[(int)T[i]];
			break;
		}
    	tmp=T[i];
        T[i]=comp[(int)T[l-i-1]];
		T[l-i-1]=comp[(int)tmp];
    }
}

int endsWith (char* base, char* str) {
    int blen = strlen(base);
    int slen = strlen(str);
    return (blen >= slen) && (0 == strcmp(base + blen - slen, str));
}



















static void
Vertex_dealloc(Vertex* self)
{
    Py_DECREF(self->edges_to);
    Py_DECREF(self->edges_from);
    Py_DECREF(self->attributes);
	Py_DECREF(self->coord_origin);
	Py_DECREF(self->coord_contig);
	Py_DECREF(self->contigid_set);
	Py_DECREF(self->originid_set);
    
    self->ob_type->tp_free((PyObject*)self);
}

static int
Vertex_init(Vertex *self, PyObject *args, PyObject *kwds)
{
    if (! PyArg_ParseTuple(args, "i", &self->id)) {
        PyErr_SetString(PyExc_TypeError, "Cannot initialize Vertex with the supplied arguments");
        return -1;
    }
    
    return 0;
}


static int
Vertex_init_c(Vertex *self, unsigned int id, int input_origin, int indexed, char * coord_origin, char * coord_contig, idx_t contig_start, idx_t contig_end, idx_t saoffset, idx_t rcsaoffset)
{
    self->id=id;
   	
	self->input_origin=input_origin;
	self->indexed=indexed;
	
	PyObject* co=Py_BuildValue("s",coord_origin);
	self->coord_origin=co;
	
	PyObject* cc=Py_BuildValue("s",coord_contig);
	self->coord_contig=cc;

    self->contig_start=contig_start;
    
    self->contig_end=contig_end;
    
    self->saoffset=saoffset;

    self->rcsaoffset=rcsaoffset;

    return 0;
}


//creates a new vertex object and increases the reference count
static PyObject *
Vertex_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
    Vertex *self;
    
    self = (Vertex *)type->tp_alloc(type, 0);
    
    if (self != NULL) {
        
        self->attributes=PyDict_New();
        if (self->attributes == NULL)
        {
            Py_DECREF(self);
            return NULL;
        }
        
        self->edges_from=PySet_New(NULL);
        if (self->edges_from == NULL)
        {
            Py_DECREF(self);
            return NULL;
        }
        
        self->edges_to=PySet_New(NULL);
        if (self->edges_to == NULL)
        {
            Py_DECREF(self);
            return NULL;
        }

		self->originid_set=PySet_New(NULL);
        if (self->originid_set == NULL)
        {
            Py_DECREF(self);
            return NULL;
        }

		self->contigid_set=PySet_New(NULL);
        if (self->contigid_set == NULL)
        {
            Py_DECREF(self);
            return NULL;
        }

        self->id=0;
        self->contig_start=0;
        self->contig_end=0;
		self->saoffset=0;
		self->rcsaoffset=0;
		self->input_origin=0;
		self->coord_origin=NULL;
		self->coord_contig=NULL;
    }
    
    return (PyObject *)self;
}

static PyObject * flip_orientation(Vertex *self, PyObject *args, PyObject *kwds) {
   
	PyObject* iterto;
	PyObject* iterfrom;
	Edge* edge;

   iterfrom=PyObject_GetIter(self->edges_from);
   if (iterfrom==NULL){
	  PyErr_SetString(PyExc_TypeError, "Unable to iterate over set of edges_from");
	  return NULL;
   }
   iterto=PyObject_GetIter(self->edges_to);
   if (iterto==NULL){
	  PyErr_SetString(PyExc_TypeError, "Unable to iterate over set of edges_to");
	  return NULL;
   }
   
   while ( (edge=(Edge*)PyIter_Next(iterfrom)) ) { //all nodes to the right
	  if (edge->orientation==1){
		 if (edge->source->id==self->id){
			edge->orientation=0;
			edge->source=edge->target;
			edge->target=self;
   		 } else {
			assert(edge->target->id==self->id);
			edge->orientation=0;
		 }
	  } else {
		 assert(edge->orientation==0);
		 //orientation == 0, source and target have no meaning so just change orientation to 1 or 2 depending on current source/target
		 if (edge->source->id==self->id){
			edge->orientation=2;
		 } else {
			assert(edge->target->id==self->id);
			edge->orientation=1;
		 }
	  }
	  Py_DECREF(edge);
   }

   while ( (edge=(Edge*)PyIter_Next(iterto)) ) { //all nodes to the right
	  if (edge->orientation==2){
		 if (edge->source->id==self->id){
			edge->orientation=0;
		 } else {
			assert(edge->target->id==self->id);
			edge->orientation=0;
			edge->target=edge->source;
			edge->source=self;
		 }
	  } else {
		 assert(edge->orientation==0);
		 //orientation == 0, source and target have no meaning so just change orientation to 1 or 2 depending on current source/target
		 if (edge->source->id==self->id){
			edge->orientation=2;
		 } else {
			assert(edge->target->id==self->id);
			edge->orientation=1;
		 }
	  }
	  Py_DECREF(edge);
   }
   
   //flip edges_from with edges_to
   PyObject* tmpe;
   tmpe=self->edges_from;
   self->edges_from=self->edges_to;
   self->edges_to=tmpe;
  
   //flip saoffset with rcsaoffset
   idx_t tmpsa;
   tmpsa=self->saoffset;
   self->saoffset=self->rcsaoffset;
   self->rcsaoffset=tmpsa;
   
   Py_DECREF(iterfrom);
   Py_DECREF(iterto);

   Py_INCREF(Py_None);
   return Py_None;
}

static PyMethodDef Vertex_methods[] = {
    {"flip_orientation", (PyCFunction) flip_orientation, METH_VARARGS, "Flips the orientation of a vertex to the graph."},
    {NULL}  /* Sentinel */
};


static PyObject *
Vertex_getid(Vertex *self, void *closure)
{
    PyObject * id=Py_BuildValue("i",self->id);
    return id;
}

static PyObject *
Vertex_getorigin(Vertex *self, void *closure)
{
  	Py_INCREF(self->originid_set);
    return self->originid_set;
}

static PyObject *
Vertex_getinputorigin(Vertex *self, void *closure)
{
    PyObject * org=Py_BuildValue("i",self->input_origin);
    return org;
}




static PyObject *
Vertex_getindexed(Vertex *self, void *closure)
{
    PyObject * idx=Py_BuildValue("i",self->indexed);
    return idx;
}
static int
Vertex_setindexed(Vertex *self, PyObject *args, void *closure)
{
	if (args== NULL) {
		PyErr_SetString(PyExc_TypeError, "Cannot set value to NULL");
		return -1;
	}

	if (! PyInt_Check(args)) {
		PyErr_SetString(PyExc_TypeError,
				"Must be an integer in the range 0 to 2");
		return -1;
	}
	
	
	self->indexed=PyInt_AS_LONG(args);

	return 0;
}






static PyObject *
Vertex_getcontigorigin(Vertex *self, void *closure)
{
	Py_INCREF(self->contigid_set);
    return self->contigid_set;
}

static PyObject *
Vertex_getcoordcontig(Vertex *self, void *closure)
{
    Py_INCREF(self->coord_contig);
    return self->coord_contig;
}
static PyObject *
Vertex_getcoordorigin(Vertex *self, void *closure)
{
    Py_INCREF(self->coord_origin);
    return self->coord_origin;
}
static PyObject *
Vertex_getcontigstart(Vertex *self, void *closure)
{
    PyObject * cs=Py_BuildValue(PYIDXTYPE,self->contig_start);
    return cs;
}
static PyObject *
Vertex_getcontigend(Vertex *self, void *closure)
{
    PyObject * ce=Py_BuildValue(PYIDXTYPE,self->contig_end);
    return ce;
}

static PyObject *
Vertex_getsaoffset(Vertex *self, void *closure)
{
    PyObject * saoffset=Py_BuildValue(PYIDXTYPE,self->saoffset);
    return saoffset;
}
static PyObject *
Vertex_getrcsaoffset(Vertex *self, void *closure)
{
    PyObject * rcsaoffset=Py_BuildValue(PYIDXTYPE,self->rcsaoffset);
    return rcsaoffset;
}

static PyObject *
Vertex_getedgesto(Vertex *self, void *closure)
{
    Py_INCREF(self->edges_to);
    return self->edges_to;
}

static PyObject *
Vertex_getedgesfrom(Vertex *self, void *closure)
{
    Py_INCREF(self->edges_from);
    return self->edges_from;
}

static PyObject *
Vertex_getattributes(Vertex *self, void *closure)
{
    Py_INCREF(self->attributes);
    return self->attributes;
}


static PyGetSetDef Vertex_getseters[] = {
    {"id",
        (getter)Vertex_getid, NULL,
        "vertex id",
        NULL},
    {"origin",
        (getter)Vertex_getorigin, NULL,
        "set of filenames this sequence originates from",
        NULL},
    {"input_origin",
        (getter)Vertex_getinputorigin, NULL,
        "id of the file this sequence originates from (0 or 1 for now)",
        NULL},

	{"indexed",
        (getter)Vertex_getindexed, (setter)Vertex_setindexed,
        "wheter the primary(0), the reverse complement (1) or both (2) sequences defined on this vertex are indexed",
        NULL},

    {"contig_origin",
        (getter)Vertex_getcontigorigin, NULL,
        "set of contig names this sequence originates from",
        NULL},
    {"coord_origin",
        (getter)Vertex_getcoordorigin, NULL,
        "Name of the file in which the contig resides that is used for locating the sequence",
        NULL},
    {"coord_contig",
        (getter)Vertex_getcoordcontig, NULL,
        "Name of the contig used for locating the sequence",
        NULL},
    {"contig_start",
        (getter)Vertex_getcontigstart, NULL,
        "start within the coord contig",
        NULL},
    {"contig_end",
        (getter)Vertex_getcontigend, NULL,
        "end within the coord contig",
        NULL},
    {"saoffset",
        (getter)Vertex_getsaoffset, NULL,
        "offset in the suffix array",
        NULL},
    {"rcsaoffset",
        (getter)Vertex_getrcsaoffset, NULL,
        "offset of the reverse complement in the suffix array",
        NULL},
    {"edges_from",
        (getter)Vertex_getedgesfrom, NULL,
        "edges leaving this node",
        NULL},
    {"edges_to",
        (getter)Vertex_getedgesto, NULL,
        "edges ending in this node",
        NULL},
    {"attributes",
        (getter)Vertex_getattributes, NULL,
        "dictionary with additional attributes for this vertex",
        NULL},
    {NULL}  /* Sentinel */
};


static PyTypeObject VertexType = {
    PyObject_HEAD_INIT(NULL)
    0,                         /*ob_size*/
    "Vertex",             /*tp_name*/
    sizeof(Vertex),             /*tp_basicsize*/
    0,                         /*tp_itemsize*/
    (destructor)Vertex_dealloc, /*tp_dealloc*/
    0,                         /*tp_print*/
    0,                         /*tp_getattr*/
    0,                         /*tp_setattr*/
    0,                         /*tp_compare*/
    0,                         /*tp_repr*/
    0,                         /*tp_as_number*/
    0,                         /*tp_as_sequence*/
    0,                         /*tp_as_mapping*/
    0,                         /*tp_hash */
    0,                         /*tp_call*/
    0,                         /*tp_str*/
    0,                         /*tp_getattro*/
    0,                         /*tp_setattro*/
    0,                         /*tp_as_buffer*/
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE, /*tp_flags*/
    "Vertex objects",           /* tp_doc */
    0,		               /* tp_traverse */
    0,		               /* tp_clear */
    0,		               /* tp_richcompare */
    0,		               /* tp_weaklistoffset */
    0,		               /* tp_iter */
    0,		               /* tp_iternext */
    Vertex_methods,                         /* tp_methods */
    0,                         /* tp_members */
    Vertex_getseters,                         /* tp_getset */
    0,                         /* tp_base */
    0,                         /* tp_dict */
    0,                         /* tp_descr_get */
    0,                         /* tp_descr_set */
    0,                         /* tp_dictoffset */
    (initproc)Vertex_init,      /* tp_init */
    0,                         /* tp_alloc */
    Vertex_new,                 /* tp_new */
};





















//EDGE

static void
Edge_dealloc(Edge* self)
{
    self->ob_type->tp_free((PyObject*)self);
}

static int
Edge_init(Edge *self, PyObject *args, PyObject *kwds)
{
    if (! PyArg_ParseTuple(args, "OOi", &self->source, &self->target, &self->orientation)) {
        PyErr_SetString(PyExc_TypeError, "Cannot initialize Edge without source and target vertices");
        return -1;
    }
    
    return 0;
}

static int
Edge_init_c(Edge *self, Vertex* source, Vertex* target, int orientation)
{
    self->source=source;
    self->target=target;
	self->orientation=orientation;
    return 0;
}


static PyObject *
Edge_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
    Edge *self;
    
    self = (Edge *)type->tp_alloc(type, 0);
	
    if (self != NULL) {
        self->source = NULL;
        self->target = NULL;
    }

    return (PyObject *)self;
}

static PyMethodDef Edge_methods[] = {
    {NULL}  /* Sentinel */
};

static PyObject *
Edge_getsource(Edge *self, void *closure)
{
    Py_INCREF(self->source);
    return (PyObject*)self->source;
}
static PyObject *
Edge_gettarget(Edge *self, void *closure)
{
    Py_INCREF(self->target);
    return (PyObject*)self->target;
}

static PyObject *
Edge_getorientation(Edge *self, void *closure)
{
    return Py_BuildValue("i", self->orientation);
}

static PyGetSetDef Edge_getseters[] = {
    {"source",
        (getter)Edge_getsource, NULL,
        "source vertex",
        NULL},
    {"target",
        (getter)Edge_gettarget, NULL,
        "target vertex",
        NULL},
     {"orientation",
        (getter)Edge_getorientation, NULL,
        "orientation same=0, inny=1, outty=2",
        NULL},
    {NULL}  /* Sentinel */
};

static PyTypeObject EdgeType = {
    PyObject_HEAD_INIT(NULL)
    0,                         /*ob_size*/
    "Edge",             /*tp_name*/
    sizeof(Edge),             /*tp_basicsize*/
    0,                         /*tp_itemsize*/
    (destructor)Edge_dealloc, /*tp_dealloc*/
    0,                         /*tp_print*/
    0,                         /*tp_getattr*/
    0,                         /*tp_setattr*/
    0,                         /*tp_compare*/
    0,                         /*tp_repr*/
    0,                         /*tp_as_number*/
    0,                         /*tp_as_sequence*/
    0,                         /*tp_as_mapping*/
    0,                         /*tp_hash */
    0,                         /*tp_call*/
    0,                         /*tp_str*/
    0,                         /*tp_getattro*/
    0,                         /*tp_setattro*/
    0,                         /*tp_as_buffer*/
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE, /*tp_flags*/
    "Edge objects",           /* tp_doc */
    0,		               /* tp_traverse */
    0,		               /* tp_clear */
    0,		               /* tp_richcompare */
    0,		               /* tp_weaklistoffset */
    0,		               /* tp_iter */
    0,		               /* tp_iternext */
    Edge_methods,                         /* tp_methods */
    0,                         /* tp_members */
    Edge_getseters,                         /* tp_getset */
    0,                         /* tp_base */
    0,                         /* tp_dict */
    0,                         /* tp_descr_get */
    0,                         /* tp_descr_set */
    0,                         /* tp_dictoffset */
    (initproc)Edge_init,      /* tp_init */
    0,                         /* tp_alloc */
    Edge_new,                 /* tp_new */
};








//TODO: implement fast BFS iterator object

//typedef struct {
//    PyObject_HEAD
//    long int m;
//    long int i;
//    
//} bfs_Iter;
//
//
//PyObject* bfs_Iter_iter(PyObject *self)
//{
//    Py_INCREF(self);
//    return self;
//}
//
//PyObject* bfs_Iter_iternext(PyObject *self)
//{
//    bfs_Iter *p = (bfs_Iter *)self;
//    
//    if (p->i < p->m) {
//        PyObject *tmp = Py_BuildValue("l", p->i);
//        (p->i)++;
//        return tmp;
//    } else {
//        /* Raising of standard StopIteration exception with empty value. */
//        PyErr_SetNone(PyExc_StopIteration);
//        return NULL;
//    }
//}









//GRAPH

static void
Graph_dealloc(Graph* self)
{
    Py_DECREF(self->edges);
    Py_DECREF(self->vertices);
    self->ob_type->tp_free((PyObject*)self);
}

static int
Graph_init(Graph *self, PyObject *args, PyObject *kwds)
{
    //nothing to be initialized for now...
    return 0;
}

static PyObject *
Graph_New(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
    Graph *self;
    
    self = (Graph *)type->tp_alloc(type, 0);
    
    if (self != NULL) {
        
        self->vertices=PyDict_New();
        if (self->vertices == NULL)
        {
            Py_DECREF(self);
            return NULL;
        }
        
        self->edges=PySet_New(NULL);
        if (self->edges == NULL)
        {
            Py_DECREF(self);
            return NULL;
        }
		
		self->origins=PySet_New(NULL);
        if (self->origins == NULL)
        {
            Py_DECREF(self);
            return NULL;
        }
		
        self->vi=0;
        self->ei=0;
    }
    
    return (PyObject *)self;
}

static PyObject *
Graph_getvertices(Graph *self, void *closure)
{
    Py_INCREF(self->vertices);
    return self->vertices;
}

static PyObject *
Graph_getedges(Graph *self, void *closure)
{
    Py_INCREF(self->edges);
    return self->edges;
}

static PyObject *
Graph_getorigins(Graph *self, void *closure)
{
    Py_INCREF(self->origins);
    return self->origins;
}

static PyGetSetDef Graph_getseters[] = {
    {"vertices",
        (getter)Graph_getvertices, NULL,
        "vertices dict",
        NULL},
    {"edges",
        (getter)Graph_getedges, NULL,
        "set of all edges in the graph",
        NULL},
    {"origins",
        (getter)Graph_getorigins, NULL,
        "set containing all the names of the samples contained in the graph",
        NULL},
   {NULL}  /* Sentinel */
};

static PyObject * add_vertex(Graph *self, PyObject *args, PyObject *kwds) {

    idx_t contig_start, contig_end, rcsaoffset, saoffset;
	int input_origin, indexed;
	char * coord_origin;
	char * coord_contig;
    if (! PyArg_ParseTuple(args,"iiss" PYIDXTYPE PYIDXTYPE PYIDXTYPE PYIDXTYPE, &input_origin, &indexed, &coord_origin, &coord_contig, &contig_start, &contig_end, &saoffset, &rcsaoffset)) {
        PyErr_SetString(PyExc_TypeError, "Cannot parse parameters");
        return NULL;
    }
    
    Vertex *v = add_vertex_c(self, input_origin, indexed, coord_origin, coord_contig, contig_start, contig_end, saoffset, rcsaoffset);
    
    return (PyObject *)v;
}


static int edge_exists(Vertex* v1, Vertex* v2, int o){
	PyObject* iterto;
	PyObject* iterfrom;
	Edge* edge;
	
	int conn=0;
	
	iterfrom=PyObject_GetIter(v1->edges_from);
	iterto=PyObject_GetIter(v1->edges_to);
	if (iterfrom==NULL){
		PyErr_SetString(PyExc_TypeError, "Unable to iterate over set of v1 edges_from");
		return -1;
	}
	if (iterto==NULL){
		PyErr_SetString(PyExc_TypeError, "Unable to iterate over set of v1 edges_from");
		return -1;
	}
	
	
	if (o==0) { //normal orientation, source and target have meaning
		
		while ( (edge=(Edge*)PyIter_Next(iterfrom)) ) {
			//edges->source == v1
			if (edge->source->id==v1->id && edge->target->id==v2->id){
				conn=1;
				break;
			}
			
			Py_DECREF(edge);
		}
		
		
	} 
	
	if (o==1) { //trying to add an inny between source and target
		
		while ( (edge=(Edge*)PyIter_Next(iterfrom)) ) {
			
			if (edge->orientation==1) {
			
				if ( (edge->source->id==v1->id && edge->target->id==v2->id) || (edge->source->id==v2->id && edge->target->id==v1->id) ){
					conn=1;
					break;
				}

			}

			Py_DECREF(edge);
		}
		
	}

	if (o==2) { //trying to add an outty between source and target
		
		while ( (edge=(Edge*)PyIter_Next(iterto)) ) {
			
			if (edge->orientation==2) {
			
				if ( (edge->source->id==v1->id && edge->target->id==v2->id) || (edge->source->id==v2->id && edge->target->id==v1->id) ){
					conn=1;
					break;
				}

			}

			Py_DECREF(edge);
		}
		
	}


	Py_DECREF(iterfrom);
	Py_DECREF(iterto);


/*
	while ( (edge=(Edge*)PyIter_Next(iter)) ) {
		if (edge->orientation==0){
			if (edge->target->id==v2->id){
				conn=1;
				break;
			}
		} else {
			if (edge->source->id==v1->id){
				if (edge->target->id==v2->id) {
					conn=1;
					break;
				}
			} else {
				if (edge->source->id==v2->id) {
					conn=1;
					break;
				}
			}
		}
		
		Py_DECREF(edge);
	}
	Py_DECREF(iter);

	iter=PyObject_GetIter(v1->edges_to);
	if (conn!=1) {
		while ( (edge=(Edge*)PyIter_Next(iter)) ) {
			if (edge->orientation==0){
				if (edge->target->id==v2->id){
					conn=1;
					break;
				}
			} else {
				if (edge->source->id==v1->id){
					if (edge->target->id==v2->id) {
						conn=1;
						break;
					}
				} else {
					if (edge->source->id==v2->id) {
						conn=1;
						break;
					}
				}
			}
			
			Py_DECREF(edge);
		}
	}
	Py_DECREF(iter);*/

	return conn;
}

static PyObject * add_edge(Graph *self, PyObject *args, PyObject *kwds) {
    PyObject* source;
    PyObject* target;
   	int orientation;
	
    if (! PyArg_ParseTuple(args,"OOi", &source, &target, &orientation)) {
        PyErr_SetString(PyExc_TypeError, "Cannot parse parameters");
        return NULL;
    }
	
	if (edge_exists((Vertex*)source, (Vertex*)target, orientation)==0) {
		Edge *e = add_edge_c(self, (Vertex*)source, (Vertex*)target, orientation);
		return (PyObject*)e;
	}
	
	Py_INCREF(Py_None);

    return Py_None;
}

static PyObject * remove_vertex(Graph *self, PyObject *args, PyObject *kwds) {
    Vertex* vertex;
    Edge* edge;
    PyObject *iter;
    int edgesremoved=0;
    
    if (! PyArg_ParseTuple(args,"O", &vertex)) {
        PyErr_SetString(PyExc_TypeError, "Cannot parse parameters");
        return NULL;
    }
    
    Py_ssize_t s;
    s=PySet_Size(vertex->edges_from);
    
    if (s>0) {
        //remove all edges that leave this vertex
        iter=PyObject_GetIter(vertex->edges_from);
        
        if (iter==NULL){
            PyErr_SetString(PyExc_TypeError, "Unable to iterate over set of edges_from");
            return NULL;
        }
        
        while ( (edge=(Edge*)PyIter_Next(iter)) )
        {
            edgesremoved++;
			
			if (edge->orientation==0){
				if (PySet_Discard(edge->target->edges_to,(PyObject*)edge)!=1){
					PyErr_SetString(PyExc_TypeError, "Unable to remove edge from incoming edges of connected vertex, orientation 0");
					return NULL;
				}
			} else {
				//edge-orientation then has to be 1
				assert(edge->orientation==1);

				//vertex can be both source and target!
				if (vertex->id == edge->source->id){
					if (PySet_Discard(edge->target->edges_from,(PyObject*)edge)!=1){
						PyErr_SetString(PyExc_TypeError, "Unable to remove edge from incoming edges of connected vertex, orientation 1");
						return NULL;
					}
				}
				else {
					if (PySet_Discard(edge->source->edges_from,(PyObject*)edge)!=1){
					PyErr_SetString(PyExc_TypeError, "Unable to remove edge from incoming edges of connected vertex, orientation 1");
					return NULL;
					}
				}
			}

            if (PySet_Discard(self->edges,(PyObject*)edge)!=1){
                PyErr_SetString(PyExc_TypeError, "Unable to remove edge from graph defined set of edges");
                return NULL;
            }

			Py_DECREF(edge);
        }
        
        Py_DECREF(iter);
    }

    s=PySet_Size(vertex->edges_to);
    
    if (s>0) {
        //remove all edges that point towards this vertex
        iter=PyObject_GetIter(vertex->edges_to);
    
        if (iter==NULL){
            PyErr_SetString(PyExc_TypeError, "Unable to iterate over set of edges_to");
            return NULL;
        }
    
        while ( (edge=(Edge*)PyIter_Next(iter)) )
        {
            edgesremoved++;
			
			if (edge->orientation==0){
				if (PySet_Discard(edge->source->edges_from,(PyObject*)edge)!=1){
					PyErr_SetString(PyExc_TypeError, "Unable to remove edge from outgoing edges of connected vertex, orientation 0");
					return NULL;
				}
			} else {
				//edge-orientation then has to be 2

				if (vertex->id == edge->target->id){
					if (PySet_Discard(edge->source->edges_to,(PyObject*)edge)!=1){
						PyErr_SetString(PyExc_TypeError, "Unable to remove edge from outgoing edges of connected vertex, orientation 2");
						return NULL;
					}
				} else {
					if (PySet_Discard(edge->target->edges_to,(PyObject*)edge)!=1){
						PyErr_SetString(PyExc_TypeError, "Unable to remove edge from outgoing edges of connected vertex, orientation 2");
						return NULL;
					}
				}
			}

            if (PySet_Discard(self->edges,(PyObject*)edge)!=1){
                PyErr_SetString(PyExc_TypeError, "Unable to remove edge from graph defined set of edges");
                return NULL;
            }
			
			Py_DECREF(edge);
        }
        
        Py_DECREF(iter);
    }
    
    PyObject* vid=Py_BuildValue("i",vertex->id);
    
    if (PyDict_DelItem(self->vertices, vid)!=0){
        PyErr_SetString(PyExc_TypeError, "Unable to remove vertex from graph dictionary");
        return NULL;
    }
    
    Py_DECREF(vid);
    
    return Py_BuildValue("i",edgesremoved);
}


static PyMethodDef Graph_methods[] = {
    {"add_vertex", (PyCFunction) add_vertex, METH_VARARGS, "Adds a vertex to the graph."},
    {"add_edge", (PyCFunction) add_edge, METH_VARARGS, "Adds an edge to the graph."},
    {"remove_vertex", (PyCFunction) remove_vertex, METH_VARARGS, "Removes a vertex from the graph, returns the number of removed edges."},
    //{"bfs", (PyCFunction) bfs, METH_VARARGS, "Returns an iterator that performs a Breadth First search."},
    {NULL}  /* Sentinel */
};


static PyTypeObject GraphType = {
    PyObject_HEAD_INIT(NULL)
    0,                         /*ob_size*/
    "Graph",             /*tp_name*/
    sizeof(Graph),             /*tp_basicsize*/
    0,                         /*tp_itemsize*/
    (destructor)Graph_dealloc, /*tp_dealloc*/
    0,                         /*tp_print*/
    0,                         /*tp_getattr*/
    0,                         /*tp_setattr*/
    0,                         /*tp_compare*/
    0,                         /*tp_repr*/
    0,                         /*tp_as_number*/
    0,                         /*tp_as_sequence*/
    0,                         /*tp_as_mapping*/
    0,                         /*tp_hash */
    0,                         /*tp_call*/
    0,                         /*tp_str*/
    0,                         /*tp_getattro*/
    0,                         /*tp_setattro*/
    0,                         /*tp_as_buffer*/
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE, /*tp_flags*/
    "Graph objects",           /* tp_doc */
    0,		               /* tp_traverse */
    0,		               /* tp_clear */
    0,		               /* tp_richcompare */
    0,		               /* tp_weaklistoffset */
    0,		               /* tp_iter */
    0,		               /* tp_iternext */
    Graph_methods,             /* tp_methods */
    0,                         /* tp_members */
    Graph_getseters,                         /* tp_getset */
    0,                         /* tp_base */
    0,                         /* tp_dict */
    0,                         /* tp_descr_get */
    0,                         /* tp_descr_set */
    0,                         /* tp_dictoffset */
    (initproc)Graph_init,      /* tp_init */
    0,                         /* tp_alloc */
    Graph_New,                 /* tp_new */
};




int compute_lcp(const unsigned char *T, const idx_t *SA, const idx_t *SAi, unsigned int *LCP, idx_t n) {
	idx_t h=0, i, j, k;
    
    for (i = 0; i < n; i++) {
        k = SAi[i];
        if (k == 0) {
            LCP[k] = 0; //not -1, because we need unsigned int
        } else {
            j = SA[k-1];
            while ((i - h < n) && (j + h < n) && (T[i+h] == T[j+h]) && (T[i+h]!='$' && T[j+h]!='$') && (T[i+h]!='N' && T[j+h]!='N') ) { ++h; } //stop comparing when a sentinel or N is encountered, so we dont find matches that span them
            LCP[k] = h;
        }
        if (h > 0) --h;
    }
    
    return 0;
}


static int
GSA_init(GSA *self, PyObject *args, PyObject *kwds)
{
    if (args==NULL) {
        return 0;
    }
    
    //static char* kwlist[]= {"input1","input2","userc", NULL};
    
	char * input1;
	char * input2;

	memset(comp,'N',128);
	comp['a']='t';
	comp['c']='g';
	comp['g']='c';
	comp['t']='a';
	comp['A']='T';
	comp['C']='G';
	comp['G']='C';
	comp['T']='A';

	comp['Y']='R';
	comp['R']='Y';
	comp['K']='M';
	comp['M']='K';

	comp['S']='S';
	comp['W']='W';

	comp['B']='V';
	comp['V']='B';
	comp['D']='H';
	comp['H']='D';
	comp['X']='X';
	comp['-']='-';


	if (PyTuple_Size(args)==3){ //index two graphs with the supplied parameters

		if (! PyArg_ParseTuple(args, "ssi", &input1, &input2, &self->rcindex)) {
			PyErr_SetString(GSAError, "Failed to parse input parameters");
			return -1;
		}

		printf("Start indexing: %s and %s\n",input1,input2);
		if (buildGSA(self, input1, input2)!=0) {
			PyErr_SetString(GSAError, "Failed to build GSA");
			return -1;
		}

		self->root=(PyObject*)self;
	}
	
	if (PyTuple_Size(args)==1){ //just load the graph of one file without indexing anything

		if (!PyArg_ParseTuple(args, "s", &input1)) {
			PyErr_SetString(GSAError, "Failed to parse input parameters");
			return -1;
		}

		idx_t tmp=0;
		
		self->graph = (Graph *)Graph_New(&GraphType, NULL, NULL);
		self->rcindex=1; //by default create T with reverse complements
		
		if (!(tmp=loadFile(self, input1, 0, tmp, self->rcindex))) {
			PyErr_SetString(GSAError, "Failed to load graph.");
			return -1;
		}
		
		self->n=tmp;
		self->orgn=tmp;

		self->root=(PyObject*)self;
	}

	return 0;
}


static PyObject* get_mums(GSA * self, PyObject* args, PyObject *kwds)
{
   unsigned int threshold;
   if (!PyArg_ParseTuple(args, "I", &threshold)) {
	  PyErr_SetString(GSAError, "Failed to parse input parameters");
	  return NULL;
   }

   if (self->LCP==0) {
	  PyErr_SetString(PyExc_TypeError, "No LCP with this index.");
	  return NULL;
   }

   PyObject *mums=PyList_New(0);

   idx_t i=0;
   idx_t aStart, bStart;

   Vertex* v1=NULL;
   Vertex* v2=NULL;

   unsigned int lb, la, l;
   int rcmatch=0;

   for (i=1; i<self->n; i++) {

	  l=self->LCP[i];

	  if (l>threshold){ //is it larger than the specified threshold

		 if ( self->SA[i] < self->sep ) {
			if ( self->SA[i-1] > self->sep ) {
			   aStart=self->SA[i];
			   bStart=self->SA[i-1];
			} else {
			   continue; //same set of sequences
			}
		 } else {
			if ( self->SA[i-1] < self->sep ) {
			   aStart=self->SA[i-1];
			   bStart=self->SA[i];
			} else {
			   continue; //same set of sequences
			}
		 }

		 if (aStart>0 && bStart>0){ //if not it has to be maximal!
			if (!((self->T[aStart-1]!=self->T[bStart-1]) || (self->T[aStart-1]=='$') || (islower(self->T[aStart-1])) )) { //is it maximal?
			   continue;
			}
		 }

		 if (i==self->n-1) { //is it the last value in the array, then only check predecessor
			lb=self->LCP[i-1];
			la=0;
		 } else {
			lb=self->LCP[i-1];
			la=self->LCP[i+1];
		 }

		 if (!((lb<l) && (la<l))) { //is it unique?
			continue;
		 }

		 v1=self->TIiv[aStart];
		 v2=self->TIiv[bStart];
		 rcmatch=0;
		 if (self->rcindex==1){
			if (v2->rcsaoffset!=v2->saoffset && bStart>=v2->rcsaoffset) {
			   bStart=v2->saoffset+((v2->contig_end-v2->contig_start)-l-(bStart-v2->rcsaoffset));
			   rcmatch=1;
			}
		 }

		 PyObject *v=Py_BuildValue("("PYIDXTYPE","PYIDXTYPE",I,O,O,i)", aStart, bStart, l, v1, v2, rcmatch);
		 PyList_Append(mums, v);
		 Py_DECREF(v);
	  }
   }

   return mums;
}



//only works when aligning two vertices within a bubble, penalty within a complexer bubble structure (multi-alignment?) is not straightforward to calculate
static PyObject* get_scored_mum(GSA * self, PyObject* args, PyObject *kwds)
{
   if (self->LCP==0) {
	  PyErr_SetString(PyExc_TypeError, "No LCP with this index.");
	  return NULL;
   }
   
   idx_t i=0;
   idx_t aStart=0, bStart=0, maxaStart=0, maxbStart=0;

   Vertex* v1=NULL;
   Vertex* v2=NULL;
   Vertex* maxv1=NULL;
   Vertex* maxv2=NULL;
   
   unsigned int lb, la, l, maxl=0;
   int rcmatch=0, maxrcmatch=0, firstMUM=1;
   
   idx_t seq1start=0, seq2start=0, seq1end=0, seq2end=0, pg1=0, pg2=0, sg1=0, sg2=0, trailing_gap=0, leading_gap=0, gappenalty=0, score=0, maxscore=0, besttrailinggap=0, bestleadinggap=0;

   int leftaligned=0, rightaligned=0, gap_open=-2, gap_extend=-1;
   static char *kwlist[] = {"left_aligned", "right_aligned", "gap_open", "gap_extend", NULL};

   if (!PyArg_ParseTupleAndKeywords(args, kwds, "|iiii", kwlist, &leftaligned, &rightaligned, &gap_open, &gap_extend))
	  return NULL;
   
   for (i=1; i<self->n; i++) {
	  l=self->LCP[i];
	  
	  if (l==0) {
		 continue;
	  }
	  
	  if ( self->SA[i] < self->sep ) {
		 if ( self->SA[i-1] > self->sep ) {
			aStart=self->SA[i];
			bStart=self->SA[i-1];
		 } else {
			continue; //same set of sequences
		 }
	  } else {
		 if ( self->SA[i-1] < self->sep ) {
			aStart=self->SA[i-1];
			bStart=self->SA[i];
		 } else {
			continue; //same set of sequences
		 }
	  }
	 
	  if (aStart>0 && bStart>0){ //if not it has to be maximal!
		 if (!((self->T[aStart-1]!=self->T[bStart-1]) || (self->T[aStart-1]=='$') || (islower(self->T[aStart-1])) )) { //is it maximal?
			continue;
		 }
	  }

	  if (i==self->n-1) { //is it the last value in the array, then only check predecessor
		 lb=self->LCP[i-1];
		 la=0;
	  } else {
		 lb=self->LCP[i-1];
		 la=self->LCP[i+1];
	  }

	  if (!((lb<l) && (la<l))) { //is it unique?
		 continue;
	  }
	  
	  //ok, so it's a maximal unique match, how well does it score?
	  v1=self->TIiv[aStart];
	  v2=self->TIiv[bStart];
	  rcmatch=0;
	  if (self->rcindex==1){
		 if (v2->rcsaoffset!=v2->saoffset && bStart>=v2->rcsaoffset) {
			bStart=v2->saoffset+((v2->contig_end-v2->contig_start)-l-(bStart-v2->rcsaoffset));
			rcmatch=1;
		 }
	  }
	  
	  if (firstMUM==1){ //constants, shouldn't change
		 seq1start=v1->saoffset;
		 seq2start=v2->saoffset;
		 seq1end=v1->saoffset+(v1->contig_end-v1->contig_start);
		 seq2end=v2->saoffset+(v2->contig_end-v2->contig_start);
	  }
	  
	  if (leftaligned){ //in case of a tip, we can only calculate penalty for either prefix or suffix
		 pg1=aStart-seq1start;//prefix gap caused by aligning this MUM on sequence 1
		 pg2=bStart-seq2start;//prefix gap caused by aligning this MUM on sequence 2
	  }
	  
	  if (rightaligned){ //in case of a tip, we can only calculate penalty for either prefix or suffix
		 sg1=seq1end-aStart-l;//suffix gap caused by aligning this MUM on sequence 1
		 sg2=seq2end-bStart-l;//suffix gap caused by aligning this MUM on sequence 2
	  }

	  leading_gap=abs(pg1-pg2)*gap_extend;
	  trailing_gap=abs(sg1-sg2)*gap_extend;
	  
	  if (leading_gap!=0){
		 leading_gap=leading_gap+gap_open;
	  }
	  if (trailing_gap!=0){
		 trailing_gap=trailing_gap+gap_open;
	  }
	  
	  gappenalty=leading_gap+trailing_gap; //TODO: for now match score for all nucleotides is 1, could add posibility of a scoring matrix
	  score=l+gappenalty;
	  
   	  if (score>maxscore || firstMUM==1){ //is it larger than the currently found max or first entry of SA
		 maxl=l;
		 maxscore=score;
		 maxaStart=aStart;
		 maxbStart=bStart;
		 maxv1=v1;
		 maxv2=v2;
		 maxrcmatch=rcmatch;
		 besttrailinggap=trailing_gap;
		 bestleadinggap=leading_gap;
		 firstMUM=0;
	  }
   }

   if (maxl>0){
	  return Py_BuildValue("("PYIDXTYPE","PYIDXTYPE",I,O,O,i,i,"PYIDXTYPE","PYIDXTYPE")", maxaStart, maxbStart, maxl, maxv1, maxv2, maxrcmatch, maxscore, bestleadinggap, besttrailinggap);
   } else {
	  Py_INCREF(Py_None);
	  return Py_None;
   }
}





static PyObject* get_mum(GSA * self, PyObject* args, PyObject *kwds)
{
   if (self->LCP==0) {
	  PyErr_SetString(PyExc_TypeError, "No LCP with this index.");
	  return NULL;
   }

   idx_t i=0;
   idx_t aStart=0, bStart=0, maxaStart=0, maxbStart=0;

   Vertex* v1=NULL;
   Vertex* v2=NULL;
   Vertex* maxv1=NULL;
   Vertex* maxv2=NULL;

   unsigned int lb, la, l, maxl=0;
   int rcmatch=0, maxrcmatch=0;

   for (i=1; i<self->n; i++) {

	  l=self->LCP[i];
	  
	  if (l>maxl){ //is it larger than the currently found max
		 if ( self->SA[i] < self->sep ) {
			if ( self->SA[i-1] > self->sep ) {
			   aStart=self->SA[i];
			   bStart=self->SA[i-1];
			} else {
			   continue; //same set of sequences
			}
		 } else {
			if ( self->SA[i-1] < self->sep ) {
			   aStart=self->SA[i-1];
			   bStart=self->SA[i];
			} else {
			   continue; //same set of sequences
			}
		 }
		 
		 if (aStart>0 && bStart>0){ //if not it has to be maximal!
			if (!((self->T[aStart-1]!=self->T[bStart-1]) || (self->T[aStart-1]=='$') || (islower(self->T[aStart-1])) )) { //is it maximal?
			   continue;
			}
		 }
		 
		 if (i==self->n-1) { //is it the last value in the array, then only check predecessor
			lb=self->LCP[i-1];
			la=0;
		 } else {
			lb=self->LCP[i-1];
			la=self->LCP[i+1];
		 }

		 if ((lb<l) && (la<l)) { //is it unique?
			v1=self->TIiv[aStart];
			v2=self->TIiv[bStart];
			rcmatch=0;
			if (self->rcindex==1){
			   if (v2->rcsaoffset!=v2->saoffset && bStart>=v2->rcsaoffset) {
				  bStart=v2->saoffset+((v2->contig_end-v2->contig_start)-l-(bStart-v2->rcsaoffset));
				  rcmatch=1;
			   }
			}
			maxl=l;
			maxaStart=aStart;
			maxbStart=bStart;
			maxv1=v1;
			maxv2=v2;
			maxrcmatch=rcmatch;
		 }
	  }
   }
   
   if (maxl>0){
	  return Py_BuildValue("("PYIDXTYPE","PYIDXTYPE",I,O,O,i)", maxaStart, maxbStart, maxl, maxv1, maxv2, maxrcmatch);
   } else {
	  Py_INCREF(Py_None);
	  return Py_None;
   }
}



static void
GSA_dealloc(GSA* self)
{
    if (self->SA!=0) {
        free(self->SA);
    }

    if (self->LCP!=0) {
        free(self->LCP);
    }

    //only free memory for these datasets when the initial index is deallocated
    //consider that sub indices keep reference to it!
    if (self->n == self->orgn) {
        
        if (self->T!=0) {
            free(self->T);
        }
        
        if (self->SAi!=0) {
            free(self->SAi);
        }
        
        if (self->TIiv!=0) {
            free(self->TIiv);
        }
        
        Py_DECREF(self->graph);
    } else {
        Py_DECREF(self->root); //should decref the root index
    }

    self->ob_type->tp_free((PyObject*)self);
}

static PyObject *
GSA_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
    GSA *self;
    
    self = (GSA *)type->tp_alloc(type, 0);
    
    if (self!=NULL) {
        self->SA=0;
        self->SAi=0;
        self->LCP=0;
        self->TIiv=0;
        self->T=0;
        self->root=0;
        self->n=0;
        self->orgn=0;
        self->sep=0;
    }
    
    return (PyObject *)self;
}

static PyMethodDef gsa_methods[] = {
    {"split", (PyCFunction) split, METH_VARARGS, "Split the index into two new GSA indices."},
    {"update", (PyCFunction) update, METH_VARARGS, "Updates the index such that the order of suffix complies after aligning."},
    {"updateVertex", (PyCFunction) updateVertex, METH_VARARGS, "Updates the index such that the order of suffix complies after extracting the specified vertex from the index."},
    {"extract", (PyCFunction) extract, METH_VARARGS, "Extract an index for the specified texts."},
    {"relabel", (PyCFunction) relabel, METH_VARARGS, "Relabel the references in TIiv for the specified vertex."},
    {"get_mum", (PyCFunction) get_mum, METH_VARARGS, "Returns the index of the largest maximal unique match between the two sets of sequences."},
    {"get_scored_mum", (PyCFunction) get_scored_mum, METH_VARARGS|METH_KEYWORDS, "Returns the index of the best scoring maximal unique match between the two sets of sequences, scoring considers gaps resulting from aligning this MUM."},
    {"get_mums", (PyCFunction) get_mums, METH_VARARGS, "Returns the indices of all maximal unique matches that are bigger then the predefined threshold."},
    //{"get_mums_fast", (PyCFunction) get_mums_fast, METH_VARARGS, "Returns the indices of all maximal unique matches that are bigger then the predefined threshold."},
    {NULL}  /* Sentinel */
};


static PyObject *
GSA_getn(GSA *self, void *closure)
{
    return Py_BuildValue(PYIDXTYPE, self->n);
}
static PyObject *
GSA_getorgn(GSA *self, void *closure)
{
    return Py_BuildValue(PYIDXTYPE, self->orgn);
}
static PyObject *
GSA_getsep(GSA *self, void *closure)
{
    return Py_BuildValue(PYIDXTYPE, self->sep);
}


static PyObject *
GSA_getrcindex(GSA *self, void *closure)
{
    return Py_BuildValue("i", self->rcindex);
}

static PyObject *
GSA_getTIiv(GSA *self, void *closure)
{
    if (self->TIiv==0) {
        PyErr_SetString(PyExc_TypeError, "No TIiv with this index.");
        return NULL;
    }
    
    PyObject *lst = PyList_New(self->orgn);
    
    if (!lst)
        return NULL;
    
    idx_t i;
    for (i = 0; i < self->orgn; i++) {
		if (self->TIiv[i]==NULL){
			Py_INCREF(Py_None);
			PyList_SET_ITEM(lst, i, Py_None);
			continue;
		}
        PyObject *num = Py_BuildValue("O", self->TIiv[i]);
        if (!num) {
            Py_DECREF(lst);
            return NULL;
        }
        PyList_SET_ITEM(lst, i, num);
    }
    return lst;
}

static PyObject *
GSA_getSA(GSA *self, void *closure)
{
    if (self->SA==0) {
        PyErr_SetString(PyExc_TypeError, "No SA with this index.");
        return NULL;
    }
    
    PyObject *lst = PyList_New(self->n);
    
    if (!lst)
        return NULL;
    
    idx_t i;
    for (i = 0; i < self->n; i++) {
        PyObject *num = Py_BuildValue(PYIDXTYPE, self->SA[i]);
        if (!num) {
            Py_DECREF(lst);
            return NULL;
        }
        PyList_SET_ITEM(lst, i, num);
    }
    
    return lst;
}

static PyObject *
GSA_getSAi(GSA *self, void *closure)
{
    if (self->SAi==0) {
        PyErr_SetString(PyExc_TypeError, "No SAi with this index.");
        return NULL;
    }
    
    PyObject *lst = PyList_New(self->orgn);
    
    if (!lst)
        return NULL;
    
    idx_t i;
    for (i = 0; i < self->orgn; i++) {
        PyObject *num = Py_BuildValue(PYIDXTYPE, self->SAi[i]);
        if (!num) {
            Py_DECREF(lst);
            return NULL;
        }
        PyList_SET_ITEM(lst, i, num);
    }
    
    return lst;
}

static PyObject *
GSA_getLCP(GSA *self, void *closure)
{
    if (self->LCP==0) {
        PyErr_SetString(PyExc_TypeError, "No LCP with this index.");
        return NULL;
    }
    
    PyObject *lst = PyList_New(self->n);
    
    if (!lst)
        return NULL;
    
    idx_t i;
    for (i = 0; i < self->n; i++) {
        PyObject *num = Py_BuildValue("i", self->LCP[i]);
        if (!num) {
            Py_DECREF(lst);
            return NULL;
        }
        PyList_SET_ITEM(lst, i, num);
    }
    
    return lst;
}

static PyObject *
GSA_getT(GSA *self, void *closure)
{
    return Py_BuildValue("s", self->T);
}

static Graph *
GSA_getGraph(GSA *self, void *closure)
{
    if (self->graph==NULL) {
        PyErr_SetString(PyExc_TypeError, "No associated graph object.");
        return NULL;
    }
    
    Py_INCREF(self->graph);
    return self->graph;
}

static int
GSA_setGraph(GSA *self, PyObject *g, void *closure)
{
	if (g == NULL) {
		Py_DECREF(self->graph);
		self->graph=NULL;
    	return -1;
	}
  
	if (self->graph!=NULL) {
		Py_DECREF(self->graph);
	} else {
		Py_INCREF(g);
		self->graph = (Graph *)g;
	}

	return 0;
}

static PyGetSetDef gsa_getseters[] = {
    {"n",
        (getter)GSA_getn, NULL,
        "Number of characters in the index.",
        NULL},
    {"orgn",
        (getter)GSA_getorgn, NULL,
        "Number of characters in the initial index.",
        NULL},
    {"sep",
        (getter)GSA_getsep, NULL,
        "Returns the index of the sentinel that seperates both sets of input sequences.",
        NULL},
    {"rcindex",
        (getter)GSA_getrcindex, NULL,
        "Whether the reverse complements are indexed as well.",
        NULL},
    {"SA",
        (getter)GSA_getSA, NULL,
        "The suffix array of the concatenation of input texts.",
        NULL},
    {"SAi",
        (getter)GSA_getSAi, NULL,
        "The inverse of the suffix array.",
        NULL},
    {"T",
        (getter)GSA_getT, NULL,
        "The concatenation of the input texts.",
        NULL},
    {"LCP",
        (getter)GSA_getLCP, NULL,
        "List specifying the length of the common prefix of consecutive values in the LCP array.",
        NULL},
   {"TIiv",
        (getter)GSA_getTIiv, NULL,
        "List specifying the vertices for every suffix.",
        NULL},
   {"graph",
        (getter)GSA_getGraph, (setter)GSA_setGraph,
        "The corresponding graph object.",
        NULL},
    {NULL}  /* Sentinel */
};

static PyTypeObject GSAType = {
    PyObject_HEAD_INIT(NULL)
    0,                         /*ob_size*/
    ""PYSONAME,                     /*tp_name*/
    sizeof(GSA),             /*tp_basicsize*/
    0,                         /*tp_itemsize*/
    (destructor)GSA_dealloc, /*tp_dealloc*/
    0,                         /*tp_print*/
    0,                         /*tp_getattr*/
    0,                         /*tp_setattr*/
    0,                         /*tp_compare*/
    0,                         /*tp_repr*/
    0,                         /*tp_as_number*/
    0,                         /*tp_as_sequence*/
    0,                         /*tp_as_mapping*/
    0,                         /*tp_hash */
    0,                         /*tp_call*/
    0,                         /*tp_str*/
    0,                         /*tp_getattro*/
    0,                         /*tp_setattro*/
    0,                         /*tp_as_buffer*/
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE, /*tp_flags*/
    "Enhanced Generalized Suffix Array",           /* tp_doc */
    0,		               /* tp_traverse */
    0,		               /* tp_clear */
    0,		               /* tp_richcompare */
    0,		               /* tp_weaklistoffset */
    0,		               /* tp_iter */
    0,		               /* tp_iternext */
    gsa_methods,             /* tp_methods */
    0,             /* tp_members */
    gsa_getseters,                         /* tp_getset */
    0,                         /* tp_base */
    0,                         /* tp_dict */
    0,                         /* tp_descr_get */
    0,                         /* tp_descr_set */
    0,                         /* tp_dictoffset */
    (initproc)GSA_init,      /* tp_init */
    0,                         /* tp_alloc */
    GSA_new,                 /* tp_new */
};




/*

 Fasta header contains vertex information and edges. Parse it and return the created Vertex.

 */
static Vertex* parseHeader(GSA* self, int input_origin, int t, char *header, PyObject *mapping, unsigned int** edges_from, unsigned int** edges_to, uint8_t** orientations, unsigned int * nrofedges, unsigned int * esize, int useRC){

	Vertex* vid=NULL;

	char *cur=header;
	char *next;

	next=strchr(cur,'|');
	next[0]='\0';
	//idx_t oldid=strtol(cur,NULL,10);
	unsigned int oldid=atoi(cur);

	int edgebuf=10000; //buffer size

	cur=next+1;
	next=strchr(cur,'|');
	next[0]='\0';
	char *el=cur;

	char *es=el;
	char *ee=strchr(es,';');
	while (ee!=NULL){
		ee[0]='\0';
		char *p=strchr(es,',');
		uint8_t o=atoi(p+1);
		p[0]='\0';
		unsigned int target=atoi(es);//strtol(es,NULL,10);

		if (o==0 || (o==2 && target>oldid)){
			if (*nrofedges==*esize){
				*esize=*esize+edgebuf;
				unsigned int* tmpul;
				uint8_t * tmpb;
				tmpul=realloc(*edges_from,*esize*sizeof(unsigned int));
				*edges_from=tmpul;

				tmpul=realloc(*edges_to,*esize*sizeof(unsigned int));
				*edges_to=tmpul;
				
				tmpb=realloc(*orientations,*esize*sizeof(uint8_t));
				*orientations=tmpb;
			}
			(*edges_from)[(*nrofedges)]=target;
			(*edges_to)[(*nrofedges)]=oldid;
			(*orientations)[(*nrofedges)]=o;
			(*nrofedges)++;
		}
		
		es=ee+1;
		ee=strchr(es,';');
	}
	
	cur=next+1;
	next=strchr(cur,'|');
	next[0]='\0';
	char * er=cur;

	es=er;
	ee=strchr(es,';');
	while (ee!=NULL){
		ee[0]='\0';
		char *p=strchr(es,',');
		uint8_t o=atoi(p+1);
		p[0]='\0';
		unsigned int target=atoi(es);//strtol(es,NULL,10);

		if (o==1 && target>oldid) {
			if (*nrofedges==*esize){
				*esize=*esize+edgebuf;

				unsigned int* tmpul;
				uint8_t * tmpb;

				tmpul=realloc(*edges_from,*esize*sizeof(unsigned int));
				*edges_from=tmpul;

				tmpul=realloc(*edges_to,*esize*sizeof(unsigned int));
				*edges_to=tmpul;
				
				tmpb=realloc(*orientations,*esize*sizeof(uint8_t));
				*orientations=tmpb;
			}

			(*edges_from)[(*nrofedges)]=oldid;
			(*edges_to)[(*nrofedges)]=target;
			(*orientations)[(*nrofedges)]=o;
			(*nrofedges)++;
		}

		es=ee+1;
		ee=strchr(es,';');
	}

	cur=next+1;
	next=strchr(cur,'|');
	next[0]='\0';
	char * origin=cur;

	cur=next+1;
	next=strchr(cur,'|');
	next[0]='\0';
	char * corigin=cur;

	cur=next+1;
	next=strchr(cur,'|');
	next[0]='\0';
	char * coordorigin=cur;

	cur=next+1;
	next=strchr(cur,'|');
	next[0]='\0';
	char * coordcontig=cur;

	cur=next+1;
	next=strchr(cur,'|');
	next[0]='\0';
	//int cstart=atoi(cur);
	idx_t cstart=strtol(cur,NULL,10);

	cur=next+1;
	//int cend=atoi(cur);
	idx_t cend=strtol(cur,NULL,10);

	vid=add_vertex_c(self->graph, input_origin, useRC==1 ? 2 : 0, coordorigin, coordcontig, cstart, cend, t, useRC==1 ? t+(cend-cstart)+1 : t);

	PyObject* o;
	char *token = strtok(origin,";");
	while(token!=NULL)
	{
		o=Py_BuildValue("s",token);
		if (PySet_Add(vid->originid_set,o)==-1){
			fprintf(stderr,"Failed to add origin of node to set");
			return NULL;
		};
		Py_DECREF(o);
		token = strtok(NULL, ";");
	}

	token = strtok(corigin,";");
	while(token!=NULL)
	{
		o=Py_BuildValue("s",token);
		if (PySet_Add(vid->contigid_set,o)==-1){
			fprintf(stderr,"Failed to add origin of node to set");
			return NULL;
		};
		Py_DECREF(o);
		token = strtok(NULL, ";");
	}

	//update the mapping table so we can relate old nodes to new nodes
	//PyObject *id = PyInt_FromLong(oldid);
	PyObject *id = Py_BuildValue("i",oldid);
	//PyObject *id = PyInt_FromInt(oldid);
	if (id == NULL){
		Py_DECREF(id);
		fprintf(stderr,"Failed to read id %u.\n",oldid);
		return NULL;
	}
	
	if (PyDict_SetItem(mapping, id, (PyObject *)vid) < 0) {
		Py_DECREF(id);
		fprintf(stderr,"Failed to add to mapping for id %u.\n",oldid);
		return NULL;
	}
	Py_DECREF(id);
	
	return vid;
}


static idx_t 
loadFile(GSA *self, char * input1, int inputid, idx_t f, int useRC){

	int isgraph1=0;
	if (endsWith(input1, GRAPHEXT) || endsWith(input1, GRAPHEXTGZ))	{
		isgraph1=1;
	}
	
	gzFile fp1;
	kseq_t *seq1;
	int l;

	//read sequences from input1
	fp1 = gzopen(input1, "r");
	if (fp1==NULL) { return -1; }
	seq1 = kseq_init(fp1);

	unsigned int nrofedges=0;
	PyObject *mapping1=PyDict_New(); //dict used to map between old and new ids for edges

	const unsigned int edgebuf=10000; //TODO: use MACRO constant

	const idx_t Tbufsize=10000; //TODO: use MACRO constant
	idx_t Tsize=f+Tbufsize;

	unsigned int esize=edgebuf;
	unsigned int *edges_from=malloc(sizeof(int)*esize);
	unsigned int *edges_to=malloc(sizeof(int)*esize);
	uint8_t *orientations=malloc(sizeof(uint8_t)*esize);
	
	self->T=realloc(self->T, Tsize*sizeof(char));
	self->TIiv=realloc(self->TIiv, Tsize*sizeof(PyObject*));

	idx_t t=f;
	idx_t i=0;

	char *tmp;
	tmp=strrchr(input1,'/');
	input1= tmp ? tmp+1 : input1; //hack to get rid of the path information in the graph

	//if gfasta, first line contains all the samples that are contained in the graph
	if (isgraph1!=0){
		l = kseq_read(seq1);
		
		char *o=seq1->name.s;
		char *oe=strchr(o,'|');
		char *oide=NULL;
		
		while (oe!=NULL) {
			oe[0]='\0';
			oide=strchr(o,';');
			oide[0]='\0';
			
			PyObject* org=Py_BuildValue("s",oide+1);
			PySet_Add(self->graph->origins, org);
			Py_DECREF(org);

			o=oe+1;
			oe=strchr(o,'|');
		}
	} else {
		PyObject* org=Py_BuildValue("s",input1);
		PySet_Add(self->graph->origins, org);
		Py_DECREF(org);
	}
	
	while ((l = kseq_read(seq1)) >= 0) { 
		//create vertex for this sequence
		Vertex* vid;

		if (isgraph1!=0){
			vid=parseHeader(self, inputid, t, seq1->name.s, mapping1, &edges_from, &edges_to, &orientations, &nrofedges, &esize, useRC);
		} else {
			if (useRC==1){
				vid=add_vertex_c(self->graph, inputid, 2, input1, seq1->name.s, 0, seq1->seq.l, t, t+seq1->seq.l+1);
			} else {
				vid=add_vertex_c(self->graph, inputid, 0, input1, seq1->name.s, 0, seq1->seq.l, t, 0);
			}

			PyObject*org=Py_BuildValue("s",input1);
			PySet_Add(vid->originid_set,org);
			Py_DECREF(org);

			//add i to vid->contigid_set
			PyObject*con=Py_BuildValue("s",seq1->name.s);
			PySet_Add(vid->contigid_set,con);
			Py_DECREF(con);
		}

		t=f+seq1->seq.l+1;
		
		if (t>=Tsize){
			while (t>Tsize){Tsize=Tsize+Tbufsize;}
			self->T=realloc(self->T, Tsize*sizeof(char));
			self->TIiv=realloc(self->TIiv, Tsize*sizeof(PyObject*));
		}
		
		for (i=0;i<seq1->seq.l;i++){
			seq1->seq.s[i]=(char)toupper(seq1->seq.s[i]);
		}
		
		if (self->T!=NULL){
		   memcpy(&self->T[f], seq1->seq.s, seq1->seq.l);
		   memcpy(&self->T[f+seq1->seq.l], "$", 1); //append $

		   if (useRC==1){
			  t=t+seq1->seq.l+1;

			  if (t>=Tsize){
				 while (t>Tsize){Tsize=Tsize+Tbufsize;}
				 self->T=realloc(self->T, Tsize*sizeof(char));
				 self->TIiv=realloc(self->TIiv, Tsize*sizeof(PyObject*));
			  }

			  rc(seq1->seq.s, seq1->seq.l);
			  memcpy(&self->T[f+l+1], seq1->seq.s, seq1->seq.l); //append the reverse complement
			  memcpy(&self->T[t-1], "$", 1); //append $
		   }

		   for (; f<t; f++) { //keep a Vertex pointer for every suffix
			  self->TIiv[f]=vid;
		   }
		} else {
			f=t;
		}
		
		if (vid==NULL) {
			PyErr_SetString(GSAError, "Failed to add vertex");
			return -1;
		}

	}
	kseq_destroy(seq1);
	gzclose(fp1);

	//Adding edges to graph
	if (isgraph1!=0){
		for (l=0; l<nrofedges; l++){ //TODO: change l here to unsigned
			PyObject *oid1 = PyInt_FromLong(edges_from[l]);
			Vertex *vFrom=(Vertex*) PyDict_GetItem(mapping1,oid1);
			if (vFrom==NULL){
				Py_DECREF(oid1);
				fprintf(stderr,"Failed to retrieve Vertex for id %u.\n",edges_from[l]);
				return -1;
			}
			Py_DECREF(oid1);

			PyObject *oid2 = PyInt_FromLong(edges_to[l]);
			Vertex *vTo=(Vertex*) PyDict_GetItem(mapping1,oid2);
			if (vTo==NULL){
				Py_DECREF(oid2);
				fprintf(stderr,"Failed to retrieve Vertex for id %u.\n",edges_to[l]);
				return -1;
			}
			Py_DECREF(oid2);

			uint8_t o=orientations[l];
			add_edge_c(self->graph, vFrom, vTo, o);
		}
		Py_DECREF(mapping1);
	}

	free(edges_from);
	free(edges_to);
	free(orientations);

	return f;
}







static int
buildGSA(GSA *self, char *input1, char* input2){

	self->graph = (Graph *)Graph_New(&GraphType, NULL, NULL);
	
	PyObject *nametuple;
	nametuple = Py_BuildValue("(s,s)",input1, input2);
	
	if (Graph_init(self->graph, nametuple, NULL)!=0){
		Py_DECREF(nametuple);
		return -1;
	}
	
	idx_t i=0,f=0,n=0;

	if (!(f=loadFile(self, input1, 0, f, 0))){
		PyErr_SetString(GSAError, "Failed to load first (g)fasta file.");
		return -1;
	}
	
	self->sep=f-1; //-1 for the splitting on the last sentinel of input1

	if (!(n=loadFile(self, input2, 1, f, self->rcindex==1 ? 1 : 0))){
		PyErr_SetString(GSAError, "Failed to load second (g)fasta file.");
		return -1;
	}
	

	self->n=n;
	self->orgn=self->n;
	self->T[self->n]='\0';
	
	//printf("Total number of bases in index: "PYIDXFORMAT"\n",(idx_t)self->n);

	//Calculate SA using SAIS
	self->SA=malloc(sizeof(idx_t)*self->n);
	if (self->SA==NULL){
		fprintf(stderr,"Failed to allocate enough memory for SA.\n");
		PyErr_SetString(GSAError, "Failed to allocate enough memory for SA.");
		return -1;
	}

	printf("building suffix array... ");
	fflush(stdout);

#if LARGEIDX == 1
	if (divsufsort64(self->T, self->SA, self->n)!=0){
#else
	if (divsufsort(self->T, self->SA, self->n)!=0){
#endif
		fprintf(stderr,"SAIS failed.\n");
		PyErr_SetString(GSAError, "SAIS-LCP failed");
		return -1;
	}
	printf("done.\n");

	printf("allocating memory for inverse of suffix array... ");
	fflush(stdout);
	self->SAi = malloc(sizeof(idx_t)*(self->n)); //inverse of SA
	if (self->SAi==NULL){
		fprintf(stderr,"Failed to allocate enough memory for SAi.\n");
		return -1;
	}
	printf("done.\n");

	printf("computing inverse of suffix array... ");
	fflush(stdout);
	//fill the inverse of SA
	for (i=0; i<self->n; i++) {
		self->SAi[self->SA[i]]=i;
	}
	printf("done.\n");

	printf("allocating memory for lcp array... ");
	fflush(stdout);
	self->LCP=malloc(sizeof(int)*self->n);
	if (self->LCP==NULL){
		fprintf(stderr,"Failed to allocate enough memory for LCP.\n");
		PyErr_SetString(GSAError, "Failed to allocate enough memory for LCP.");
		return -1;
	}
	printf("done.\n");

	printf("computing lcp array... ");
	fflush(stdout);
	compute_lcp(self->T, self->SA, self->SAi, self->LCP, self->n);
	//compute_lcp2(self->T, self->SA, self->LCP, self->n);

	self->LCP[0]=0;
	
	printf("done.\n");
	fflush(stdout);

	return 0;
}

static Vertex * add_vertex_c(Graph *self, int input_origin, int indexed, char * coord_origin, char * coord_contig, idx_t contig_start, idx_t contig_end, idx_t saoffset, idx_t rcsaoffset)
{
    Vertex *v = (Vertex *)Vertex_new(&VertexType, NULL, NULL);
    
    if (Vertex_init_c(v, self->vi, input_origin, indexed, coord_origin, coord_contig, contig_start, contig_end, saoffset, rcsaoffset)!=0){
        Py_DECREF(v);
        return NULL;
    }
    
    PyObject *k=Py_BuildValue("i", self->vi);
    
    if (k==NULL) {
        Py_DECREF(v);
        Py_DECREF(k);
        return NULL;
    }
    
    PyDict_SetItem(self->vertices, k, (PyObject *)v);
    
    Py_DECREF(k);
    
    self->vi++;
    
    return v;
}

static Edge * add_edge_c(Graph *self, Vertex* source, Vertex* target, uint8_t orientation)
{
    Edge *e;
    
    e = (Edge *)Edge_new(&EdgeType, NULL, NULL);
    
    if (Edge_init_c(e, source, target, orientation)!=0){
        return NULL;
    }
    
	if (orientation==0) { // --> -->
    	PySet_Add(source->edges_from,(PyObject*)e);
    	PySet_Add(target->edges_to,(PyObject*)e);
    	PySet_Add(self->edges,(PyObject*)e);
	}
	
	if (orientation==1) { // --> <--
    	PySet_Add(source->edges_from,(PyObject*)e);
    	PySet_Add(target->edges_from,(PyObject*)e);
    	PySet_Add(self->edges,(PyObject*)e);
	}
	
	if (orientation==2) { // <-- -->
    	PySet_Add(source->edges_to,(PyObject*)e);
    	PySet_Add(target->edges_to,(PyObject*)e);
    	PySet_Add(self->edges,(PyObject*)e);
	}

    self->ei++;
    
    return e;
}




static void updateInterval(GSA* self, idx_t iStart, idx_t length){
	
	idx_t i, l;
	
	//TODO: for now scan over all SA values, instead maybe lookup only values that are in specific affected interval
	for (i=0; i<self->n; i++) {
	
		if (((self->SA[i]>iStart) && (self->SA[i]<iStart+length)) ){
			continue;
		}
	
		if ((self->SA[i]<iStart) && (self->SA[i]+self->LCP[i]>iStart)){
			l=i;
			idx_t tmpSA=self->SA[i];
			idx_t tmpl=self->LCP[i];
						
			while ((self->LCP[l] >= iStart-tmpSA) && (l>0) ) {
				self->SAi[self->SA[l-1]]=l;
				self->SA[l]=self->SA[l-1];
				self->LCP[l]=self->LCP[l-1];
				l--;
			}
			self->SAi[tmpSA]=l;
			self->SA[l]=tmpSA;
			self->LCP[l+1]=iStart-tmpSA;
			
			if ((tmpl<self->LCP[i+1]) && (i<self->n-1)){
				self->LCP[i+1]=tmpl;
			}
			continue;
		}

		if (i<self->n-1){
			if ((self->SA[i]<iStart) && (self->SA[i]+self->LCP[i+1]>iStart)){
				self->LCP[i+1]=iStart-self->SA[i];
				continue;
			}
		}
	}
}


static PyObject* updateVertex(GSA* self, PyObject* args, PyObject* kwds)
{
	idx_t length;
	Vertex * v;

	if (!PyArg_ParseTuple(args, "O", &v)) {
		PyErr_SetString(GSAError, "Failed to parse input parameters, expected Vertex.");
		return NULL;
	}
	
	length=v->contig_end-v->contig_start;
	updateInterval(self, v->saoffset,length);
	
	if (v->rcsaoffset!=v->saoffset){
		updateInterval(self, v->rcsaoffset, length);
	}

	Py_INCREF(Py_None);
	return Py_None;
}


/*
		TODO: use SAi to find which values have to be resorted. Since we know that we don't have to 
		reorder all values of SA. Just the values preceding aStart and bStart with a max of the max LCP value!
*/

static PyObject* update(GSA* self, PyObject* args, PyObject* kwds)
{
	idx_t aStart, length, i;
	Vertex * v;
	if (!PyArg_ParseTuple(args, PYIDXTYPE PYIDXTYPE "O", &aStart, &length, &v)) {
		PyErr_SetString(GSAError, "Failed to parse input parameters.");
		return NULL;
	}

	for (i=0; i<length; i++){
		self->T[aStart+i]=tolower(self->T[aStart+i]);
		self->TIiv[aStart+i]=v; //make the interval point to this merged node
	}

	updateInterval(self, aStart, length);
	
	Py_INCREF(Py_None);
	return Py_None;
}

static PyObject* relabel(GSA *self, PyObject *args){
	Vertex* v;
	
	if (! PyArg_ParseTuple(args, "O", &v)) {
		PyErr_SetString(GSAError, "Failed to parse input parameters, expect a vertex.");
    	return NULL;
    }
	
	idx_t i;
	
	for (i=v->saoffset; i<v->saoffset+(v->contig_end-v->contig_start); i++){
		self->TIiv[i]=v;
	}
	
	if (self->rcindex==1 && v->saoffset!=v->rcsaoffset) {
		for (i=v->rcsaoffset; i<v->rcsaoffset+(v->contig_end-v->contig_start); i++){
			self->TIiv[i]=v;
		}
	}
	Py_INCREF(Py_None);
	return Py_BuildValue("O",Py_None);
}


/*
 
 Returns a new GSA instance based on the suffixes that correspond to the supplied set of Vertices
 and update the original index. 

 TODO: check for empty set!
 
 */
static GSA* extract(GSA *self, PyObject *args)
{
    PyObject* textids;
    
    if (! PyArg_ParseTuple(args, "O", &textids)) {
        PyErr_SetString(GSAError, "Failed to parse input parameters, expect set of text ids.");
        return NULL;
    }

    //create a new index
    GSA * nGSA=(GSA *)GSA_new(&GSAType, NULL, NULL);
 
    //iterate over all specified text ids
    PyObject *iter=PyObject_GetIter(textids);
    
    if (iter==NULL){
		printf("iter==NULL\n");
        PyErr_SetString(PyExc_TypeError, "Unable to iterate over set of text ids");
        Py_DECREF(nGSA);
        return NULL;
    }
    
    //PyObject *tup;
    
    //tmp array that indicates which values of SA have to be extracted
	//TODO: can use less space by reusing TIiv with additional datastructure here
    Vertex **tmp=(Vertex **) calloc(self->n,sizeof(Vertex *));
    
    idx_t n=0,i=0,x=0;

	//iterate over set of tuples
	//each tuple has a vertex and the relative orientation wrt the merged neighbor
	Vertex *v;
	while ( (v=(Vertex*)PyIter_Next(iter)) )
	{
		//int o;

		/*if (! PyArg_ParseTuple(tup, "O", &v)) {
			PyErr_SetString(PyExc_TypeError, "Unable to parse tuple of vertex and orientation");
			free(tmp);
			Py_DECREF(iter);
			Py_DECREF(nGSA);
			return NULL;
		}*/
		
		//orientation is the same, or we want both
		//if (o==1 || o==2) {
		if (v->indexed==0 || v->indexed==2) {

			//printf("%lu --> Extracting normal from: %d to %d length: %d\n",v->id,v->saoffset,v->saoffset+(v->contig_end-v->contig_start),v->contig_end-v->contig_start);

			assert(v->contig_start>=0);
			assert(v->contig_end>=0);
			assert(v->saoffset>=0);
			assert(v->contig_end-v->contig_start>=0);

			for (i=v->saoffset; i<v->saoffset+(v->contig_end-v->contig_start); i++){
				assert(i<self->orgn);
				assert(i>=0);
				if (self->SAi[i]>=self->n){
					//printf("n=%lu sai=%d i=%d\n",self->n,self->SAi[i], i);
					printf("Vertex %u could not be extracted from index.\n",v->id);
					PyErr_SetString(PyExc_TypeError, "Vertex not indexed.");
					return NULL;
				}
				assert(self->SAi[i]<self->n);
				assert(self->SAi[i]>=0);

				if (tmp[self->SAi[i]]!=0){
					//fprintf(stderr, "ERROR tmp[self->SAi[i]]!=0 --> %u i="PYIDXFORMAT"...\n", tmp[self->SAi[i]]->id, i);
					PyErr_SetString(PyExc_TypeError, "Index corrupt. SAi for SA not unique.");
					return NULL;
				}

				tmp[self->SAi[i]]=v;

				if (v==NULL){
					PyErr_SetString(PyExc_TypeError, "Vertex reference lost!.");
					return NULL;
				}

				x++;
			}

			n=n+(v->contig_end - v->contig_start);
		}

		//opposite orientation, or we want to have both
		if (v->indexed==1 || v->indexed==2) {

			assert(self->rcindex==1);
			assert(v->rcsaoffset>0);

			for (i=v->rcsaoffset; i<v->rcsaoffset+(v->contig_end-v->contig_start); i++) {
				assert(i<self->orgn);
				assert(i>=0);
				if (self->SAi[i]>self->n){
					printf("Vertex %u (rc) could not be extracted from index.",v->id);									
					PyErr_SetString(PyExc_TypeError, "Index corrupt. self->SAi[i]>=self->n");
					return NULL;
				}
				assert(self->SAi[i]<self->n);
				assert(self->SAi[i]>=0);
				if (tmp[self->SAi[i]]!=0){
					//fprintf(stderr, "ERROR tmp[self->SAi[i]]!=0 --> %u i="PYIDXFORMAT"...\n", tmp[self->SAi[i]]->id, i);
					PyErr_SetString(PyExc_TypeError, "Index corrupt. SAi for SA not unique.");
					return NULL;
				}

				tmp[self->SAi[i]]=v;

				if (v==NULL){
					PyErr_SetString(PyExc_TypeError, "Vertex reference lost!.");
					return NULL;
				}

				x++;
			}
			assert(v->contig_end>v->contig_start);
			n=n+(v->contig_end-v->contig_start);
		}

		Py_DECREF(v);
	}

	Py_DECREF(iter);

	assert(n==x);

	nGSA->n=n;
    
    nGSA->SA=malloc(n * sizeof(idx_t));
    if (nGSA->SA==NULL) {
        free(tmp);
		//printf("Failed trying to allocate "PYIDXFORMAT" integers of memory.\n",n);
        PyErr_SetString(PyExc_TypeError, "Could not allocate enough memory for new index.");
        return NULL;
    }
   
    nGSA->LCP=malloc(n * sizeof(int));
    if (nGSA->LCP==NULL) {
        free(tmp);
		//printf("Failed trying to allocate "PYIDXFORMAT" integers of memory.\n",n);
        PyErr_SetString(PyExc_TypeError, "Could not allocate enough memory for new index.");
        return NULL;
    }
    
    idx_t ni=(self->n)-n;
    
    idx_t *_SA=malloc(ni * sizeof(idx_t));
    if (_SA==NULL) {
        free(tmp);
        PyErr_SetString(PyExc_TypeError, "Could not allocate enough memory for new index _SA");
        return NULL;
    }
   
    unsigned int *_LCP=malloc(ni * sizeof(int));
    if (_LCP==NULL) {
        free(tmp);
        PyErr_SetString(PyExc_TypeError, "Could not allocate enough memory for new index LCPi.");
        return NULL;
    }

    idx_t c=0,last=0,trace=0;
    idx_t ci=0,lasti=0,tracei=0;
	
	unsigned int minlcp=0, minlcpi=0;
	
    for (i=0; i<self->n; i++) {
		if (tmp[i]!=0) {
			if ((self->LCP[i]<minlcpi) || tracei==1) {
				minlcpi=self->LCP[i];
				tracei=0;
			}

			//copy value to new SA
			nGSA->SA[c]=self->SA[i];
			//determine new value for LCP
			if (c==0) { //if it's the first record in the new LCP array start with -1
				assert(minlcp==0);
				nGSA->LCP[c]=minlcp;
				trace=1;
			} else {
				if (i-last>1) { //if values were not consecutive
					if (self->LCP[i]<minlcp) {
						minlcp=self->LCP[i];
					}
					nGSA->LCP[c]=minlcp;
					trace=1;
				} else { //if values were consecutive
					nGSA->LCP[c]=self->LCP[i];
				}
			}
			last=i;
			//update SAi with c
			self->SAi[self->SA[i]]=c;

			//update TIiv with reference to new vertices
			self->TIiv[self->SA[i]]=tmp[i];
			c++;
		} else {
			if ((self->LCP[i]<minlcp) || trace==1) {
				minlcp=self->LCP[i];
				trace=0;
			}
				
			//keep track of the other values as well
			_SA[ci]=self->SA[i];
			//determine new value for LCP
			if (ci==0) { //if it's the first record in the new LCP array start with -1
				assert(minlcpi==0);
				_LCP[ci]=minlcp;
				tracei=1;
			} else {
				if (i-lasti>1) { //if values were not consecutive
					if (self->LCP[i]<minlcpi) {
						minlcpi=self->LCP[i];
					}
					_LCP[ci]=minlcpi;
					tracei=1;
				} else { //if values were consecutive
					_LCP[ci]=self->LCP[i];
				}
			}
			lasti=i;
			//update SAi with ci
			self->SAi[self->SA[i]]=ci;
			
			ci++;
		}
    }

	//Keep reference to the SAi and TIiv for subindices
	nGSA->orgn=self->orgn;
	nGSA->SAi=self->SAi;
	nGSA->TIiv=self->TIiv;
	nGSA->T=self->T;
	nGSA->rcindex=self->rcindex;
	nGSA->sep=self->sep;
	nGSA->root=self->root; //keep reference to the root index, so we can decref at dealloc
	Py_INCREF(self->root); //make sure SAi and TIiv are not deallocated
	
	free(self->SA);
	free(self->LCP);
	self->SA=_SA;
	self->LCP=_LCP;
	self->n=ni;

	free(tmp);
	
	return nGSA;
}


//TODO: remove dependency on TI for this function --> only works for aligning a set of max two sequences, so use self->sep to distinguish between them!
//cannot be used to segment sets of sequences!!!
static PyObject * split(GSA *self, PyObject *args){
    GSA * gsa_a;
    GSA * gsa_b;
    
    //create two new index objects
    gsa_a=(GSA *)GSA_new(&GSAType, NULL, NULL);
    gsa_b=(GSA *)GSA_new(&GSAType, NULL, NULL);
    
    idx_t idx,ml;
    
    if (! PyArg_ParseTuple(args, "ii", &idx, &ml)) {
        PyErr_SetString(GSAError, "Failed to parse input parameters");
        return NULL;
    }
    
    idx_t startInT1,startInT2,endInT1,endInT2,ss,lcp,tid,tid1,tid2,i,bi=0,ai=0,minlcp1=-1,minlcp0=-1,last_set=-1;
    startInT1 = self->SA[idx];
    startInT2 = self->SA[idx-1];
    endInT1 = self->SA[idx]+ml;
    endInT2 = self->SA[idx-1]+ml;
    
    //allocate memory for new index datastructures
    if (startInT2>startInT1) {
        gsa_a->n=(startInT1+(startInT2 - self->sep));
        gsa_a->SA=malloc(gsa_a->n * sizeof(int));
        gsa_a->LCP=malloc(gsa_a->n * sizeof(int));

        gsa_b->n=(self->sep-endInT1) + (self->n-endInT2);
        gsa_b->SA=malloc(gsa_b->n * sizeof(int));
        gsa_b->LCP=malloc(gsa_b->n * sizeof(int));
    } else {
        gsa_a->n=(startInT2 + (startInT1 - self->sep));
        gsa_a->SA=malloc(gsa_a->n * sizeof(int));
        gsa_a->LCP=malloc(gsa_a->n * sizeof(int));
        
        gsa_b->n=(self->sep-endInT2) + (self->n-endInT1);
        gsa_b->SA=malloc(gsa_b->n * sizeof(int));
        gsa_b->LCP=malloc(gsa_b->n * sizeof(int));
    }
    
    tid1=self->TIiv[self->SA[idx]]->id;
    tid2=self->TIiv[self->SA[idx-1]]->id;
    
    if (tid1==tid2) {
        PyErr_SetString(GSAError, "Invalid match. Match within instead of between sequences.");
        return NULL;
    }

    //fprintf(stderr,"startInT1=%d endInT1=%d startInT2=%d endInT2=%d\n",startInT1,endInT1,startInT2,endInT2);
    
    for (i=0; i<self->n; i++) {
        ss=self->SA[i];
        lcp=self->LCP[i];
        tid=self->TIiv[self->SA[i]]->id;
        
        //fprintf(stderr,"%d %d %d %d %d %d %d %d %d %d\n",i,ss,lcp,tid,tid1,tid2,minlcp0,minlcp1,ai,bi);
        
        if ( ((tid==tid1) && (ss<startInT1)) || ((tid==tid2) && (ss<startInT2)) ) {
            //printf("--> a\n");
            
            if (ai==0){ //first record
                gsa_a->LCP[ai]=-1;
                minlcp0=-2;
            } else {
                if (last_set==0) { //consecutive values, no need to recalc LCP
                    gsa_a->LCP[ai]=lcp;
                    minlcp0=-2;
                } else { //not consecutive, use minlcp
                    if ((lcp<minlcp0) || (minlcp0<0)) {
                        minlcp0=lcp;
                    }
                    gsa_a->LCP[ai]=minlcp0;
                    minlcp0=-2;
                }
            }
            gsa_a->SA[ai]=ss;
            //gsa_a->TI[ai]=tid;
            last_set=0;
            ai++;
        }
        
        if (((tid==tid1) && (ss>=endInT1)) || ((tid==tid2) && (ss>=endInT2))) {
            //printf("--> b\n");
            
            if (bi==0){ //first record
                gsa_b->LCP[bi]=-1;
                minlcp1=-2;
            } else {
                if (last_set==1){ //consecutive values, no need to recalc LCP
                    gsa_b->LCP[bi]=lcp;
                    minlcp1=-2;
                } else { //not consecutive, use minlcp
                    if ((lcp<minlcp1) || (minlcp1<0)) {
                        minlcp1=lcp;
                    }
                    gsa_b->LCP[bi]=minlcp1;
                    minlcp1=-2;
                }
            }
            gsa_b->SA[bi]=ss;
            //gsa_b->TI[bi]=tid;
            last_set=1;
            bi++;
        }
        
        if (((tid==tid1) && ((ss>=startInT1) && (ss<endInT1))) || ((tid==tid2) && ((ss>=startInT2) && (ss<endInT2))) ) {
            //printf("--> NONE\n");
            //check min for both minlcp0 and minlcp1
            if ((lcp<minlcp0) || ((minlcp0==-2) || (minlcp0==-1))) { minlcp0=lcp;}
            if ((lcp<minlcp1) || ((minlcp1==-2) || (minlcp1==-1))) { minlcp1=lcp;}
            last_set=2;
            continue;
        }
    
    
        if ((last_set==1) && (minlcp0!=-1)) {
            //check min for minlcp0
            if ((lcp<minlcp0) || (minlcp0==-2)) {minlcp0=lcp;}
        }
    
        if ((last_set==0) && (minlcp1!=-1)) {
            //check min for minlcp1
            if ((lcp<minlcp1) || (minlcp1==-2)) {minlcp1=lcp;}
        }
    }
    
    return Py_BuildValue("OO", gsa_a, gsa_b);
}



#ifndef PyMODINIT_FUNC	/* declarations for DLL import/export */
#define PyMODINIT_FUNC void
#endif
PyMODINIT_FUNC
initGSA(void)
{
    PyObject* m;
    
    if (PyType_Ready(&GSAType) < 0)
        return;
    
    if (PyType_Ready(&GraphType) < 0)
        return;
    
    if (PyType_Ready(&VertexType) < 0)
        return;
    
    if (PyType_Ready(&EdgeType) < 0)
        return;
    
    m = Py_InitModule3(PYSONAME, gsa_methods,
                       "Module that ...");
    
    GSAError = PyErr_NewException("GSA.error", NULL, NULL);
    Py_INCREF(GSAError);
    PyModule_AddObject(m, "error", GSAError);

    Py_INCREF(&GraphType);
    PyModule_AddObject(m, "Graph", (PyObject *)&GraphType);
    
    Py_INCREF(&VertexType);
    PyModule_AddObject(m, "Vertex", (PyObject *)&VertexType);
    
    Py_INCREF(&EdgeType);
    PyModule_AddObject(m, "Edge", (PyObject *)&EdgeType);
    
    Py_INCREF(&GSAType);
    PyModule_AddObject(m, "index", (PyObject *)&GSAType);
}

PyMODINIT_FUNC
initGSA_64(void)
{
	initGSA();
}
