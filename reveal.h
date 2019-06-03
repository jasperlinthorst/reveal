#ifndef REVEAL
#define REVEAL
#define QUEUE_BUF 10000

void *aligner(void *arg);

#ifdef SA64
#define saidx_t int64_t
#define lcp_t uint32_t
#else
#define saidx_t int32_t
#define lcp_t int32_t
#endif

#define diff(a,b) (a > b) ? a - b : b - a;

typedef struct
{
    PyObject_HEAD
    char             * T;   //initial input Text
    saidx_t          * SA;  //Suffix array
    saidx_t          * SAi; //Inverse of suffix array
    lcp_t            * LCP; //LCP array, lcp[0]==-1
    uint16_t         * SO;  //Array indicating for each suffix which sample it originated from (max 2**16 samples!)
    saidx_t            n;   //length of SA and LCP after extract
    saidx_t            nT;   //length of the original input text T
    saidx_t          * nsep;    //array of integers pointing to sentinels that seperate samples within the input T
    int depth; //depth within the hierarchical alignment tree
    int nsamples; //number of samples in T
    char* safile;
    char* lcpfile;
    int rc; //whether the index use the reverse complement of the qry part
    int cache; //whether to store the suffix array and lcp array after construction, so they can be reused
    PyObject * main; //main index
    PyObject * samples; //list of sample names that are contained in the index
    PyObject * nodes; //list of intervals in T that are associated with this index
    PyObject * left_node;
    PyObject * right_node;
    PyObject * skipmums;
} RevealIndex;

RevealIndex* newIndex(void);

typedef struct
{
    int threadid;
    PyObject * mumpicker; //callback function that return the best exact match from a list of exact matches
    PyObject * graphalign; //callback that updates the interval tree and graph for the alignment
    int wpen;
    int wscore;
    int minl;
    int minn;
} RevealWorker;

typedef struct
{
    lcp_t l; //length of the exact match
    saidx_t *sp; //array of starting positions
    //int *ni; //array of sample index to which sp corresponds
    int n;   //number of samples in which the exact match occurs
    long long score; //can be negative
    unsigned long long penalty;
} RevealMultiMUM;

int getbestmum(RevealIndex *index, RevealMultiMUM *mum, int w_penalty, int w_score);
int getbestmultimum(RevealIndex *index, RevealMultiMUM *mmum, int min_n);

PyObject * getmultimums(RevealIndex *index, PyObject *args, PyObject *keywds);
PyObject * getmultimems(RevealIndex *index, PyObject *args, PyObject *keywds);
PyObject * getmums(RevealIndex *index, PyObject *args, PyObject *keywds);
PyObject * splitindex(RevealIndex *index, PyObject *args, PyObject *keywds);
PyObject * extract(RevealIndex *index, PyObject *args, PyObject *keywds);
PyObject * getscoredmums(RevealIndex *index, PyObject *args, PyObject *keywds);

#endif