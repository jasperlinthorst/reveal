#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <pthread.h>
#include <semaphore.h>
#include <unistd.h>
#include <time.h>
#include "Python.h"
#include "reveal.h"

/* The mutex lock */
extern pthread_mutex_t mutex, python;

extern RevealIndex **index_queue;
extern int maxqsize,qsize,qstart,aw,nmums,err_flag,die;

/* pops an index of the queue */
RevealIndex* pop_index() {
    //fprintf(stderr,"qend=%d qstart=%d qsize=%d\n",qsize,qstart,qsize-qstart);
    //if(qsize > qstart) { FIFO
    if(qsize > 0) { //LIFO
        return index_queue[--qsize]; //LIFO
        //return index_queue[qstart++]; //FIFO
    }
    else { /* Error buffer empty */
        return NULL;
    }
}

/* puts an index on the queue */
int push_index(RevealIndex *idx) {
    if(qsize == maxqsize) {
        RevealIndex **newq;//realloc the index_queue to be able to hold more indices
        newq=realloc(index_queue,(maxqsize+QUEUE_BUF)*sizeof(RevealIndex *));
        if (newq==NULL) {
            fprintf(stderr,"Failed to allocate memory for index queue.\n");
            return -1; //out of memory?
        } else {
            index_queue=newq;
        }
        maxqsize=maxqsize+QUEUE_BUF;
    }
    index_queue[qsize] = idx;
    qsize++;
    return 0;
}

PyObject * getmums(RevealIndex *index, PyObject *args, PyObject *keywds){
    uint32_t minl=0;
    
    if (args!=NULL) {
        if (!PyArg_ParseTuple(args, "i", &minl))        
            return NULL;
    }
    
    saidx_t i=0,aStart,bStart;
    uint32_t lb,la;
    PyObject *mums=PyList_New(0);

    for (i=1;i<index->n;i++){
        if (index->LCP[i]<minl){
            continue;
        }
	if (index->SA[i]>index->nsep[0] == index->SA[i-1]>index->nsep[0]){ //repeat
	    continue;
	}
	if (index->SA[i]<index->SA[i-1]) {
	    aStart=index->SA[i];
	    bStart=index->SA[i-1];
	} else {
	    aStart=index->SA[i-1];
	    bStart=index->SA[i];
	}
	if (aStart>0 && bStart>0){ //if not it has to be maximal!
	    if (!((index->T[aStart-1]!=index->T[bStart-1]) || (index->T[aStart-1]=='N') || (index->T[aStart-1]=='$') || (islower(index->T[aStart-1])) )) {
		continue; //not maximal
	    }
	}
	if (i==index->n-1) { //is it the last value in the array, then only check predecessor
	    lb=index->LCP[i-1];
	    la=0;
	} else {
	    lb=index->LCP[i-1];
	    la=index->LCP[i+1];
	}
	if (lb>=index->LCP[i] || la>=index->LCP[i]){
	    continue;//not unique
	}

	//match is not a repeat and is maximally unique
        
#ifdef SA64
	PyObject *mum=Py_BuildValue("I,i,(L,L)",index->LCP[i],2,aStart,bStart);
#else
	PyObject *mum=Py_BuildValue("I,i,(i,i)",index->LCP[i],2,aStart,bStart);
#endif
	
        if (PyList_Append(mums,mum)==0){
	    Py_DECREF(mum);
	} else {
            Py_DECREF(mum); //append increments reference count!
	    return NULL;
	}
    }
    return mums;
}

PyObject * getmems(RevealIndex *index, PyObject *args, PyObject *keywds){
    
    uint32_t minl=0;
    
    if (args!=NULL) {
        if (!PyArg_ParseTuple(args, "i", &minl))
            return NULL;
    }
    
    saidx_t i=0,aStart,bStart;
    uint32_t lb,la,uniq;
    PyObject *mems=PyList_New(0);
    for (i=1;i<index->n;i++){
        if (index->LCP[i]<minl){
            continue;
        }
	if (index->SA[i]>index->nsep[0] == index->SA[i-1]>index->nsep[0]){ //repeat within same genome
	    continue;
	}
	if (index->SA[i]<index->SA[i-1]) {
	    aStart=index->SA[i];
	    bStart=index->SA[i-1];
	} else {
	    aStart=index->SA[i-1];
	    bStart=index->SA[i];
	}
	if (aStart>0 && bStart>0){ //if not it has to be maximal!
	    if (!((index->T[aStart-1]!=index->T[bStart-1]) || (index->T[aStart-1]=='N') || (index->T[aStart-1]=='$') || (islower(index->T[aStart-1])) )) {
		continue; //not maximal
	    }
	}
	if (i==index->n-1) { //is it the last value in the array, then only check predecessor
	    lb=index->LCP[i-1];
	    la=0;
	} else {
	    lb=index->LCP[i-1];
	    la=index->LCP[i+1];
	}
	if (lb>=index->LCP[i] || la>=index->LCP[i]){
            uniq=0;
	} else {
            uniq=1;
        }
        
	//match is not a repeat and is maximal

#ifdef SA64
	PyObject *mem=Py_BuildValue("I,i,(L,L),i",index->LCP[i],2,aStart,bStart,uniq);
#else
	PyObject *mem=Py_BuildValue("I,i,(i,i),i",index->LCP[i],2,aStart,bStart,uniq);
#endif

        if (PyList_Append(mems,mem)==0){
	    Py_DECREF(mem);
	} else {
            Py_DECREF(mem); //append increments reference count!
	    return NULL;
	}
    }
    return mems;
}

PyObject * getscoredmums(RevealIndex *index, PyObject *args, PyObject *keywds){
    unsigned long long i=0,aStart,bStart;
    unsigned long long start1=-1,end1=-1,start2=-1,end2=-1; //TODO: add to keywds
    long long score;
    
    uint32_t lb,la;
    int w_score=2;
    int w_penalty=1;
    int penalize;
    
    if (PyList_Size(index->nodes)==2){
        penalize=1;

        PyObject *node1=PyList_GetItem(index->nodes,0);
        start1=PyInt_AS_LONG(PyTuple_GetItem(node1,0));
        end1=PyInt_AS_LONG(PyTuple_GetItem(node1,1));

        PyObject *node2=PyList_GetItem(index->nodes,1);
        start2=PyInt_AS_LONG(PyTuple_GetItem(node2,0));
        end2=PyInt_AS_LONG(PyTuple_GetItem(node2,1));

    } else {
        penalize=0;
    }
    
    PyObject *mums=PyList_New(0);
    for (i=1;i<index->n;i++){
	if (index->SA[i]>index->nsep[0] == index->SA[i-1]>index->nsep[0]){ //repeat
	    continue;
	}
	if (index->SA[i]<index->SA[i-1]) {
	    aStart=index->SA[i];
	    bStart=index->SA[i-1];
	} else {
	    aStart=index->SA[i-1];
	    bStart=index->SA[i];
	}
	if (aStart>0 && bStart>0){ //if not it has to be maximal!
	    if (!((index->T[aStart-1]!=index->T[bStart-1]) || (index->T[aStart-1]=='N') || (index->T[aStart-1]=='$') || (islower(index->T[aStart-1])) )) {
		continue; //not maximal
	    }
	}
	if (i==index->n-1) { //is it the last value in the array, then only check predecessor
	    lb=index->LCP[i-1];
	    la=0;
	} else {
	    lb=index->LCP[i-1];
	    la=index->LCP[i+1];
	}
	if (lb>=index->LCP[i] || la>=index->LCP[i]){
	    continue;//not unique
	}

        long long penalty=0;
        unsigned long long lpenalty=0;
        unsigned long long tpenalty=0;
        
        if (index->depth>0 && penalize==1){
            
            if (start1<=aStart && end1>=aStart){
                
                lpenalty=diff(aStart,start1);
                //lpenalty= abs( (aStart-start1)-(bStart-start2) ); //leading penalty
                tpenalty=diff((end1-(aStart+index->LCP[i])),(end2-(bStart+index->LCP[i])));
                //tpenalty=abs((end1-(aStart+index->LCP[i]))-(end2-(bStart+index->LCP[i]))); //trailing penalty
                
                assert(lpenalty>=0);
                assert(tpenalty>=0);
                if (lpenalty>tpenalty){
                    penalty=tpenalty;
                } else {
                    penalty=lpenalty;
                }
                //penalty=tpenalty+lpenalty;
                assert(penalty>=0);
            } else {
                
                lpenalty=diff((bStart-start1),(aStart-start2));
                //lpenalty=abs((bStart-start1)-(aStart-start2)); //leading penalty
                tpenalty=diff((end1-(bStart+index->LCP[i])), (end2-(aStart+index->LCP[i])));
                //tpenalty=abs((end1-(bStart+index->LCP[i]))-(end2-(aStart+index->LCP[i]))); //trailing penalty
                
                assert(lpenalty>=0);
                assert(tpenalty>=0);
                if (lpenalty>tpenalty){
                    penalty=tpenalty;
                } else {
                    penalty=lpenalty;
                }
                //penalty=tpenalty+lpenalty;
                assert(penalty>=0);
            }
        }
        
        score=(w_score*index->LCP[i])-(w_penalty*penalty);
        
#ifdef SA64
	PyObject *mum=Py_BuildValue("I,i,(L,L),L,K,K",index->LCP[i],2,aStart,bStart,score,lpenalty,tpenalty);
#else
	PyObject *mum=Py_BuildValue("I,i,(i,i),L,K,K",index->LCP[i],2,aStart,bStart,score,lpenalty,tpenalty);
#endif
        if (PyList_Append(mums,mum)==0){
	    Py_DECREF(mum);
	} else {
            Py_DECREF(mum); //append increments reference count!
	    return NULL;
	}
    }
    return mums;
}

int getlongestmum(RevealIndex *index, RevealMultiMUM *mum){
    saidx_t i=0,aStart,bStart;
    uint32_t lb,la;
    mum->u=1;
    mum->l=0;
    mum->score=0;
    mum->penalty=0;
    mum->n=2;
    for (i=1;i<index->n;i++){
        if (index->LCP[i]>mum->l){
            if (index->SA[i]>index->nsep[0] == index->SA[i-1]>index->nsep[0]){ //repeat
               continue;
            }
            if (index->SA[i]<index->SA[i-1]) {
                aStart=index->SA[i];
                bStart=index->SA[i-1];
            } else {
                aStart=index->SA[i-1];
                bStart=index->SA[i];
            }
            if (aStart>0 && bStart>0){ //if not it has to be maximal!
                if (!((index->T[aStart-1]!=index->T[bStart-1]) || (index->T[aStart-1]=='N') || (index->T[aStart-1]=='$') || (islower(index->T[aStart-1])) )) {
                    continue; //not maximal
                }
            }
            if (i==index->n-1) { //is it the last value in the array, then only check predecessor
                lb=index->LCP[i-1];
                la=0;
            } else {
                lb=index->LCP[i-1];
                la=index->LCP[i+1];
            }
            if (lb>=index->LCP[i] || la>=index->LCP[i]){
                continue;//not unique
            }
            //match is not a repeat and is maximally unique
            mum->l=index->LCP[i];
            mum->score=(mum->l*mum->n)-mum->penalty;
            mum->sp[0]=aStart;
            mum->sp[1]=bStart;
        }
    }
    return 0;
}

int getbestmum(RevealIndex *index, RevealMultiMUM *mum, int w_penalty, int w_score){
    saidx_t i=0;
    saidx_t aStart,bStart;
    saidx_t start1=-1,end1=-1,start2=-1,end2=-1;
    long long score;
    
    uint32_t lb,la,penalize;
    
    if (PyList_Size(index->nodes)==2){
        penalize=1;
        PyObject *node1=PyList_GetItem(index->nodes,0);

#ifdef SA64
        start1=PyLong_AsLong(PyTuple_GetItem(node1,0));
        end1=PyLong_AsLong(PyTuple_GetItem(node1,1));
        //fprintf(stderr,"start1=%lld end1=%lld\n",start1,end1);
#else 
        start1=PyInt_AS_LONG(PyTuple_GetItem(node1,0));
        end1=PyInt_AS_LONG(PyTuple_GetItem(node1,1));
        //fprintf(stderr,"start1=%d end1=%d\n",start1,end1);
#endif
        
        assert(start1>=0);
        assert(end1>=0);
        
        PyObject *node2=PyList_GetItem(index->nodes,1);

#ifdef SA64
        start2=PyLong_AsLong(PyTuple_GetItem(node2,0));
        end2=PyLong_AsLong(PyTuple_GetItem(node2,1));
        //fprintf(stderr,"start2=%lld end2=%lld\n",start2,end2);
#else 
        start2=PyInt_AS_LONG(PyTuple_GetItem(node2,0));
        end2=PyInt_AS_LONG(PyTuple_GetItem(node2,1));
        //fprintf(stderr,"start2=%d end2=%d\n",start2,end2);
#endif
        
        assert(start2>=0);
        assert(end2>=0);

    } else {
        penalize=0;
    }
    
    mum->u=1;
    mum->l=0;
    mum->score=-2147483648;
    mum->penalty=0;
    mum->n=2;

    for (i=1;i<index->n;i++){
        
        if (index->LCP[i]*w_score > mum->score){
            if (index->SA[i]>index->nsep[0] == index->SA[i-1]>index->nsep[0]){ //repeat
               continue;
            }
            if (index->SA[i]<index->SA[i-1]) {
                aStart=index->SA[i];
                bStart=index->SA[i-1];
            } else {
                aStart=index->SA[i-1];
                bStart=index->SA[i];
            }
            if (aStart>0 && bStart>0){ //if not it has to be maximal!
                if (!((index->T[aStart-1]!=index->T[bStart-1]) || (index->T[aStart-1]=='N') || (index->T[aStart-1]=='$') || (islower(index->T[aStart-1])) )) {
                    continue; //not maximal
                }
            }
            if (i==index->n-1) { //is it the last value in the array, then only check predecessor
                lb=index->LCP[i-1];
                la=0;
            } else {
                lb=index->LCP[i-1];
                la=index->LCP[i+1];
            }
            if (lb>=index->LCP[i] || la>=index->LCP[i]){
                continue;//not unique
            }
            //match is not a repeat and is maximally unique
            
            //penalize for the indel it creates in leading and trailing part
            long long penalty=0;
            long long lpenalty=0;
            long long tpenalty=0;

            if (index->depth>0 && penalize==1){
                
                if (start1<=aStart && end1>=aStart){

                    lpenalty=diff((aStart-start1),(bStart-start2));
                    tpenalty=diff((end1-(aStart+index->LCP[i])),(end2-(bStart+index->LCP[i])));
                    
                    //lpenalty=abs((aStart-start1)-(bStart-start2)); //leading penalty
                    //tpenalty=abs( (end1-(aStart+index->LCP[i]))-(end2-(bStart+index->LCP[i]))); //trailing penalty
                    
                    assert(lpenalty>=0);
                    assert(tpenalty>=0);

                    if (lpenalty>tpenalty){
                        penalty=tpenalty;
                    } else {
                        penalty=lpenalty;
                    }
                    //penalty=tpenalty+lpenalty;
                    assert(penalty>=0);
                } else {
                    //fprintf(stderr,"lcp=%d bStart=%d, start1=%lld, aStart=%d start2=%lld\n",index->LCP[i],bStart,start1,aStart,start2);
                    
                    lpenalty=diff((bStart-start1),(aStart-start2));

                    //fprintf(stderr,"a=%lld b=%lld\n",(bStart-start1),(aStart-start2));
                    //fprintf(stderr,"lpen=%lld\n",lpenalty);

                    tpenalty=diff((end1-(bStart+index->LCP[i])), (end2-(aStart+index->LCP[i])));
                    
                    //lpenalty=abs((bStart-start1)-(aStart-start2)); //leading penalty
                    //tpenalty=abs((end1-(bStart+index->LCP[i]))-(end2-(aStart+index->LCP[i]))); //trailing penalty
                    
                    assert(lpenalty>=0);
                    assert(tpenalty>=0);
                    if (lpenalty>tpenalty){
                        penalty=tpenalty;
                    } else {
                        penalty=lpenalty;
                    }
                    //penalty=tpenalty+lpenalty;
                    assert(penalty>=0);
                }
            }
            
            score=(w_score*index->LCP[i])-sqrt(w_penalty*penalty);
             
            if (score > mum->score){
                mum->score=score;
                mum->penalty=w_penalty*penalty;
                mum->l=index->LCP[i];
                mum->sp[0]=aStart;
                mum->sp[1]=bStart;
            }
        }
    }
    return 0;
}

int getbestmultimum(RevealIndex *index, RevealMultiMUM *mum, int min_n){
    mum->u=1;
    mum->l=0;
    mum->score=0;
    mum->penalty=0;
    mum->n=0;
    unsigned long long i=0;
    uint32_t j=0,la,lb;
    int *flag_so=calloc(index->nsamples, sizeof(int));
    int flag_maximal=0;
    int flag_unique=0;
    int flag_lcp;
    int x;
    
    for (i=index->nsamples-1;i<index->n;i++){
        if (index->LCP[i]>mum->l && (i==index->n-1 || index->LCP[i+1]<index->LCP[i])){
            int so;
            flag_maximal=0;
            flag_unique=0;
            flag_lcp=index->LCP[i];
            int n=0;
            for (j=0;j<(index->nsamples-1);j++){ //scan back n places in SA
                if (j==0){
                    assert(index->SO[index->SA[i]]>=0 && index->SO[index->SA[i]]<index->nsamples);
                    flag_so[index->SO[index->SA[i]]]=1;
                    n=1;
                }
                assert(i-j-1>=0);
                so=index->SO[index->SA[i-j-1]];
                assert(so>=0 && so<=index->nsamples);
                if (flag_so[so]==0){
                    flag_so[so]=1;
                } else {
                    j--;
                    break;
                }
                if (!flag_maximal){
                    if (index->SA[i-j]>0 && index->SA[i-j-1]>0){
                        assert((index->SA[i-j]-1)>=0);
                        assert((index->SA[i-j-1]-1)>=0);
                        if (index->T[index->SA[i-j]-1] != index->T[index->SA[i-j-1]-1] || islower(index->T[index->SA[i-j]-1]) || index->T[index->SA[i-j]-1]=='$' ){
                            flag_maximal=1;
                        }
                    } else {
                        flag_maximal=1;
                    }
                }
                if (index->LCP[i-j]<flag_lcp){
                    assert((i-j)>=0 && (i-j)<index->n);
                    flag_lcp=index->LCP[i-j];
                }
                n++;
            }
            
            if (n>=min_n && flag_maximal && flag_lcp>mum->l){
                //is it unique within the entire index?
                if (i==index->n-1) { //is it the last domain in the array, then only check predecessor
                    assert(i-(n-1)>=0);
                    lb=index->LCP[i-(n-1)];
                    la=0;
                } else if ((i-(n-1))==0) { //is it the first domain in the array, then only check successor
                    lb=0;
                    assert(i+1<index->n);
                    la=index->LCP[i+1];
                } else {
                    assert(i-(n-1)>=0);
                    lb=index->LCP[i-(n-1)];
                    assert(i+1<index->n);
                    la=index->LCP[i+1];
                }
                if (flag_lcp>lb && flag_lcp>la){
                    flag_unique=1;
                    mum->l=flag_lcp;
                    mum->n=n;
                    for (x=0;x<n;x++){
                        assert(i-x>=0);
                        mum->sp[x]=index->SA[i-x];
                    }
                }
            }
            for (j=0;j<index->nsamples;j++){ //memset?
                flag_so[j]=0; //reset so flags
            }
            flag_maximal=0; //reset maximal flag
        }
    }
    free(flag_so);
    return 0;
}

int ismultimum(RevealIndex * idx, uint32_t l, int lb, int ub, int * flag_so) {
    if (l>0){
        int j;
        memset(flag_so,0,((RevealIndex *) idx->main)->nsamples * sizeof(int));
        
        if (((RevealIndex *) idx->main)->nsamples==2){ //dont need SO in case of only two samples
            if ( (idx->SA[ub]>idx->nsep[0]) == (idx->SA[lb]>idx->nsep[0]) ){
                return 0;
            }
        } else {
            for (j=lb; j<ub+1; j++) { //has to occur in all samples once
                if (flag_so[idx->SO[idx->SA[j]]]==0){
                    flag_so[idx->SO[idx->SA[j]]]=1;
                } else {
                    return 0;
                }
            }
        }
        
        for (j=lb; j<ub; j++){ //check maximal
            if (idx->SA[j]==0){

                return 1; //success
            }
            if (idx->SA[j+1]==0){

                return 1; //success
            }
            if (idx->T[idx->SA[j]-1]!=idx->T[idx->SA[j+1]-1] || idx->T[idx->SA[j]-1]=='N' || idx->T[idx->SA[j]-1]=='$' || islower(idx->T[idx->SA[j]-1])){ //#has to be maximal
                return 1; //success
            }
        }
    }
    return 0;
}

PyObject * getmultimums(RevealIndex *index, PyObject *args, PyObject *keywds) {

    uint32_t minl=0;
    
    if (args!=NULL) {
        if (!PyArg_ParseTuple(args, "i", &minl))        
            return NULL;
    }
    
    PyObject * multimums;
    multimums=PyList_New(0);
    if (index==NULL){
        fprintf(stderr,"No valid index object.\n");
        return NULL;
    }
    RevealIndex * main = (RevealIndex *) index->main;
    int maxdepth=1000;
    int *flag_so=calloc(main->nsamples,sizeof *flag_so);
    uint32_t *stack_lcp=malloc(maxdepth * sizeof *stack_lcp);
    uint32_t *stack_lb=malloc(maxdepth * sizeof *stack_lb);
    uint32_t *stack_ub=malloc(maxdepth * sizeof *stack_ub);    
    int depth=0;
    saidx_t i,lb,i_lb,i_ub;
    uint32_t i_lcp;
    
    //int minl=10;

    stack_lcp[0]=0;
    stack_lb[0]=0;
    stack_ub[0]=0;
    for (i=1;i<index->n;i++){
        lb = i-1;
        assert(depth>=0);
        while (index->LCP[i] < stack_lcp[depth]) {
            stack_ub[depth]=i-1; //assign
            i_lcp = stack_lcp[depth];
            i_lb = stack_lb[depth];
            i_ub = stack_ub[depth];
            depth--;
            int n=(i_ub-i_lb)+1;
            
            if (i_lcp>minl){
                if (n<=main->nsamples){
                    if (ismultimum(index, i_lcp, i_lb, i_ub, flag_so)==1){
                        int x;
                        PyObject *sp=PyList_New(n);
                        for (x=0;x<n;x++) {
#ifdef SA64
                            PyObject *v = Py_BuildValue("L", index->SA[i_lb+x]);
#else
                            PyObject *v = Py_BuildValue("i", index->SA[i_lb+x]);
#endif
                            PyList_SET_ITEM(sp, x, v);
                        }

                        PyObject *multimum=Py_BuildValue("I,i,O",i_lcp,n,sp);
                        PyList_Append(multimums,multimum);
                        Py_DECREF(multimum);
                    }
                }
            }

            assert(depth>=0);
            lb = i_lb;
        }
        if (index->LCP[i] > stack_lcp[depth]){
            depth++;
            if (depth>=maxdepth){
                maxdepth=maxdepth+1000;
                stack_lcp=realloc(stack_lcp,maxdepth * sizeof *stack_lcp);
                if (stack_lcp==NULL){
                    fprintf(stderr,"Failed to allocate memory for stack_lcp.\n");
                    return NULL;
                }
                stack_lb=realloc(stack_lb,maxdepth * sizeof *stack_lb);
                if (stack_lb==NULL){
                    fprintf(stderr,"Failed to allocate memory for stack_lb.\n");
                    return NULL;
                }
                stack_ub=realloc(stack_ub,maxdepth * sizeof *stack_ub);
                if (stack_ub==NULL){
                    fprintf(stderr,"Failed to allocate memory for stack_ub.\n");
                    return NULL;
                }
            }
            stack_lcp[depth]=index->LCP[i];
            stack_lb[depth]=lb;
            stack_ub[depth]=0; //initialize
        }
    }

    while (depth>=0) {
        stack_ub[depth]=index->n-1;
        i_lcp = stack_lcp[depth];
        i_lb = stack_lb[depth];
        i_ub = stack_ub[depth];
        depth--;
        
        int n=(i_ub-i_lb)+1;
        if (n<=main->nsamples){
            if (ismultimum(index, i_lcp, i_lb, i_ub, flag_so)==1){
                int x;
                PyObject *sp=PyList_New(n);
                for (x=0;x<n;x++) {
#ifdef SA64
                    PyObject *v = Py_BuildValue("L", index->SA[i_lb+x]);
#else
                    PyObject *v = Py_BuildValue("i", index->SA[i_lb+x]);
#endif
                    PyList_SET_ITEM(sp, x, v);
                }

                PyObject *multimum=Py_BuildValue("I,i,O",i_lcp,n,sp);
                PyList_Append(multimums,multimum);
                Py_DECREF(multimum);
            }
        }
    }
    
    free(stack_lcp);
    free(stack_lb);
    free(stack_ub);
    free(flag_so);
    return multimums;
}

void split(RevealIndex *idx, uint8_t *D, RevealIndex *i_leading, RevealIndex *i_trailing, RevealIndex *i_par){
    saidx_t i=0,ip=0,il=0,it=0,lastp=0,lastl=0,lastt=0;
    uint32_t minlcpp=0,minlcpl=0,minlcpt=0;

    for (i=0; i<idx->n; i++){
        if (D[i]==1){ //write to leading
            assert(il<i_leading->n);
            i_leading->SA[il]=idx->SA[i];
            if (il==0){
                i_leading->LCP[il]=0;
            } else {
                i_leading->LCP[il]=minlcpl;
            }
            assert(il<idx->n);
            idx->SAi[idx->SA[i]]=il; //update inverse
            il++;
            lastl=i;
        } else if (D[i]==2){ //write to trailing
            assert(it<i_trailing->n);
            i_trailing->SA[it]=idx->SA[i];
            if (it==0){
                i_trailing->LCP[it]=0;
            } else {
                i_trailing->LCP[it]=minlcpt;
            }
            assert(it<idx->n);
            idx->SAi[idx->SA[i]]=it; //update inverse
            it++;
            lastt=i;
        } else {
            if (D[i]==3){ //suffixes that have been matched
                //fprintf(stderr,"MUM! %d\n",idx->SA[i]);
            } else{
                if (D[i]!=4){
                    //assert(idx->T[idx->SA[i]]=='$'); //can only happen after first alignment step
                    //fprintf(stderr,"D=%d i=%d n=%d\n",D[i],i,i_par->n);
                    continue;
                }
                
                assert(ip<i_par->n);
                i_par->SA[ip]=idx->SA[i];
                if (ip==0){
                    i_par->LCP[ip]=0;
                } else {
                    i_par->LCP[ip]=minlcpp;
                }
                assert(ip<idx->n);
                idx->SAi[idx->SA[i]]=ip; //update inverse
                ip++;
                lastp=i;
            }
        }

        if (i==idx->n-1){
            break;
        }

        if (i==lastt){
            minlcpt=idx->LCP[i+1];
        } else {
            if (idx->LCP[i+1]<minlcpt){
                minlcpt=idx->LCP[i+1];
            }
        }
        
        if (i==lastl){
            minlcpl=idx->LCP[i+1];
        } else {
            if (idx->LCP[i+1]<minlcpl){
                minlcpl=idx->LCP[i+1];
            }
        }
        
        if (i==lastp){
            minlcpp=idx->LCP[i+1];
        } else {
            if (idx->LCP[i+1]<minlcpp){
                minlcpp=idx->LCP[i+1];
            }
        }
    }
}

void bubble_sort(RevealIndex* idx, saidx_t* sp, int n, uint32_t l){
    int j=0;
    saidx_t i=0,x,tmpSA;
    uint32_t tmpLCP;
    
    for (j=0; j<n; j++) {
        
        for (i=0; i<idx->n; i++){
            if (idx->SA[i]<sp[j] && idx->SA[i]+idx->LCP[i]>sp[j]){ // if match overlaps the start position of sp[j]

                x=i;
                tmpSA=idx->SA[i];
                tmpLCP=idx->LCP[i];
                
                while ((idx->LCP[x] > sp[j]-tmpSA) && (x>0)){
                    assert(x<idx->n);
                    idx->SAi[idx->SA[x-1]]=x;
                    idx->SA[x]=idx->SA[x-1];
                    idx->LCP[x]=idx->LCP[x-1]; //!
                    x--;
                }

		assert(x<idx->n);
		idx->SAi[tmpSA]=x;

		idx->SA[x]=tmpSA;
		idx->LCP[x+1]=sp[j]-tmpSA;
                
		if (i<idx->n-1){
		    if (tmpLCP < idx->LCP[i+1]){
			idx->LCP[i+1]=tmpLCP; //!
		    }
		}
            } else {
                if (i<idx->n-1){
                    if (idx->SA[i]<sp[j] && idx->SA[i]+idx->LCP[i+1]>sp[j]){
                        if (idx->LCP[i+1] > idx->LCP[i]){
                            idx->LCP[i+1]= sp[j]-idx->SA[i];
                        }
                    }
                }
            }
        }
    }
    
    for (i=0; i<n; i++){
        for (j=sp[i];j<sp[i]+l;j++){
            idx->T[j]=tolower(idx->T[j]); //mark corresponding intervals in T
        }
    }
}

/* Alignment Thread */
void *aligner(void *arg) {

    RevealWorker *rw = arg;
    PyGILState_STATE gstate;
    RevealIndex * idx;
    while(1) {
	int hasindex=0;
        int i=0;
        
        pthread_mutex_lock(&mutex);/* acquire the mutex lock */
        aw++;
        idx=pop_index();
        if (idx==NULL) {
            hasindex=0;
            aw--;
        } else {
            hasindex=1;
        }
        pthread_mutex_unlock(&mutex);/* release the mutex lock */
        
        if (die==1 || (rw->threadid==-1 && hasindex==0)){
            break;
        }
        
        if (hasindex==1){
            //fprintf(stderr,"Starting alignment cycle (%d)...\n", rw->threadid);
            //fprintf(stderr,"samples=%d\n",idx->nsamples);
            //fprintf(stderr,"depth=%d\n",idx->depth);
            //fprintf(stderr,"n=%d\n",idx->n);
            
            assert(idx->nsamples>0);
            
            RevealMultiMUM mmum;
            mmum.sp=(saidx_t *) malloc(idx->nsamples*sizeof(saidx_t));
            mmum.l=0;
            mmum.score=0;
            mmum.penalty=0;
            
            PyObject *result;
            
            pthread_mutex_lock(&python);
            gstate = PyGILState_Ensure();

            if (rw->mumpicker!=Py_None){
                if (!PyCallable_Check(rw->mumpicker)) {
                    PyErr_SetString(PyExc_TypeError, "**** mumpicker isn't callable");
                    err_flag=1;
                    Py_DECREF(idx);
                    PyGILState_Release(gstate);
                    pthread_mutex_unlock(&python);

                    free(mmum.sp);
                    break;
                }
                
                PyObject *multimums;
                if (((RevealIndex *) idx->main)->nsamples>2){
                    multimums = getmultimums(idx,NULL,NULL);
                } else {
                    multimums = getmums(idx,NULL,NULL);
                }
                 
                PyObject *arglist = Py_BuildValue("(O,O)", multimums, idx);

                PyObject *mum = PyEval_CallObject(rw->mumpicker, arglist); //mumpicker returns intervals
                
                if (mum==Py_None){
                    //TODO 1: NO MORE MUMS, call prune nodes here!
                    Py_DECREF(idx);
                    Py_DECREF(arglist);
                    Py_DECREF(multimums);
                    Py_DECREF(mum);
                    PyGILState_Release(gstate);
                    pthread_mutex_unlock(&python);

                    pthread_mutex_lock(&mutex);
                    aw--;
                    pthread_mutex_unlock(&mutex);

                    free(mmum.sp);
                    continue;
                }
                
                if (!PyTuple_Check(mum)){
                    PyErr_SetString(PyExc_TypeError, "**** call to mumpicker failed");
                    err_flag=1;
                    Py_DECREF(idx);
                    Py_DECREF(arglist);
                    Py_DECREF(multimums);
                    PyGILState_Release(gstate);
                    pthread_mutex_unlock(&python);

                    free(mmum.sp);
                    break;
                }
                
                PyObject *sp=NULL;
                PyObject *tidx=NULL;

                PyArg_ParseTuple(mum,"IOiLOK", &mmum.l, &tidx, &mmum.n, &mmum.score, &sp, &mmum.penalty);
                
                for (i=0; i<mmum.n; i++){
                    PyObject * pos=PyList_GetItem(sp,i);
                    if (pos==NULL){
                        fprintf(stderr,"**** invalid results from mumpicker\n");
                        mmum.sp[i]=0;
                        continue;
                    }
#ifdef SA64
                    mmum.sp[i]=PyLong_AsLongLong(pos); //TODO: check what do to with this..
#else
                    mmum.sp[i]=PyInt_AS_LONG(pos);
#endif
                }

                result = PyEval_CallObject(rw->graphalign, mum);
                Py_DECREF(arglist);
                Py_DECREF(mum);
                Py_DECREF(multimums);
            } else {

                if (idx->nsamples>2){
                    getbestmultimum(idx, &mmum, idx->nsamples);
                } else { //simpler method (tiny bit more efficient), but essentially should produce the same results
                    //getlongestmum(idx, &mmum);
                    getbestmum(idx, &mmum, rw->wpen, rw->wscore);
                }
                
                if (mmum.l==0){
                    //TODO 2: NO MORE MUMS, call prune nodes here!

                    Py_DECREF(idx);
                    PyGILState_Release(gstate);
                    pthread_mutex_unlock(&python);

                    pthread_mutex_lock(&mutex);
                    aw--;
                    pthread_mutex_unlock(&mutex);

                    free(mmum.sp);
                    continue;                
                }
                
                PyObject *sps;
                PyObject *arglist;
                
                //build a list of start positions in the index
                sps=PyList_New(mmum.n);
                PyObject *pos;
                for (i=0;i<mmum.n;i++){
#ifdef SA64
                    pos=Py_BuildValue("L",mmum.sp[i]);
#else
                    pos=Py_BuildValue("i",mmum.sp[i]);
#endif
                    PyList_SetItem(sps,i,pos);
                }
                
                arglist = Py_BuildValue("(I,O,i,L,O,K)", mmum.l, idx, mmum.n, mmum.score, sps, mmum.penalty);
                
                if (arglist==NULL){
                    PyErr_SetString(PyExc_TypeError, "***** Building argument list failed");
                    err_flag=1;
                    Py_DECREF(idx);
                    PyGILState_Release(gstate);
                    pthread_mutex_unlock(&python);

                    free(mmum.sp);
                    break;
                }
                
                if (!PyCallable_Check(rw->graphalign)) {
                    PyErr_SetString(PyExc_TypeError, "***** graphalign isn't callable");
                    err_flag=1;
                    PyGILState_Release(gstate);
                    pthread_mutex_unlock(&python);
                    free(mmum.sp);
                    break;
                }
                
                result = PyEval_CallObject(rw->graphalign,arglist);
                Py_DECREF(arglist);
            }
            
            if (result==NULL){
                fprintf(stderr,"***** Python callback raised exception, stopping all threads!\n");
                err_flag=1;
                Py_DECREF(idx);
                PyGILState_Release(gstate);
                pthread_mutex_unlock(&python);
                free(mmum.sp);
                break;
            }
            
            PyObject *leading_intervals;
            PyObject *trailing_intervals;
            PyObject *rest;
            PyObject *merged;
            PyObject *newleft;
            PyObject *newright;
            
            if (result==Py_None){
                //TODO 3: NO MORE MUMS, call prune nodes here!

                Py_DECREF(idx);
                PyGILState_Release(gstate);
                pthread_mutex_unlock(&python);

                pthread_mutex_lock(&mutex);
                aw--;
                pthread_mutex_unlock(&mutex);

                free(mmum.sp);
                continue;                
            }
            
            if (!PyArg_ParseTuple(result, "OOOOOO", &leading_intervals, &trailing_intervals, &rest, &merged, &newleft, &newright)){
                fprintf(stderr,"Failed to parse result of call to graph_align!\n");
                //no tuple returned by python call, apparently we're done...
                Py_DECREF(idx);
                PyGILState_Release(gstate);
                pthread_mutex_unlock(&python);

                pthread_mutex_lock(&mutex);
                aw--;
                pthread_mutex_unlock(&mutex);
                free(mmum.sp);
                continue;
            }

            uint8_t *D=calloc(idx->n,sizeof(uint8_t));
            int *flag_so=calloc( ((RevealIndex *) idx->main)->nsamples,sizeof(int));
            
            saidx_t i,j,begin,end,trailingn=0,leadingn=0,parn=0;

            int nintv_leading=0, nintv_trailing=0, nintv_par=0;
            int leadingsamples=0, trailingsamples=0, parsamples=0;
            
            nintv_leading=PyList_Size(leading_intervals);
            for (i=0;i<nintv_leading;i++){
                PyObject *tup;
                tup=PyList_GetItem(leading_intervals,i);
#ifdef SA64
                PyArg_ParseTuple(tup,"LL",&begin,&end);
#else
                PyArg_ParseTuple(tup,"ii",&begin,&end);
#endif
                for (j=begin; j<end; j++){
                    D[idx->SAi[j]]=1; //leading  ************
                    leadingn++;
                }
                if (((RevealIndex *) idx->main)->nsamples > 2){
                    if (flag_so[idx->SO[begin]]==0){
                        flag_so[idx->SO[begin]]=1;
                        leadingsamples++;
                    }
                } else {
                    if (begin < idx->nsep[0] && flag_so[0]==0){
                        flag_so[0]=1;
                        leadingsamples++;
                    }
                    if (begin > idx->nsep[0] && flag_so[1]==0){
                        flag_so[1]=1;
                        leadingsamples++;
                    }
                }
            }
            memset(flag_so,0,((RevealIndex *) idx->main)->nsamples * sizeof(int));

            nintv_trailing=PyList_Size(trailing_intervals);
            for (i=0;i<nintv_trailing;i++){
                PyObject *tup;
                tup=PyList_GetItem(trailing_intervals,i);
#ifdef SA64
                PyArg_ParseTuple(tup,"LL",&begin,&end);
#else
                PyArg_ParseTuple(tup,"ii",&begin,&end);
#endif
                for (j=begin; j<end; j++){
                    D[idx->SAi[j]]=2; //trailing   *****
                    trailingn++;
                }
                if (((RevealIndex *) idx->main)->nsamples > 2){
                    if (flag_so[idx->SO[begin]]==0){
                        flag_so[idx->SO[begin]]=1;
                        trailingsamples++;
                    }
                } else {
                    if (begin < idx->nsep[0] && flag_so[0]==0){
                        flag_so[0]=1;
                        trailingsamples++;
                    }
                    if (begin > idx->nsep[0] && flag_so[1]==0){
                        flag_so[1]=1;
                        trailingsamples++;
                    }
                }

            }
            memset(flag_so,0,((RevealIndex *) idx->main)->nsamples * sizeof(int));
            
            nintv_par=PyList_Size(rest);
            for (i=0;i<nintv_par;i++){
                PyObject *tup;
                tup=PyList_GetItem(rest,i);
#ifdef SA64
                PyArg_ParseTuple(tup,"LL",&begin,&end);
#else
                PyArg_ParseTuple(tup,"ii",&begin,&end);
#endif
                
                for (j=begin; j<end; j++){
                    D[idx->SAi[j]]=4; //parallel paths  **********
                    parn++;
                }

                if (((RevealIndex *) idx->main)->nsamples > 2){
                    if (flag_so[idx->SO[begin]]==0){
                        flag_so[idx->SO[begin]]=1;
                        parsamples++;
                    }
                } else {
                    if (begin < idx->nsep[0] && flag_so[0]==0){
                        flag_so[0]=1;
                        parsamples++;
                    }
                    if (begin > idx->nsep[0] && flag_so[1]==0){
                        flag_so[1]=1;
                        parsamples++;
                    }
                }
            }

            free(flag_so);
            
            for (i=0;i<mmum.n;i++){
                for (j=mmum.sp[i];j<mmum.sp[i]+mmum.l;j++){
                    D[idx->SAi[j]]=3; //matching     **************
                }
            }
            
            /*fprintf(stderr,"Trailingsamples %d\n",trailingsamples);
            fprintf(stderr,"Leadingsamples %d\n",leadingsamples);
            fprintf(stderr,"Trailingn %lld\n",trailingn);
            fprintf(stderr,"Leadingn %lld\n",leadingn);            
            fprintf(stderr,"Parsamples %d\n",parsamples);
            fprintf(stderr,"mmum l %d\n",mmum.l);
            fprintf(stderr,"mmum n %d\n",mmum.n);
            fprintf(stderr,"Index n %lld\n",idx->n);*/
            
            //assert(parn==idx->n-(trailingn+leadingn+(mmum.l*mmum.n)));
            
            int newdepth=idx->depth+1; //update depth in recursion tree
            
            assert(newdepth>0);
            
            assert(leadingn>=0);
            Py_INCREF(leading_intervals);
            
            //fprintf(stderr,"Allocating leading (%zd nodes) %lld\n", PyList_Size(leading_intervals), leadingn);
            
            RevealIndex *i_leading=newIndex();
            
            i_leading->SA=malloc(leadingn*sizeof(saidx_t));
            i_leading->LCP=malloc(leadingn*sizeof(uint32_t));
            i_leading->nodes=leading_intervals;
            i_leading->depth=newdepth;
            i_leading->n=leadingn;
            i_leading->SAi=idx->SAi;
            i_leading->T=idx->T;
            i_leading->SO=idx->SO;
            i_leading->nsamples=leadingsamples;
            i_leading->nsep=idx->nsep;
            i_leading->main=idx->main;
            Py_INCREF(idx->left);
            i_leading->left=idx->left; //interval that is bounding on the left
            //Py_INCREF(merged);
            //i_leading->right=merged; //interval that is bounding on the right
            Py_INCREF(newright);
            i_leading->right=newright; //interval that is bounding on the right
            
            assert(trailingn>=0);
            Py_INCREF(trailing_intervals);
            //fprintf(stderr,"Allocating trailing (%zd nodes) %llu\n", PyList_Size(trailing_intervals), trailingn);
            RevealIndex *i_trailing=newIndex();

            i_trailing->SA=malloc(trailingn*sizeof(saidx_t));
            i_trailing->LCP=malloc(trailingn*sizeof(uint32_t));
            i_trailing->nodes=trailing_intervals;
            i_trailing->depth=newdepth;
            i_trailing->n=trailingn;
            i_trailing->SAi=idx->SAi;
            i_trailing->T=idx->T;
            i_trailing->SO=idx->SO;
            i_trailing->nsamples=trailingsamples;
            i_trailing->nsep=idx->nsep;
            i_trailing->main=idx->main;
            //Py_INCREF(merged);
            //i_trailing->left=merged; //interval that is bounding on the left
            Py_INCREF(newleft);
            i_trailing->left=newleft; //interval that is bounding on the left
            Py_INCREF(idx->right);
            i_trailing->right=idx->right; //interval that is bounding on the right

            assert(parn>=0);
            Py_INCREF(rest);
            //fprintf(stderr,"Allocating parallel (%zd nodes) %llu %d %d %llu\n", PyList_Size(rest), parn, mmum.l, mmum.n, idx->n);
            RevealIndex *i_parallel=newIndex();
            
            i_parallel->SA=malloc(parn*sizeof(saidx_t));
            i_parallel->LCP=malloc(parn*sizeof(uint32_t));
            i_parallel->nodes=rest;
            i_parallel->depth=newdepth;
            i_parallel->n=parn;
            i_parallel->SAi=idx->SAi;
            i_parallel->T=idx->T;
            i_parallel->SO=idx->SO;
            i_parallel->nsamples=parsamples;//idx->nsamples-(mmum.n);
            i_parallel->nsep=idx->nsep;
            i_parallel->main=idx->main;
            Py_INCREF(idx->left);
            i_parallel->left=idx->left; //interval that is bounding on the left
            Py_INCREF(idx->right);
            i_parallel->right=idx->right; //interval that is bounding on the right

            PyGILState_Release(gstate);
            pthread_mutex_unlock(&python);
            
	    //fprintf(stderr,"splitting...\n");
            split(idx, D, i_leading, i_trailing, i_parallel);
	    //fprintf(stderr,"done.\n");

            /*for (i=0; i<trailingn; i++){
                assert(i_trailing->LCP[i]<i_trailing->n);
            }

            for (i=0; i<leadingn; i++){
                assert(i_leading->LCP[i]<i_leading->n);
            }*/            
 
            bubble_sort(i_leading, mmum.sp, mmum.n, mmum.l);

            free(D);
            free(mmum.sp);
            if (idx->depth==0){
                free(idx->SA);
                idx->SA=NULL;
                free(idx->LCP);
                idx->LCP=NULL;
            }
            
            pthread_mutex_lock(&python); 
            gstate=PyGILState_Ensure();
            
            Py_DECREF(result);
            Py_DECREF(idx); //trigger gc for subindex
            
            pthread_mutex_lock(&mutex);
            nmums++;
            //add resulting indices to the queue 
            if (i_leading->n>0){
		     //fprintf(stderr,"push leading\n");
                if (!(push_index(i_leading)==0)){
                    fprintf(stderr,"Failed to add leading index to queue.\n");
                    Py_DECREF(i_leading);
                    err_flag=1;
                    pthread_mutex_unlock(&mutex);
                    break;
                }
            } else {
                Py_DECREF(i_leading);
            }

            if (i_trailing->n>0){
		     //fprintf(stderr,"push trailing\n");
                if (!(push_index(i_trailing)==0)){
                    fprintf(stderr,"Failed to add trailing index to queue.\n");
                    Py_DECREF(i_trailing);
                    err_flag=1;
                    pthread_mutex_unlock(&mutex);
                    break;
                }
            } else {
                Py_DECREF(i_trailing);
            }
            
            if (i_parallel->n>0){
    		     //fprintf(stderr,"push parallel\n");
                if (!(push_index(i_parallel)==0)){
                    fprintf(stderr,"Failed to add parallel paths index to queue.\n");
                    Py_DECREF(i_parallel);
                    err_flag=1;
                    pthread_mutex_unlock(&mutex);
                    break;
                }
            } else {
                Py_DECREF(i_parallel);
            }
            
            PyGILState_Release(gstate);
            pthread_mutex_unlock(&python);
            aw--;
            pthread_mutex_unlock(&mutex);
        }
        else {
            usleep(1);
        }
    }
    //fprintf(stderr,"Stopping alignment thread %d.\n",rw->threadid);
    free(rw);
    return NULL;
}

