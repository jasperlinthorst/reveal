#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <pthread.h>
#include <semaphore.h>
#include <unistd.h>
#include <time.h>
#include "Python.h"
#include "reveal.h"
#include <math.h>

/* The mutex lock */
extern pthread_mutex_t mutex, python;
extern RevealIndex **index_queue;
extern int maxqsize,qsize,qstart,aw,nmums,err_flag,die;

/* pops an index of the queue */
RevealIndex* pop_index(void) {
    //fprintf(stderr,"qend=%d qstart=%d qsize=%d\n",qsize,qstart,qsize-qstart);
    //if(qsize > qstart) { FIFO
    if(qsize > 0) { //LIFO
#ifdef REVEALDEBUG
        fprintf(stderr,"POPPED: Number of indices left on queue: %d\n",qsize-1);
#endif
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

#ifdef REVEALDEBUG
    fprintf(stderr,"PUSHED: Number of indices on queue: %d\n",qsize);
#endif
    return 0;
}

PyObject * getmums(RevealIndex *index, PyObject *args, PyObject *keywds){
    int minl=0;
    lcp_t lb,la;

    if (args!=NULL) {
        if (!PyArg_ParseTuple(args, "i", &minl))
            return NULL;
    }
    
    saidx_t i=0,aStart,bStart;
    PyObject *mums=PyList_New(0);

    for (i=1;i<index->n;i++){
        if (index->LCP[i]<minl){
            continue;
        }
        if (((index->SA[i])>(index->nsep[0])) == ((index->SA[i-1])>(index->nsep[0]))) { //repeat
            continue;
        }
        if ((index->SA[i])<(index->SA[i-1])) {
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

        if (index->rc==1){
            bStart=index->nsep[0] + (index->nT - bStart - index->LCP[i]);
        }

#ifdef SA64
        PyObject *mum=Py_BuildValue("I,(L,L),i",index->LCP[i],aStart,bStart,index->rc);
#else
        PyObject *mum=Py_BuildValue("I,(i,i),i",index->LCP[i],aStart,bStart,index->rc);
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


PyObject * getmums_rem(RevealIndex *index, PyObject *args, PyObject *keywds){
    int minl=0;
    lcp_t lb,la;

    if (args!=NULL) {
        if (!PyArg_ParseTuple(args, "i", &minl))
            return NULL;
    }
    
    saidx_t i=0,aStart,bStart;
    PyObject *mums=PyList_New(0);

    for (i=1;i<index->n;i++){
        if (index->LCP[i]<minl){
            continue;
        }
        if (((index->SA[i])>(index->nsep[0])) == ((index->SA[i-1])>(index->nsep[0]))) { //repeat
            continue;
        }
        if ((index->SA[i])<(index->SA[i-1])) {
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

        if (index->rc==1){
            bStart=index->nsep[0] + (index->n - bStart - index->LCP[i]);
        }

#ifdef SA64
        PyObject *mum=Py_BuildValue("I,i,((i:L),(i:L))",index->LCP[i],2,0,aStart,1,bStart);
#else
        PyObject *mum=Py_BuildValue("I,i,((i:i),(i:i))",index->LCP[i],2,0,aStart,1,bStart);
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
    lcp_t lb,la;
    mum->l=0;
    mum->score=0;
    mum->penalty=0;
    mum->n=2;
    for (i=1;i<index->n;i++){
        if (index->LCP[i]>mum->l){
            if ((index->SA[i]>index->nsep[0]) == (index->SA[i-1]>index->nsep[0])) { //repeat
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

int ismultimum(RevealIndex * idx, lcp_t l, int lb, int ub, int * flag_so) {
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

int ismultimem(RevealIndex * idx, lcp_t l, int lb, int ub, int * flag_so) {
    if (l>0){
        int j;
        memset(flag_so,0,((RevealIndex *) idx->main)->nsamples * sizeof(int));
        
        if (((RevealIndex *) idx->main)->nsamples==2){ //dont need SO in case of only two samples
            flag_so[(idx->SA[ub]>idx->nsep[0]) == (idx->SA[lb]>idx->nsep[0])]++;
            // if ( (idx->SA[ub]>idx->nsep[0]) == (idx->SA[lb]>idx->nsep[0]) ){
            //     return 0;
            // }
        } else {
            for (j=lb; j<ub+1; j++) { //has to occur in all samples (but may exist more than once)
                flag_so[idx->SO[idx->SA[j]]]++;
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

PyObject * getmultimems(RevealIndex *index, PyObject *args, PyObject *keywds) {
    lcp_t minl=0;
    int minn=2;
    static char *kwlist[] = {"minlength","minn", NULL};
    
    if (args!=NULL) {
        if (!PyArg_ParseTupleAndKeywords(args, keywds, "|ii", kwlist, &minl, &minn))
            return NULL;
    }
    
    PyObject * multimems;
    multimems=PyList_New(0);
    if (index==NULL){
        fprintf(stderr,"No valid index object.\n");
        return NULL;
    }
    RevealIndex * mainidx = (RevealIndex *) index->main;
    int maxdepth=1000;
    int *flag_so=calloc(mainidx->nsamples,sizeof *flag_so);
    lcp_t *stack_lcp=malloc(maxdepth * sizeof *stack_lcp);
    lcp_t *stack_lb=malloc(maxdepth * sizeof *stack_lb);
    lcp_t *stack_ub=malloc(maxdepth * sizeof *stack_ub);
    lcp_t i_lcp;
    int depth=0;
    saidx_t i,lb,i_lb,i_ub;
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
            
            if (i_lcp>=minl){
                if (n>=minn){
                    if (ismultimem(index, i_lcp, i_lb, i_ub, flag_so)==1){
                        int c=0,ci=0;
                        for (ci=0; ci<((RevealIndex *) index->main)->nsamples; ci++){
                            if (flag_so[ci]>0){
                                c++;
                            }
                        }
                        if (c<minn){
                            continue;
                        }
                        int x;
                        PyObject *crdlst = PyTuple_New(n);
                        for (x=0;x<n;x++) {
#ifdef SA64
                            PyObject *v = Py_BuildValue("(i,L)", index->SO[index->SA[i_lb+x]], index->SA[i_lb+x]);
#else
                            PyObject *v = Py_BuildValue("(i,i)", index->SO[index->SA[i_lb+x]], index->SA[i_lb+x]);
#endif
                            PyTuple_SetItem(crdlst, x, v);
                        }
                        PyObject *multimem=Py_BuildValue("I,i,O",i_lcp,c,crdlst);
                        Py_DECREF(crdlst);
                        PyList_Append(multimems,multimem);
                        Py_DECREF(multimem);
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
        if (i_lcp>=minl){
            if (n>=minn){
                if (ismultimem(index, i_lcp, i_lb, i_ub, flag_so)==1){
                    int c=0,ci=0;
                    for (ci=0; ci<((RevealIndex *) index->main)->nsamples; ci++){
                        if (flag_so[ci]>0){
                            c++;
                        }
                    }
                    if (c<minn){
                        continue;
                    }
                    int x;
                    PyObject *crdlst = PyTuple_New(n);
                    for (x=0;x<n;x++) {
#ifdef SA64
                        PyObject *v = Py_BuildValue("(i,L)", index->SO[index->SA[i_lb+x]], index->SA[i_lb+x]);
#else
                        PyObject *v = Py_BuildValue("(i,i)", index->SO[index->SA[i_lb+x]], index->SA[i_lb+x]);
#endif
                        PyTuple_SetItem(crdlst, x, v);
                    }
                    PyObject *multimem=Py_BuildValue("I,i,O",i_lcp,c,crdlst);
                    Py_DECREF(crdlst);
                    PyList_Append(multimems,multimem);
                    Py_DECREF(multimem);
                }
            }
        }
    }
    free(stack_lcp);
    free(stack_lb);
    free(stack_ub);
    free(flag_so);
    return multimems;
}

PyObject * getmultimums(RevealIndex *index, PyObject *args, PyObject *keywds) {
    lcp_t minl=0;
    int minn=2;
    static char *kwlist[] = {"minlength","minn", NULL};
    
    if (args!=NULL) {
        if (!PyArg_ParseTupleAndKeywords(args, keywds, "|ii", kwlist, &minl, &minn))
            return NULL;
    }
    
    PyObject * multimums;
    multimums=PyList_New(0);
    if (index==NULL){
        fprintf(stderr,"No valid index object.\n");
        return NULL;
    }

    RevealIndex * mainidx = (RevealIndex *) index->main;
    int maxdepth=1000;
    int *flag_so=calloc(mainidx->nsamples,sizeof *flag_so);
    lcp_t *stack_lcp=malloc(maxdepth * sizeof *stack_lcp);
    lcp_t *stack_lb=malloc(maxdepth * sizeof *stack_lb);
    lcp_t *stack_ub=malloc(maxdepth * sizeof *stack_ub);
    lcp_t i_lcp;
    int depth=0;
    saidx_t i,lb,i_lb,i_ub;
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
            
            if (i_lcp>=minl){
                if (n<=mainidx->nsamples && n>=minn){
                    if (ismultimum(index, i_lcp, i_lb, i_ub, flag_so)==1){
                        int x;
                        // PyObject *crdmap = PyDict_New();
                        PyObject *crdmap = PyTuple_New(n);
                        for (x=0;x<n;x++) {
                            // PyObject *s = Py_BuildValue("i", index->SO[index->SA[i_lb+x]]);
#ifdef SA64
                            // PyObject *v = Py_BuildValue("L", index->SA[i_lb+x]);
                            PyObject *v = Py_BuildValue("(i,L)", index->SO[index->SA[i_lb+x]], index->SA[i_lb+x]);
#else
                            // PyObject *v = Py_BuildValue("i", index->SA[i_lb+x]);
                            PyObject *v = Py_BuildValue("(i,i)", index->SO[index->SA[i_lb+x]], index->SA[i_lb+x]);
#endif
                            PyTuple_SetItem(crdmap, x, v);

                            // PyDict_SetItem(crdmap, s, v);
                            // Py_DECREF(s);
                            // Py_DECREF(v);
                        }
                        PyObject *multimum=Py_BuildValue("I,i,O",i_lcp,n,crdmap);
                        Py_DECREF(crdmap);
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
        if (i_lcp>=minl){
            if (n<=mainidx->nsamples && n>=minn){
                if (ismultimum(index, i_lcp, i_lb, i_ub, flag_so)==1){
                    int x;
                    // PyObject *crdmap = PyDict_New();
                    PyObject *crdmap = PyTuple_New(n);
                    for (x=0;x<n;x++) {
                        // PyObject *s = Py_BuildValue("i", index->SO[index->SA[i_lb+x]]);
#ifdef SA64
                        // PyObject *v = Py_BuildValue("L", index->SA[i_lb+x]);
                        PyObject *v = Py_BuildValue("(i,L)", index->SO[index->SA[i_lb+x]], index->SA[i_lb+x]);
#else
                        // PyObject *v = Py_BuildValue("i", index->SA[i_lb+x]);
                        PyObject *v = Py_BuildValue("(i,i)", index->SO[index->SA[i_lb+x]], index->SA[i_lb+x]);
#endif
                        PyTuple_SetItem(crdmap, x, v);

                        // PyDict_SetItem(crdmap, s, v);
                        // Py_DECREF(s);
                        // Py_DECREF(v);
                    }
                    PyObject *multimum=Py_BuildValue("I,i,O",i_lcp,n,crdmap);
                    Py_DECREF(crdmap);
                    PyList_Append(multimums,multimum);
                    Py_DECREF(multimum);
                }
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
    
    lcp_t minlcpp=0,minlcpl=0,minlcpt=0;

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

void bubble_sort(RevealIndex* idx, PyObject* matching_intervals){

    lcp_t tmpLCP;
    saidx_t i=0,x,tmpSA,begin,end;
    PyObject *iter;
    PyObject *tup;

    iter=PyObject_GetIter(matching_intervals);
    while ((tup=PyIter_Next(iter))){
        
        #ifdef SA64
        PyArg_ParseTuple(tup,"LL",&begin,&end);
        #else
        PyArg_ParseTuple(tup,"ii",&begin,&end);
        #endif

        // fprintf(stderr,"BUBBLE SORT: %d-%d: (%.100s)-(%.100s)\n",begin,end,idx->T+begin,idx->T+end-100);

        for (i=0; i<idx->n; i++) { // for each suffix

            if ( (idx->SA[i] < begin) && ((idx->SA[i]+idx->LCP[i]) > begin) ){ // if match overlaps the start position
                x=i;
                tmpSA=idx->SA[i];
                tmpLCP=idx->LCP[i];
                
                while ((idx->LCP[x] >= begin-tmpSA) && (x>0)){
                    assert(x<idx->n);
                    idx->SAi[idx->SA[x-1]]=x;
                    idx->SA[x]=idx->SA[x-1];
                    idx->LCP[x]=idx->LCP[x-1]; //!
                    
                    x--;
                }
                assert(x<idx->n);
                idx->SAi[tmpSA]=x;
                idx->SA[x]=tmpSA;
                idx->LCP[x+1]=begin-tmpSA;
                        
                if (i<idx->n-1){ //if not last entry of LCP
                    if (tmpLCP < idx->LCP[i+1]){ 
                        idx->LCP[i+1]=tmpLCP;
                        //Check if T[SA[i+1]+LCP[i+1]-1] is not lower()
                    }
                }
            } else {

                if (i<idx->n-1){

                    if ((idx->SA[i] < begin) && ((idx->SA[i]+idx->LCP[i+1]) > begin ) ){
                        if (idx->LCP[i+1] > idx->LCP[i]) {
                            idx->LCP[i+1]=begin-idx->SA[i];
                        }
                    }

                }
            }
        }

        Py_DECREF(tup);
    }
    Py_DECREF(iter);
}


/* Alignment Thread */
void *aligner(void *arg) {

    #ifdef REVEALDEBUG
    time_t t0,t1;
    #endif

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
        
        if (hasindex==1) {

            #ifdef REVEALDEBUG
            // fprintf(stderr,"Initial.\n");
            //checkindex(idx);
            fprintf(stderr,"Starting alignment cycle (threadid=%d)\n", rw->threadid);
            fprintf(stderr,"samples=%d\n",idx->nsamples);
            fprintf(stderr,"depth=%d\n",idx->depth);
            fprintf(stderr,"n=%d\n",idx->n);
            fprintf(stderr,"minl=%d\n", rw->minl);
            // fprintf(stderr,"wpen=%d...\n", rw->wpen);
            // fprintf(stderr,"wscore=%d...\n", rw->wscore);
            #endif
            assert(idx->nsamples>0);

            RevealMultiMUM mmum;
            mmum.sp=(saidx_t *) malloc(idx->nsamples*sizeof(saidx_t));
            mmum.l=0;
            mmum.score=0;
            mmum.penalty=0;
            
            PyObject *result;
            
            pthread_mutex_lock(&python);
            gstate = PyGILState_Ensure();

            // if (rw->mumpicker!=Py_None){
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
            PyObject *mumobject;
            // PyObject *sp=NULL;
            PyObject *spd=NULL;            
            PyObject *skipmumsleft=NULL;
            PyObject *skipmumsright=NULL;
            PyObject *precomputed=NULL;

            if (PyList_Size(idx->skipmums)==0){
                Py_INCREF(Py_False);
                precomputed=Py_False;
                #ifdef REVEALDEBUG
                time(&t0);
                fprintf(stderr,"Extracting new mums... %d\n",rw->minn);
                #endif
                if (((RevealIndex *) idx->main)->nsamples>2){
                    PyObject *args = PyTuple_New(0);
                    PyObject *kwargs = Py_BuildValue("{s:i, s:i}", "minlength", rw->minl, "minn", rw->minn);
                    
                    multimums = getmultimums(idx,args,kwargs);
                    // multimums = getmultimems(idx,args,kwargs);

                    Py_DECREF(kwargs);
                    Py_DECREF(args);
                } else {
                    PyObject *args = Py_BuildValue("(i)", rw->minl);
                    multimums = getmums_rem(idx,args,NULL);
                    Py_DECREF(args);
                }

#ifdef REVEALDEBUG
                time(&t1);
                fprintf(stderr,"Done (took %.f seconds).\n",difftime(t1,t0));
#endif
                
            }
            else {
#ifdef REVEALDEBUG
                fprintf(stderr,"Using precomputed mum...\n");
#endif
                multimums=idx->skipmums;
                Py_INCREF(Py_True);
                precomputed=Py_True;
            }

            PyObject *keywds = Py_BuildValue("{s:O, s:i}", "precomputed", precomputed, "minlength", rw->minl);
            // PyObject *keywds = Py_BuildValue("{s:O}", "prevchain", multimums);
            Py_DECREF(precomputed);
            // Py_DECREF(idx->skipmums);
            
            PyObject *arglist = Py_BuildValue("(O,O)", multimums, idx);
            Py_DECREF(multimums);
            
#ifdef REVEALDEBUG
            time(&t0);
            fprintf(stderr,"Selecting best mum (python callback)...\n");
#endif
            PyObject *pickresult = PyEval_CallObjectWithKeywords(rw->mumpicker, arglist, keywds); //mumpicker returns intervals
            Py_DECREF(arglist);

#ifdef REVEALDEBUG
            time(&t1);
            fprintf(stderr,"Done (took %.f seconds).\n",difftime(t1,t0));
#endif
            
            if (!PyTuple_Check(pickresult)){
                PyErr_SetString(PyExc_TypeError, "**** call to mumpicker failed");
                err_flag=1;
                Py_DECREF(idx);
                PyGILState_Release(gstate);
                pthread_mutex_unlock(&python);

                free(mmum.sp);
                break;
            }

            if (PyTuple_Size(pickresult)==0){
                //TODO 1: NO MORE MUMS, call prune nodes here!
                Py_DECREF(idx);
                Py_DECREF(pickresult);

                PyGILState_Release(gstate);
                pthread_mutex_unlock(&python);

                pthread_mutex_lock(&mutex);
                aw--;
                pthread_mutex_unlock(&mutex);
                
                free(mmum.sp);
                continue;
            }
            
#ifdef REVEALDEBUG
            fprintf(stderr,"Parsing mum tuple...\n");
#endif
            PyArg_ParseTuple(pickresult,"OOO", &mumobject, &skipmumsleft, &skipmumsright);
            
            Py_INCREF(mumobject);
            Py_INCREF(skipmumsleft);
            Py_INCREF(skipmumsright);
            
            Py_DECREF(pickresult);

            if (!PyTuple_Check(mumobject)) {
                fprintf(stderr,"Invalid mum tuple...\n");
            }

            PyArg_ParseTuple(mumobject,"IiO", &mmum.l, &mmum.n, &spd);

#ifdef REVEALDEBUG
            fprintf(stderr,"Done.\n");
#endif

#ifdef REVEALDEBUG
            fprintf(stderr,"Convert PyList[%d] to c...\n",mmum.n);
#endif

            for (i=0; i<mmum.n; i++){
                PyObject * tup=PyTuple_GetItem(spd,i);
                PyObject * pos=PyTuple_GetItem(tup,1);

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

            // Py_DECREF(spd);

#ifdef REVEALDEBUG
            fprintf(stderr,"Done.\n");
#endif

#ifdef REVEALDEBUG
            fprintf(stderr,"Graphalign (python callback)...\n");
#endif

            PyObject *tmp =Py_BuildValue("(O,O)", idx, mumobject);

            result = PyEval_CallObject(rw->graphalign, tmp);

            Py_DECREF(tmp);
            Py_DECREF(mumobject);

            

#ifdef REVEALDEBUG
            fprintf(stderr,"Done.\n");
#endif

            // Py_DECREF(arglist);
            // Py_DECREF(pickresult);
            // Py_DECREF(multimums);
            
            PyObject *leading_intervals;
            PyObject *trailing_intervals;
            PyObject *matching_intervals;
            PyObject *rest;
            PyObject *merged;
            PyObject *newleftnode;
            PyObject *newrightnode;

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
            
            if (!PyTuple_Check(result)){
                fprintf(stderr,"**** call to graphalign failed\n");
                PyErr_SetString(PyExc_TypeError, "**** call to graphalign failed");
                err_flag=1;
                Py_DECREF(idx);
                PyGILState_Release(gstate);
                pthread_mutex_unlock(&python);
                free(mmum.sp);
                break;
            }

            if (!PyArg_ParseTuple(result, "OOOOOOO", &leading_intervals, &trailing_intervals, &matching_intervals, &rest, &merged, &newleftnode, &newrightnode)) {
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

#ifdef REVEALDEBUG
            fprintf(stderr,"Parsing done.\n");
#endif

            uint8_t *D=calloc(idx->n,sizeof(uint8_t));
            int *flag_so=calloc( ((RevealIndex *) idx->main)->nsamples,sizeof(int));
            
            saidx_t i,j,begin,end,trailingn=0,leadingn=0,parn=0;

            // int nintv_leading=0, nintv_trailing=0, nintv_par=0;
            int leadingsamples=0, trailingsamples=0, parsamples=0;
            
            PyObject *iter;
            PyObject *tup;

            // nintv_leading=PySet_Size(leading_intervals);
            iter=PyObject_GetIter(leading_intervals);
            while ((tup=PyIter_Next(iter))){
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
                Py_DECREF(tup);
            }
            Py_DECREF(iter);
            memset(flag_so,0,((RevealIndex *) idx->main)->nsamples * sizeof(int));

            // nintv_trailing=PySet_Size(trailing_intervals);
            iter=PyObject_GetIter(trailing_intervals);
            while ((tup=PyIter_Next(iter))){
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
                Py_DECREF(tup);
            }
            Py_DECREF(iter);
            memset(flag_so,0,((RevealIndex *) idx->main)->nsamples * sizeof(int));
            
            // nintv_par=PySet_Size(rest);
            iter=PyObject_GetIter(rest);
            while ((tup=PyIter_Next(iter))){
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
                Py_DECREF(tup);
            }
            Py_DECREF(iter);
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
            
            RevealIndex *i_leading=NULL;
            if (leadingn>0){
                //fprintf(stderr,"Allocating leading (%zd nodes) %lld\n", PyList_Size(leading_intervals), leadingn);
                i_leading=newIndex();
                i_leading->SA=malloc(leadingn*sizeof(saidx_t));
                i_leading->LCP=malloc(leadingn*sizeof(lcp_t));
                Py_INCREF(leading_intervals);
                i_leading->nodes=leading_intervals;
                i_leading->depth=newdepth;
                i_leading->n=leadingn;
                i_leading->SAi=idx->SAi;
                i_leading->T=idx->T;
                i_leading->SO=idx->SO;
                i_leading->nsamples=leadingsamples;
                i_leading->nsep=idx->nsep;
                i_leading->main=idx->main;
                Py_INCREF(idx->left_node);
                i_leading->left_node=idx->left_node; //interval that is bounding on the left
                Py_INCREF(newrightnode);
                i_leading->right_node=newrightnode; //interval that is bounding on the right
                Py_INCREF(skipmumsleft);
                i_leading->skipmums=skipmumsleft;
            }

            assert(trailingn>=0);
            RevealIndex *i_trailing=NULL;
            if (trailingn>0){
                //fprintf(stderr,"Allocating trailing (%zd nodes) %llu\n", PyList_Size(trailing_intervals), trailingn);
                i_trailing=newIndex();
                i_trailing->SA=malloc(trailingn*sizeof(saidx_t));
                i_trailing->LCP=malloc(trailingn*sizeof(lcp_t));
                Py_INCREF(trailing_intervals);
                i_trailing->nodes=trailing_intervals;
                i_trailing->depth=newdepth;
                i_trailing->n=trailingn;
                i_trailing->SAi=idx->SAi;
                i_trailing->T=idx->T;
                i_trailing->SO=idx->SO;
                i_trailing->nsamples=trailingsamples;
                i_trailing->nsep=idx->nsep;
                i_trailing->main=idx->main;
                Py_INCREF(newleftnode);
                i_trailing->left_node=newleftnode; //interval that is bounding on the left
                Py_INCREF(idx->right_node);
                i_trailing->right_node=idx->right_node; //interval that is bounding on the right
                Py_INCREF(skipmumsright);
                i_trailing->skipmums=skipmumsright;
            }

            RevealIndex *i_parallel=NULL;
            assert(parn>=0);
            if (parn>0){
                //fprintf(stderr,"Allocating parallel (%zd nodes) %llu %d %d %llu\n", PyList_Size(rest), parn, mmum.l, mmum.n, idx->n);
                i_parallel=newIndex();
                i_parallel->SA=malloc(parn*sizeof(saidx_t));
                i_parallel->LCP=malloc(parn*sizeof(lcp_t));
                Py_INCREF(rest);
                i_parallel->nodes=rest;
                i_parallel->depth=newdepth;
                i_parallel->n=parn;
                i_parallel->SAi=idx->SAi;
                i_parallel->T=idx->T;
                i_parallel->SO=idx->SO;
                i_parallel->nsamples=parsamples;//idx->nsamples-(mmum.n);
                i_parallel->nsep=idx->nsep;
                i_parallel->main=idx->main;
                Py_INCREF(idx->left_node);
                i_parallel->left_node=idx->left_node; //interval that is bounding on the left
                Py_INCREF(idx->right_node);
                i_parallel->right_node=idx->right_node; //interval that is bounding on the right
                i_parallel->skipmums=PyList_New(0);
            }

            PyGILState_Release(gstate);
            pthread_mutex_unlock(&python);
            
#ifdef REVEALDEBUG
            time(&t0);
            fprintf(stderr,"Splitting SA... ");
#endif

            split(idx, D, i_leading, i_trailing, i_parallel);

#ifdef REVEALDEBUG
            time(&t1);
            fprintf(stderr,"Done (took %.f seconds).\n",difftime(t1,t0));
#endif

#ifdef REVEALDEBUG
            time(&t0);
            fprintf(stderr,"Marking intervals in T...");
#endif
            
            //mark corresponding intervals in T
            for (j=0; j<mmum.n; j++){
                for (i=mmum.sp[j];i<mmum.sp[j]+mmum.l;i++){
                    idx->T[i]=tolower(idx->T[i]);
                }
            }

#ifdef REVEALDEBUG
            time(&t0);
            fprintf(stderr,"done.\n");
#endif


#ifdef REVEALDEBUG
            time(&t0);
            fprintf(stderr,"Bubble sorting leading SA...");
#endif

            pthread_mutex_lock(&python);
            gstate=PyGILState_Ensure();

            if (leadingn>0){
                bubble_sort(i_leading, matching_intervals);
            }

#ifdef REVEALDEBUG
            time(&t0);
            fprintf(stderr,"done.\n");
#endif

#ifdef REVEALDEBUG
            time(&t1);
            fprintf(stderr,"Done (took %.f seconds).\n",difftime(t1,t0));
            if (trailingn>0) {
                fprintf(stderr,"Check trailing.\n");
                checkindex(i_trailing);
            }
            if (parn>0) {
                fprintf(stderr,"Check parallel.\n");
                checkindex(i_parallel);
            }
            if (leadingn>0) {
                fprintf(stderr,"Check leading.\n");
                checkindex(i_leading);
            }
#endif

            free(D);
            free(mmum.sp);

            if (idx->depth==0){
                free(idx->SA);
                idx->SA=NULL;
                free(idx->LCP);
                idx->LCP=NULL;
            }
            
            // pthread_mutex_lock(&python);
            // gstate=PyGILState_Ensure();
            
            Py_DECREF(result);
            Py_DECREF(idx); //trigger gc for subindex
            
            pthread_mutex_lock(&mutex);
            nmums++;

            //add resulting indices to the queue 
            if (parn>0){
                if (!(push_index(i_parallel)==0)){
                    fprintf(stderr,"Failed to add parallel paths index to queue.\n");
                    Py_DECREF(i_parallel);
                    err_flag=1;
                    pthread_mutex_unlock(&mutex);
                    break;
                }
            }
            
            if (leadingn>0){
                if (!(push_index(i_leading)==0)){
                    fprintf(stderr,"Failed to add leading index to queue.\n");
                    Py_DECREF(i_leading);
                    err_flag=1;
                    pthread_mutex_unlock(&mutex);
                    break;
                }
            }
            
            if (trailingn>0){
                if (!(push_index(i_trailing)==0)){
                    fprintf(stderr,"Failed to add trailing index to queue.\n");
                    Py_DECREF(i_trailing);
                    err_flag=1;
                    pthread_mutex_unlock(&mutex);
                    break;
                }
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






void checkindex(RevealIndex* idx){
    saidx_t i=0;
    int l=0, j=0;

    for (i=0; i<idx->n; i++) {
        l=idx->LCP[i];
        assert(l>=0);
        if (l==0){
            continue;
        }

        j=l-1;

        //for (j=0; j<l; j++){

            if (!(idx->T[idx->SA[i]+j]<=90 && idx->T[idx->SA[i]+j]>64)){
                #ifdef SA64
                fprintf(stderr,"i=%lld; l=%d j=%d --> %c %c %c\n",(long long)i,l,j,idx->T[idx->SA[i]+j-1],idx->T[idx->SA[i]+j],idx->T[idx->SA[i]+j+1]);
                #else
                fprintf(stderr,"i=%d; l=%d j=%d --> %c %c %c\n",i,l,j,idx->T[idx->SA[i]+j-1], idx->T[idx->SA[i]+j], idx->T[idx->SA[i]+j+1]);

                fprintf(stderr,"SA[%d]=%d %d %.100s\n",i-1,idx->SA[i-1],idx->LCP[i-1],idx->T+(idx->SA[i-1]));
                fprintf(stderr,"SA[%d]=%d %d %.100s\n",i,idx->SA[i],idx->LCP[i],idx->T+(idx->SA[i]));
                fprintf(stderr,"SA[%d]=%d %d %.100s\n",i+1,idx->SA[i+1],idx->LCP[i+1],idx->T+(idx->SA[i+1]));
                // fprintf(stderr,"%d\n",idx->SO[idx->SA[i]+j]-1);
                // fprintf(stderr,"%d\n",idx->SO[idx->SA[i]+j]);
                // fprintf(stderr,"%d\n",idx->SO[idx->SA[i]+j]+1);

                #endif
            }

            // i=4179; l=5 j=4 --> t a t

            assert(idx->T[idx->SA[i]+j]<=90); //check it wasn't matched
            assert(idx->T[idx->SA[i]+j]>64); //check it does not contain sentinel
        //}
    }
}


PyObject* extract(RevealIndex* idx, PyObject* args, PyObject *keywds) {
    PyObject *intervals;
    PyArg_ParseTuple(args,"O",&intervals);
    uint8_t *D=calloc(idx->n,sizeof(uint8_t));

    PyObject *iter;
    PyObject *tup;
    saidx_t i=0,j=0,begin,end,_begin,_end;

    saidx_t matching=0;

    iter=PyObject_GetIter(intervals);

    int x=0;
    while ((tup=PyIter_Next(iter))){

        #ifdef SA64
            PyArg_ParseTuple(tup,"LL",&begin,&end);
        #else
            PyArg_ParseTuple(tup,"ii",&begin,&end);
        #endif

        // fprintf(stderr,"before remap %d-%d: (%.100s)-(%.100s)\n",begin,end,idx->T+begin,idx->T+end-100);
        // fprintf(stderr,"nsep=%d n=%d nT=%d\n",idx->nsep[0],idx->n, idx->nT);

        if (idx->rc==1 && begin>idx->nsep[0]) { //map qry coordinates back to correct intervals
            _begin= idx->nsep[0]+(idx->nT - begin - (end-begin));
            _end=   idx->nsep[0]+(idx->nT - begin);
            begin=_begin;
            end=_end;
            assert(begin<end);

            #ifdef SA64
                PyObject *tup_=Py_BuildValue("(L,L)",begin,end);
            #else
                PyObject *tup_=Py_BuildValue("(i,i)",begin,end);
            #endif

            PyList_SetItem(intervals,x,tup_);
        }

        // fprintf(stderr,"before marking %d-%d: (%.100s)-(%.100s)\n",begin,end,idx->T+begin,idx->T+end-100);

        for (j=begin; j<end; j++){
            D[idx->SAi[j]]=3; //matching  ************
            idx->T[j]=tolower(idx->T[j]); //mark suffixes in T
            matching++;
        }

        // fprintf(stderr,"after marking %d-%d: (%.100s)-(%.100s)\n",begin,end,idx->T+begin,idx->T+end-100);

        x++;
        Py_DECREF(tup);
    }

    Py_DECREF(iter);

    // fprintf(stderr,"Allocate new SA and LCP...%d %d\n",idx->n,matching);
    saidx_t *_SA=malloc((idx->n-matching)*sizeof(saidx_t));    
    lcp_t *_LCP=malloc((idx->n-matching)*sizeof(lcp_t));
    
    lcp_t minlcp=0;
    j=1;

    // fprintf(stderr,"Mark matching suffixes... %d %d\n",idx->n,matching);

    _LCP[0]=0;
    for (i=1; i<idx->n; i++){
        assert(j<=i);
        if (D[i]!=3){ //not a matching suffix, add to the new SA
            _SA[j]=idx->SA[i];
            idx->SAi[_SA[j]]=j;

            if (D[i-1]==3){

                if (minlcp<idx->LCP[i]){
                    _LCP[j]=minlcp;
                } else {
                    _LCP[j]=idx->LCP[i];
                }
                
            } else {
                _LCP[j]=idx->LCP[i];
            }

            j++;
        } else {
            if (D[i-1]!=3){ //first match suffix
                minlcp=idx->LCP[i];
            } else {
                if (idx->LCP[i]<minlcp) {
                    minlcp=idx->LCP[i];
                }
            }
        }
    }

    // fprintf(stderr,"Free up old SA and LCP...\n");
    if (idx->SA!=NULL){
        free(idx->SA);
    }
    
    if (idx->LCP!=NULL){
        free(idx->LCP);
    }

    idx->SA=_SA;
    idx->LCP=_LCP;
    idx->n=idx->n-matching;

    // bubble_sort(idx, tmpIntervals);
    bubble_sort(idx, intervals);

    // fprintf(stderr,"Checkindex.\n");
    // checkindex(idx);
    // fprintf(stderr,"Done.\n");

    Py_INCREF(Py_None);
    return Py_None;
}









PyObject* splitindex(RevealIndex* idx, PyObject* args, PyObject *keywds) {
    
    PyObject *leading_intervals;
    PyObject *trailing_intervals;
    PyObject *matching_intervals;
    PyObject *rest;
    PyObject *merged;
    PyObject *newleftnode;
    PyObject *newrightnode;
    PyObject *skipmumsleft;
    PyObject *skipmumsright;

    PyArg_ParseTuple(args,"OOOOOOOOO",&leading_intervals,&trailing_intervals,&matching_intervals,&rest,&merged,&newleftnode,&newrightnode,&skipmumsleft,&skipmumsright);

    uint8_t *D=calloc(idx->n,sizeof(uint8_t));
    
    saidx_t j,begin,end,trailingn=0,leadingn=0,parn=0;

    int leadingsamples=0, trailingsamples=0, parsamples=0;

    PyObject *iter;
    PyObject *tup;

    int *flag_so=calloc( ((RevealIndex *) idx->main)->nsamples,sizeof(int));

    iter=PyObject_GetIter(leading_intervals);
    while ((tup=PyIter_Next(iter))){
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

        Py_DECREF(tup);
    }
    Py_DECREF(iter);

    memset(flag_so,0,((RevealIndex *) idx->main)->nsamples * sizeof(int));
    iter=PyObject_GetIter(trailing_intervals);
    while ((tup=PyIter_Next(iter))){
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
        Py_DECREF(tup);
    }
    Py_DECREF(iter);

    memset(flag_so,0,((RevealIndex *) idx->main)->nsamples * sizeof(int));
    iter=PyObject_GetIter(matching_intervals);
    while ((tup=PyIter_Next(iter))){
    #ifdef SA64
        PyArg_ParseTuple(tup,"LL",&begin,&end);
    #else
        PyArg_ParseTuple(tup,"ii",&begin,&end);
    #endif
        for (j=begin; j<end; j++){
            D[idx->SAi[j]]=3; //matching  **********
            idx->T[j]=tolower(idx->T[j]); //mark suffixes in T
        }
        Py_DECREF(tup);
    }
    Py_DECREF(iter);
    
    memset(flag_so,0,((RevealIndex *) idx->main)->nsamples * sizeof(int));
    iter=PyObject_GetIter(rest);
    while ((tup=PyIter_Next(iter))){
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
        Py_DECREF(tup);
    }
    Py_DECREF(iter);

    free(flag_so);

    int newdepth=idx->depth+1; //update depth in recursion tree

    assert(newdepth>0);
    assert(leadingn>=0);

    RevealIndex *i_leading=NULL;
    if (leadingn>0){
        i_leading=newIndex();
        i_leading->SA=malloc(leadingn*sizeof(saidx_t));
        i_leading->LCP=malloc(leadingn*sizeof(lcp_t));
        Py_INCREF(leading_intervals);
        i_leading->nodes=leading_intervals;
        i_leading->depth=newdepth;
        i_leading->n=leadingn;
        i_leading->SAi=idx->SAi;
        i_leading->T=idx->T;
        i_leading->SO=idx->SO;
        i_leading->nsamples=leadingsamples;
        i_leading->nsep=idx->nsep;
        i_leading->main=idx->main;
        Py_INCREF(idx->left_node);
        i_leading->left_node=idx->left_node; //interval that is bounding on the left
        Py_INCREF(newrightnode);
        i_leading->right_node=newrightnode; //interval that is bounding on the right
        Py_INCREF(skipmumsleft);
        i_leading->skipmums=skipmumsleft;
    } else{
        Py_INCREF(Py_None);
        i_leading=(RevealIndex *) Py_None;
    }

    assert(trailingn>=0);
    RevealIndex *i_trailing=NULL;
    if (trailingn>0){
        //fprintf(stderr,"Allocating trailing (%zd nodes) %llu\n", PyList_Size(trailing_intervals), trailingn);
        i_trailing=newIndex();
        i_trailing->SA=malloc(trailingn*sizeof(saidx_t));
        i_trailing->LCP=malloc(trailingn*sizeof(lcp_t));
        Py_INCREF(trailing_intervals);
        i_trailing->nodes=trailing_intervals;
        i_trailing->depth=newdepth;
        i_trailing->n=trailingn;
        i_trailing->SAi=idx->SAi;
        i_trailing->T=idx->T;
        i_trailing->SO=idx->SO;
        i_trailing->nsamples=trailingsamples;
        i_trailing->nsep=idx->nsep;
        i_trailing->main=idx->main;
        Py_INCREF(newleftnode);
        i_trailing->left_node=newleftnode; //interval that is bounding on the left
        Py_INCREF(idx->right_node);
        i_trailing->right_node=idx->right_node; //interval that is bounding on the right
        Py_INCREF(skipmumsright);
        i_trailing->skipmums=skipmumsright;
    } else{
        Py_INCREF(Py_None);
        i_trailing=(RevealIndex *) Py_None;
    }

    RevealIndex *i_parallel=NULL;
    assert(parn>=0);
    if (parn>0){
        //fprintf(stderr,"Allocating parallel (%zd nodes) %llu %d %d %llu\n", PyList_Size(rest), parn, mmum.l, mmum.n, idx->n);
        i_parallel=newIndex();
        i_parallel->SA=malloc(parn*sizeof(saidx_t));
        i_parallel->LCP=malloc(parn*sizeof(lcp_t));
        Py_INCREF(rest);
        i_parallel->nodes=rest;
        i_parallel->depth=newdepth;
        i_parallel->n=parn;
        i_parallel->SAi=idx->SAi;
        i_parallel->T=idx->T;
        i_parallel->SO=idx->SO;
        i_parallel->nsamples=parsamples;//idx->nsamples-(mmum.n);
        i_parallel->nsep=idx->nsep;
        i_parallel->main=idx->main;
        Py_INCREF(idx->left_node);
        i_parallel->left_node=idx->left_node; //interval that is bounding on the left
        Py_INCREF(idx->right_node);
        i_parallel->right_node=idx->right_node; //interval that is bounding on the right
        i_parallel->skipmums=PyList_New(0);
    } else{
        Py_INCREF(Py_None);
        i_parallel=(RevealIndex *) Py_None;
    }

    split(idx, D, i_leading, i_trailing, i_parallel);

    if (leadingn>0){
        bubble_sort(i_leading, matching_intervals);
    }

    tup=Py_BuildValue("(O,O,O)",i_leading,i_trailing,i_parallel);

    return tup;
}


