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
    //fprintf(stderr,"qsize=%d qstart=%d\n",qsize,qstart);
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

int getlongestmum(RevealIndex *index, RevealMultiMUM *mum){
    int i=0,aStart,bStart;
    int lb,la;
    mum->u=1;
    mum->l=0;
    mum->score=0;
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
            mum->sp[0]=aStart;
            mum->sp[1]=bStart;
        }
    }
    return 0;
}

int getbestmum(RevealIndex *index, RevealMultiMUM *mum){
    int i=0,aStart,bStart;
    int lb,la;
    int penalize;
    if (PyList_Size(index->nodes)==2){
        penalize=1;
    } else {
        penalize=0;
    }
    mum->u=1;
    mum->l=0;
    mum->score=-2147483648;
    mum->n=2;
    for (i=1;i<index->n;i++){
        if (index->LCP[i]>mum->score){
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
            int penalty=0;
            int lpenalty=0;
            int tpenalty=0;
            int start1=-1,end1=-1,start2=-1,end2=-1;

            if (index->depth>0 && penalize==1){
                PyObject *node1=PyList_GetItem(index->nodes,0);
                start1=PyInt_AS_LONG(PyTuple_GetItem(node1,0));
                end1=PyInt_AS_LONG(PyTuple_GetItem(node1,1));
    
                PyObject *node2=PyList_GetItem(index->nodes,1);
                start2=PyInt_AS_LONG(PyTuple_GetItem(node2,0));
                end2=PyInt_AS_LONG(PyTuple_GetItem(node2,1));
                
                if (start1<=aStart && end1>=aStart){
                    lpenalty=abs((aStart-start1)-(bStart-start2)); //leading penalty
                    tpenalty=abs((end1-(aStart+index->LCP[i]))-(end2-(bStart+index->LCP[i]))); //trailing penalty
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
                    lpenalty=abs((bStart-start1)-(aStart-start2)); //leading penalty
                    tpenalty=abs((end1-(bStart+index->LCP[i]))-(end2-(aStart+index->LCP[i]))); //trailing penalty
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

	    //fprintf(stderr,"l=%d p=%d score=%d start1=%d aStart=%d start2=%d bStart=%d lpenalty=%d, tpenalty=%d\n",index->LCP[i],penalty,mum->score, start1, aStart, start2, bStart, lpenalty, tpenalty);
            if (index->LCP[i]-penalty > mum->score){
                mum->score=index->LCP[i]-penalty;
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
    mum->n=0;
    int i=0,j=0,la,lb;
    int *flag_so=calloc( index->nsamples,sizeof(int));
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

//PyObject * getmultimums(RevealIndex *index){
//    
//    PyObject * multimums;
//    multimums=PyList_New(0);
//    
//    if (index==NULL){
//        fprintf(stderr,"No valid index object.\n");
//        return NULL;
//    }
//    
////    int trace=0;
////    if (index->n==1515){
////        trace=1;    
////    }
//    
//    RevealIndex * main = (RevealIndex *) index->main;
//    int i=0,j=0,la,lb;
//    int *flag_so=calloc(main->nsamples,sizeof(int));
//    int flag_maximal=0;
//    int flag_unique=0;
//    int flag_lcp;
//    int x;
//    int min_n=2;
//    
//    for (i=min_n-1;i<index->n;i++){
//        if (i==index->n-1 || index->LCP[i+1]<index->LCP[i]){
//            int so;
//            flag_maximal=0;
//            flag_unique=0;
//            flag_lcp=index->LCP[i];
//            int n=0;
//
//            for (j=0;j<(main->nsamples-1) && i-j-1>=0;j++){ //scan back n places in SA
//                if (j==0){
//                    assert(index->SO[index->SA[i]]>=0 && index->SO[index->SA[i]]<main->nsamples);
//                    flag_so[index->SO[index->SA[i]]]=1;
//                    n=1;
//                }
//                assert(i-j-1>=0);
//                so=index->SO[index->SA[i-j-1]];
//                assert(so>=0 && so<=main->nsamples);
//                if (flag_so[so]==0){
//                    flag_so[so]=1;
//                } else {
//                    j--;
//                    break;
//                }
//                if (!flag_maximal){
//                    if (index->SA[i-j]>0 && index->SA[i-j-1]>0){
//                        assert((index->SA[i-j]-1)>=0);
//                        assert((index->SA[i-j-1]-1)>=0);
//                        if (index->T[index->SA[i-j]-1] != index->T[index->SA[i-j-1]-1] || islower(index->T[index->SA[i-j]-1]) || index->T[index->SA[i-j]-1]=='$'){
//                            flag_maximal=1;
//                        }
//                    } else {
//                        flag_maximal=1;
//                    }
//                }
//                
//                if (index->LCP[i-j]<flag_lcp){
//                    break;
//                    //found another multimum!
//                    
//                    assert((i-j)>=0 && (i-j)<index->n);
//                    flag_lcp=index->LCP[i-j];
//                }
//                n++;
//            }
//            
////            if (trace)
////                fprintf(stderr,"%d lcp=%d maximal=%d so[0]=%d so[1]=%d, n=%d min_n=%d\n",i,flag_lcp,flag_maximal,flag_so[0], flag_so[1], n, min_n);
//            
//            if (n>=min_n && flag_maximal){
//                
//                //is it unique within the entire index?
//                if (i==index->n-1) { //is it the last domain in the array, then only check predecessor
//                    assert(i-(n-1)>=0);
//                    lb=index->LCP[i-(n-1)];
//                    la=0;
//                } else if ((i-(n-1))==0) { //is it the first domain in the array, then only check successor
//                    lb=0;
//                    assert(i+1<index->n);
//                    la=index->LCP[i+1]; //USE THIS TO SCAN FASTER! LCP[i+1] has to be smaller than LCP[i] otherwise not optimal interval!
//                } else {
//                    assert(i-(n-1)>=0);
//                    lb=index->LCP[i-(n-1)];
//                    assert(i+1<index->n);
//                    la=index->LCP[i+1];
//                }
//                
////                if (trace)
////                    fprintf(stderr,"%d lb=%d la=%d\n",i,lb,la);
//                
//                if (flag_lcp>lb && flag_lcp>la){
//                    flag_unique=1;
//
//                    PyObject *sp=PyList_New(n);
//                    for (x=0;x<n;x++) {
//                        assert(i-x>=0);
//                        PyObject *v = Py_BuildValue("i", index->SA[i-x]);
//                        PyList_SET_ITEM(sp, x, v);
//                    }
//                    PyObject *multimum=Py_BuildValue("i,i,O",flag_lcp,n,sp);
//                    PyList_Append(multimums,multimum);
//                }
//            }
//            for (j=0;j<main->nsamples;j++){ //memset?
//                flag_so[j]=0; //reset so flags
//            }
//            flag_maximal=0; //reset maximal flag
//        }
//    }
//    free(flag_so);
//
//    return multimums;
//}


int ismultimum(RevealIndex * idx, int l, int lb, int ub, int * flag_so) {
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

PyObject * getmultimums(RevealIndex *index) {
    PyObject * multimums;
    multimums=PyList_New(0);
    if (index==NULL){
        fprintf(stderr,"No valid index object.\n");
        return NULL;
    }
    RevealIndex * main = (RevealIndex *) index->main;
    int maxdepth=1000;
    int *flag_so=calloc(main->nsamples,sizeof *flag_so);
    int *stack_lcp=malloc(maxdepth * sizeof *stack_lcp);
    int *stack_lb=malloc(maxdepth * sizeof *stack_lb);
    int *stack_ub=malloc(maxdepth * sizeof *stack_ub);    
    int depth=0;
    int i,lb,i_lb,i_ub,i_lcp;
    stack_lcp[0]=0;
    stack_lb[0]=0;
    stack_ub[0]=0;
    for (i=1;i<index->n;i++){
        lb = i-1;
        assert(depth>=0);
        while (index->LCP[i] < stack_lcp[depth]) {
            stack_ub[depth]=i-1;
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
                        PyObject *v = Py_BuildValue("i", index->SA[i_lb+x]);
                        PyList_SET_ITEM(sp, x, v);
                    }
                    PyObject *multimum=Py_BuildValue("i,i,O",i_lcp,n,sp);
                    PyList_Append(multimums,multimum);
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
                stack_lb=realloc(stack_lb,maxdepth * sizeof *stack_lb);
                stack_ub=realloc(stack_ub,maxdepth * sizeof *stack_ub);
            }
            stack_lcp[depth]=index->LCP[i];
            stack_lb[depth]=lb;
            stack_ub[depth]=0;
        }
    }
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
                PyObject *v = Py_BuildValue("i", index->SA[i_lb+x]);
                PyList_SET_ITEM(sp, x, v);
            }
            PyObject *multimum=Py_BuildValue("i,i,O",i_lcp,n,sp);
            PyList_Append(multimums,multimum);
        }
    }
    free(stack_lcp);
    free(stack_lb);
    free(stack_ub);
    free(flag_so);
    return multimums;
}


void split(RevealIndex *idx, uint8_t *D, RevealIndex *i_leading, RevealIndex *i_trailing, RevealIndex *i_par){
    int i=0,ip=0,il=0,it=0,minlcpp=0,minlcpl=0,minlcpt=0,lastp=0,lastl=0,lastt=0;
    
    for (i=0; i<idx->n; i++){
        if (D[i]==1){ //write to leading
            i_leading->SA[il]=idx->SA[i];
            if (il==0){
                i_leading->LCP[il]=0;
            } else {
                i_leading->LCP[il]=minlcpl;
            }
            idx->SAi[idx->SA[i]]=il; //update inverse
            il++;
            lastl=i;
        } else if (D[i]==2){ //write to trailing
            i_trailing->SA[it]=idx->SA[i];
            if (it==0){
                i_trailing->LCP[it]=0;
            } else {
                i_trailing->LCP[it]=minlcpt;
            }
            idx->SAi[idx->SA[i]]=it; //update inverse
            it++;
            lastt=i;
        } else {
            if (D[i]==3){ //suffixes that have been matched
                //fprintf(stderr,"MUM! %d\n",idx->SA[i]);
            } else{
                assert(D[i]==0); //parallel paths
                
                if (ip>=i_par->n)
                    fprintf(stderr,"ip=%d n=%d\n",ip,i_par->n);
                
                assert(ip<i_par->n);
                i_par->SA[ip]=idx->SA[i];
                if (ip==0){
                    i_par->LCP[ip]=0;
                } else {
                    i_par->LCP[ip]=minlcpp;
                }
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

void bubble_sort(RevealIndex* idx, int* sp, int n, int l){
    int i=0,j=0,x,tmpSA,tmpLCP;
    for (i=0; i<n; i++){
        for (j=sp[i];j<sp[i]+l;j++){
            idx->T[j]=tolower(idx->T[j]);
        }
    }
    for (j=0; j<n; j++) {
        for (i=0; i<idx->n; i++){
            if (idx->SA[i]<sp[j] && idx->SA[i]+idx->LCP[i]>sp[j]){ //
                x=i;
                tmpSA=idx->SA[i];
                tmpLCP=idx->LCP[i];
                while ( (idx->LCP[x] > sp[j]-tmpSA) && (x>0) ){
                    idx->SAi[idx->SA[x-1]]=x;
                    idx->SA[x]=idx->SA[x-1];
                    idx->LCP[x]=idx->LCP[x-1]; //!
                    x--;
                }
    			idx->SAi[tmpSA]=x;
	    		idx->SA[x]=tmpSA;
		    	idx->LCP[x+1]=sp[j]-tmpSA;
			    if ((tmpLCP < idx->LCP[i+1]) && (i < idx->n-1)){
				    idx->LCP[i+1]=tmpLCP; //!
			    }
                continue;
            }
            if (i<idx->n-1){ 
                if (idx->SA[i]<sp[j] && idx->SA[i]+idx->LCP[i+1]>sp[j]){
                    if (idx->LCP[i+1] > idx->LCP[i]){
                        idx->LCP[i+1]= sp[j]-idx->SA[i]; 
                        continue;
                    }
                }
            }
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
        idx=pop_index();
        if (idx==NULL){
            hasindex=0;
        } else {
            hasindex=1;
            aw++;
        }
        pthread_mutex_unlock(&mutex);/* release the mutex lock */
        
        if (die==1 || (rw->threadid==-1 && hasindex==0)){
	    //if (hasindex==0) fprintf(stderr,"no more indices\n");
	    //if (die==1) fprintf(stderr,"die\n");
            break;
        }

        if (hasindex==1){
            RevealMultiMUM mmum;
            mmum.sp=(int *) malloc(idx->nsamples*sizeof(int));
            mmum.l=0;
            mmum.score=0;
            
            PyObject *result;
            
            if (rw->mumpicker!=Py_None){
                pthread_mutex_lock(&python);      //LOCK PYTHON            
                gstate = PyGILState_Ensure();

                if (!PyCallable_Check(rw->mumpicker)) {
                    PyErr_SetString(PyExc_TypeError, "**** mumpicker isn't callable");
                    err_flag=1;
                    Py_DECREF(idx);
                    PyGILState_Release(gstate);
                    pthread_mutex_unlock(&python);          //UNLOCK PYTHON
                    break;
                }
                
                PyObject *multimums = getmultimums(idx);
                
                PyObject *arglist = Py_BuildValue("(O,O)", multimums, idx);
                
                PyObject *mum = PyEval_CallObject(rw->mumpicker, arglist); //mumpicker calls align and returns domains
                
                if (mum==Py_None){
                    Py_DECREF(idx);
                    Py_DECREF(arglist);
                    Py_DECREF(multimums);
                    Py_DECREF(mum);
                    PyGILState_Release(gstate);
                    pthread_mutex_unlock(&python);      //UNLOCK PYTHON
                    pthread_mutex_lock(&mutex);
                    aw--;
                    pthread_mutex_unlock(&mutex);
                    continue;
                }
                
                if (!PyTuple_Check(mum)){
                    PyErr_SetString(PyExc_TypeError, "**** call to mumpicker failed");
                    err_flag=1;
                    Py_DECREF(idx);
                    Py_DECREF(arglist);
                    Py_DECREF(multimums);
                    PyGILState_Release(gstate);
                    pthread_mutex_unlock(&python);          //UNLOCK PYTHON
                    break;
                }
                                
                PyObject *sp=NULL;
                PyObject *tmp=NULL;
                
                PyArg_ParseTuple(mum,"iOiiO",&mmum.l,&tmp,&mmum.n,&mmum.score,&sp);
                
                for (i=0; i<mmum.n; i++){
                    PyObject * pos=PyList_GetItem(sp,i);
                    if (pos==NULL){
                        PyErr_SetString(PyExc_TypeError, "**** invalid results from mumpicker");
                        err_flag=1;
                        Py_DECREF(idx);
                        PyGILState_Release(gstate);
                        pthread_mutex_unlock(&python);          //UNLOCK PYTHON
                        break;
                    }
                    mmum.sp[i]=PyInt_AS_LONG(pos);
                }
                
                result = PyEval_CallObject(rw->graphalign, mum); //mumpicker calls align and returns domains
                
                Py_DECREF(arglist);
                Py_DECREF(mum);
                Py_DECREF(multimums);
            } else {
                
                if (idx->nsamples>2){
                    getbestmultimum(idx, &mmum, idx->nsamples);
                } else { //simpler method (tiny bit more efficient), but essentially should produce the same results
                    //getlongestmum(idx, &mmum);
                    getbestmum(idx, &mmum);
                }
                
                pthread_mutex_lock(&python);      //LOCK PYTHON            
                gstate = PyGILState_Ensure();
                
                PyObject *sps;
                PyObject *arglist;
                
                sps=PyList_New(mmum.n);
                PyObject *pos;
                for (i=0;i<mmum.n;i++){
                    pos=Py_BuildValue("i",mmum.sp[i]);
                    PyList_SetItem(sps,i,pos);
                }
                
                arglist = Py_BuildValue("(i,O,i,i,O)", mmum.l, idx, mmum.n, mmum.score, sps);
                if (arglist==NULL){
                    PyErr_SetString(PyExc_TypeError, "***** Building argument list failed");
                    err_flag=1;
                    Py_DECREF(idx);
                    PyGILState_Release(gstate);
                    pthread_mutex_unlock(&python);      //UNLOCK PYTHON
                    break;
                }
                
                if (!PyCallable_Check(rw->graphalign)) {
                    PyErr_SetString(PyExc_TypeError, "***** graphalign isn't callable");
                    err_flag=1;
                    PyGILState_Release(gstate);
                    pthread_mutex_unlock(&python);          //UNLOCK PYTHON
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
                pthread_mutex_unlock(&python);      //UNLOCK PYTHON
                break;
            }
            
            PyObject *leading_intervals;
            PyObject *trailing_intervals;
            PyObject *rest;
            //TODO: pycheck result is tuple
            
            if (result==Py_None){
                Py_DECREF(idx);
                PyGILState_Release(gstate);
                pthread_mutex_unlock(&python);      //UNLOCK PYTHON
                pthread_mutex_lock(&mutex);
                aw--;
                pthread_mutex_unlock(&mutex);
                continue;                
            }
            
            if (!PyArg_ParseTuple(result, "OOO", &leading_intervals, &trailing_intervals, &rest)){
                fprintf(stderr,"Failed to parse result of call to graph_align, expecting a tuple with two lists!\n");
                //no tuple returned by python call, apparently we're done...
                Py_DECREF(idx);
                PyGILState_Release(gstate);
                pthread_mutex_unlock(&python);      //UNLOCK PYTHON
                
                pthread_mutex_lock(&mutex);
                aw--;
                pthread_mutex_unlock(&mutex);
                continue;
            }
            
            uint8_t *D=calloc(idx->n,sizeof(uint8_t));
            
            int i,j,begin,end,trailingn=0,leadingn=0,parn=0;
            for (i=0;i<PyList_Size(leading_intervals);i++){
                PyObject *tup;
                tup=PyList_GetItem(leading_intervals,i);
                PyArg_ParseTuple(tup,"ii",&begin,&end);
                for (j=begin; j<end; j++){
                    D[idx->SAi[j]]=1;
                    leadingn++;
                }
            }
            
            for (i=0;i<PyList_Size(trailing_intervals);i++){
                PyObject *tup;
                tup=PyList_GetItem(trailing_intervals,i);
                PyArg_ParseTuple(tup,"ii",&begin,&end);
                for (j=begin; j<end; j++){
                    D[idx->SAi[j]]=2;
                    trailingn++;
                }
            }
            
//            for (i=0;i<PyList_Size(rest);i++){
//                PyObject *tup;
//                tup=PyList_GetItem(rest,i);
//                PyArg_ParseTuple(tup,"ii",&begin,&end);
//                for (j=begin; j<end; j++){
//                    D[idx->SAi[j]]=0;
//                    parn++;
//                }
//            }
            
            for (i=0;i<mmum.n;i++){
                for (j=mmum.sp[i];j<mmum.sp[i]+mmum.l;j++){
                    D[idx->SAi[j]]=3;
                }
            }
//            fprintf(stderr,"Trailingn %d\n",trailingn);
//            fprintf(stderr,"Leadingn %d\n",leadingn);
//            fprintf(stderr,"mmum l %d\n",mmum.l);
//            fprintf(stderr,"mmum n %d\n",mmum.n);
//            fprintf(stderr,"Index n %d\n",idx->n);
            parn=idx->n-(trailingn+leadingn+(mmum.l*mmum.n));
            //assert that parn equals length of all intervals in rest
//            fprintf(stderr,"parn %d\n",parn);

            int newdepth=idx->depth+1; //update depth in recursion tree
            
            assert(newdepth>0);
            
            Py_INCREF(leading_intervals);
            //fprintf(stderr,"Allocating leading (%zd nodes) %d\n", PyList_Size(leading_intervals), leadingn);
            RevealIndex *i_leading=newIndex();
            i_leading->SA=malloc(leadingn*sizeof(int));
            i_leading->LCP=malloc(leadingn*sizeof(int));
            i_leading->nodes=leading_intervals;
            i_leading->depth=newdepth;
            i_leading->n=leadingn;
            i_leading->SAi=idx->SAi;
            i_leading->T=idx->T;
            i_leading->SO=idx->SO;
            i_leading->nsamples=mmum.n;
            i_leading->nsep=idx->nsep;
            i_leading->main=idx->main;
            
            Py_INCREF(trailing_intervals);
            //fprintf(stderr,"Allocating trailing (%zd nodes) %d\n", PyList_Size(trailing_intervals), trailingn);
            RevealIndex *i_trailing=newIndex();
            i_trailing->SA=malloc(trailingn*sizeof(int));
            i_trailing->LCP=malloc(trailingn*sizeof(int));
            i_trailing->nodes=trailing_intervals;
            i_trailing->depth=newdepth;
            i_trailing->n=trailingn;
            i_trailing->SAi=idx->SAi;
            i_trailing->T=idx->T;
            i_trailing->SO=idx->SO;
            i_trailing->nsamples=mmum.n;
            i_trailing->nsep=idx->nsep;
            i_trailing->main=idx->main;
            
            Py_INCREF(rest);
            //fprintf(stderr,"Allocating parallel (%zd nodes) %d %d %d %d\n", PyList_Size(rest), parn, mmum.l, mmum.n, idx->n);
            RevealIndex *i_parallel=newIndex();
            i_parallel->SA=malloc(parn*sizeof(int));
            i_parallel->LCP=malloc(parn*sizeof(int));
            i_parallel->nodes=rest;
            i_parallel->depth=newdepth;
            i_parallel->n=parn;
            i_parallel->SAi=idx->SAi;
            i_parallel->T=idx->T;
            i_parallel->SO=idx->SO;
            i_parallel->nsamples=idx->nsamples;
            i_parallel->nsep=idx->nsep;
            i_parallel->main=idx->main;
            
            Py_DECREF(result);

            PyGILState_Release(gstate);
            pthread_mutex_unlock(&python);       //UNLOCK PYTHON
            
            split(idx, D, i_leading, i_trailing, i_parallel);
            
            bubble_sort(i_leading, mmum.sp, mmum.n, mmum.l);
            //bubble_sort(i_trailing, mmum.sp, mmum.n, mmum.l);
            //bubble_sort(i_parallel, mmum.sp, mmum.n, mmum.l);
            
            free(D);
            free(mmum.sp);

            //pthread_mutex_lock(&python);       //LOCK PYTHON
            //PyGILState_Ensure();
            Py_DECREF(idx); //trigger gc for subindex
            //PyGILState_Release(gstate);
            //pthread_mutex_unlock(&python);       //UNLOCK PYTHON

            pthread_mutex_lock(&mutex);
            nmums++;
            //add resulting indices to the queue 
            if (i_leading->n>0){
                if (!(push_index(i_leading)==0)){
                    fprintf(stderr,"Failed to add leading index to queue.\n");
                    Py_DECREF(i_leading);
                    err_flag=1;
                    pthread_mutex_unlock(&mutex);
                    break;
                }
            }
            if (i_trailing->n>0){
                if (!(push_index(i_trailing)==0)){
                    fprintf(stderr,"Failed to add trailing index to queue.\n");
                    Py_DECREF(i_trailing);
                    err_flag=1;
                    pthread_mutex_unlock(&mutex);
                    break;
                }
            }
            if (i_parallel->n>0){
                if (!(push_index(i_parallel)==0)){
                    fprintf(stderr,"Failed to add parallel paths index to queue.\n");
                    Py_DECREF(i_parallel);
                    err_flag=1;
                    pthread_mutex_unlock(&mutex);
                    break;
                }
            }
            aw--;
            pthread_mutex_unlock(&mutex);
        } else {
            sleep(1);
        }
    }
    //fprintf(stderr,"Stopping alignment thread %d.\n",rw->threadid);
    free(rw);
    return NULL;
}

