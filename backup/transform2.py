import reveallib
import reveallib64
from utils import *
from multiprocessing.pool import Pool
import signal
try:
    from matplotlib import pyplot as plt
    from matplotlib import patches as patches
except:
    pass

import os
import math


def transform(args):
    logging.debug("Extracting mums.")

    if args.output==None:
        pref=[]
        for f in [os.path.basename(args.reference),os.path.basename(args.contigs)]:
            bn=os.path.basename(f)
            if '.' in bn:
                pref.append(bn[:bn.find('.')])
            else:
                pref.append(bn)
        args.output="_".join(pref)

    mums=dict()
    blocks=dict()

    dseq={}
    do={}

    o=0
    for name,seq in fasta_reader(args.reference):
        dseq[name]=seq
        do[name]=o
        o+=len(seq)
        mums[name]=dict()
        blocks[name]=dict()

    o=0
    for name,seq in fasta_reader(args.contigs):
        assert(name not in dseq) #name for contigs in query and reference should be unique!
        do[name]=o
        o+=len(seq)
        dseq[name]=seq
        for ref in mums:
            mums[ref][name]=[[],[]]
            blocks[ref][name]=[[],[]]

    os.system("mummer -qthreads 3 -F -c -l 20 -b -n %s %s -mum > mums"%(args.reference,args.contigs))

    dl=[]
    with open("mums") as mumfile:
        for line in mumfile:
            if line.startswith(">"):
                line=line.rstrip()
                if line.endswith("Reverse"):
                    qry=line.split()[1]
                    revcomp=1
                else:
                    qry=line.split()[1]
                    revcomp=0
            else:
                v=line.strip().split()
                mums[v[0]][qry][revcomp].append((int(v[1]),int(v[2]),int(v[3])))
                dl.append(int(v[3]))

    maxdist=10000
    dd=[]

    for ref in mums:
        ro=do[ref]
        for qry in mums[ref]:
            qo=do[qry]
            for revcomp in range(2):
                t=mums[ref][qry][revcomp]
                if len(t)==0:
                    continue

                t.sort(key=lambda l: l[0]) #sort by reference position
                
                # tree=kdtree(t, 2)
                # print tree

                block=[t[0]]
                for i,mum in enumerate(t[1:]):
                    
                    assert(mum[0]>t[i][0])

                    d=abs( mum[0] - (t[i][0]+t[i][2]) )
                    dd.append(d)

                    if d>maxdist: #end of a block
                        blocks[ref][qry][revcomp].append(block)
                        block=[mum]
                    else:
                        block.append(mum)

                blocks[ref][qry][revcomp].append(block)

                print ref,qry,revcomp,len(t),len(blocks[ref][qry][revcomp])
    
    # plt.hist(dd,bins=100)
    # plt.show()

    # plt.clf()

    colstr="bgrcmk"

    # for ref in blocks:
    #     # ro=do[ref]
    #     ro=0
    #     for qry in blocks[ref]:
            
    #         plt.clf()

    #         # qo=do[qry]
    #         qo=0
    #         for revcomp in range(2):
    #             print ref,qry

    #             for i,block in enumerate(blocks[ref][qry][revcomp]):
                    
    #                 print "Block consists of",len(block),"anchors"
    #                 print block
    #                 c=colstr[i%(len(colstr)-1)]

    #                 for mum in block:
    #                     if revcomp==0:
    #                         plt.plot((ro+mum[0],ro+mum[0]+mum[2]),(qo+mum[1],qo+mum[1]+mum[2]),'%s-'%c)
    #                     else:
    #                         plt.plot((ro+mum[0],ro+mum[0]+mum[2]),(qo+mum[1]+mum[2],qo+mum[1]),'%s-'%c)

    #                 # plt.plot( (ro+block[0][0],ro+block[-1][0]) , (qo+block[0][1], qo+block[-1][1]) ,'b-')

    #         plt.show()

    d=[]

    for ref in blocks:
        # ro=do[ref]
        ro=0
        for qry in blocks[ref]:
            
            # plt.clf()

            # qo=do[qry]
            qo=0
            for revcomp in range(2):
                print ref,qry

                for i,block in enumerate(blocks[ref][qry][revcomp]):
                    

                    chain(block)


                    # print "Block consists of",len(block),"anchors"
                    # print block
                    c=colstr[i%(len(colstr)-1)]

                    for mum in block:
                        if revcomp==0:

                            x=(ro+mum[0],ro+mum[0]+mum[2])
                            y=(qo+mum[1],qo+mum[1]+mum[2])

                            d.append(x[0]-y[0])
                            dl.append(mum[2])
                            
                            # rx=int(x[0]*math.cos(math.radians(360-45))-y[0]*math.sin(math.radians(360-45)))
                            # ry=int(x[0]*math.sin(math.radians(360-45))+y[0]*math.cos(math.radians(360-45)))
                            
                            x=(x[0]+y[0],x[1]+y[1])
                            y=(x[0]-y[0],x[1]-y[1])

                            # print "after rotation",x,y

                            plt.plot(x,y,'r-')

                # plt.ylim(-5000000,5000000)
                # plt.xlim(0,5000000)
                plt.show()

    plt.clf()
    plt.hist(d,bins=100)
    plt.show()

    print len(mums), len(blocks)


def cluster(F):
    F.sort(key=lambda f: f[0]-f[1]) #sort by diagonal

    

def chain(F):

    F.sort(key=lambda f: f[0]-f[1]) #sort by diagonal

    F.sort(key=lambda f: f[0]+f[1]) #sort by anti-diagonal
    
    score_at=[0]*len(F) #initialize score at with 0

    q=[] #priority queue

    for i in xrange(len(F)):
        print F[i], F[i][0]+F[i][1], F[i][0]-F[i][1]
        if i>=1:
            assert(F[i][0]>=F[i-1][0])
            assert(F[i][1]>=F[i-1][1])

        score_at[i]=0


def gap_cost(f1,f2,lmbda=2,eps=1):
    pass



























































