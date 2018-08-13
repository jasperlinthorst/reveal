import os
import sys
from matplotlib import pyplot as plt
import subprocess
import pickle
import numpy as np
import networkx as nx
from reveal import utils
from ete3 import Tree
import uuid
import datetime
import random
import argparse
import logging

def mut(seq, seqids, idoffset, rate=0.0001, indelfrac=0.2, zipfd=1.7, maxindellength=2000,):

    p=[1/float(len(seq))]*len(seq) #uniform probability per position, TODO: generate regions that are more or less variable

    npos=int(rate*len(seq)) #number of positions to sample

    positions=np.random.choice(len(seq), npos, replace=False, p=p)
    positions=sorted(positions)

    mutations=[]
    
    mutseq=""
    nseqids=[]
    offset=0

    for pos in positions:

        if pos<offset: #position was previously removed by indel
            continue
        
        if np.random.random()<indelfrac:
            if np.random.random()<0.5:
                muttype="ins"
                length=np.random.zipf(zipfd)
                if length>maxindellength:
                    length=maxindellength
                insseq="".join(["ACGT"[np.random.randint(0,3)] for x in range(length)])
            else:
                muttype="del"
                length=np.random.zipf(zipfd)
                if length>maxindellength:
                    length=maxindellength
                insseq=None
        else:
            length=1
            muttype="snp"
            insseq="".join(["ACGT".replace(seq[pos],"")[np.random.randint(0,2)] for x in range(length)])

        if muttype=="ins":
            mutseq+=seq[offset:pos]+insseq
            nseqids+=seqids[offset:pos]+range(idoffset,idoffset+length)
            idoffset+=length
            offset=pos
        elif muttype=="del":
            mutseq+=seq[offset:pos]
            nseqids+=seqids[offset:pos]
            offset=pos+length
        else:
            mutseq+=seq[offset:pos]+insseq
            nseqids+=seqids[offset:pos]+[idoffset]
            
            idoffset+=length
            offset=pos+length

        mutations.append((muttype,pos,length))

    mutseq+=seq[offset:]
    nseqids+=seqids[offset:]

    return mutseq,nseqids,mutations,idoffset

def reveal(run,fastas,m=20,minconf=0,nproc=1,check=False):
    
    with open("run_reveal.sh",'w') as runfile:
        cmd1="reveal rem %s -m %d -o %s &> %s.rem.out"%(" ".join(fastas),m,run,run)
        cmd2="reveal unzip -u10 %s.gfa &> %s.unzip.out"%(run,run)
        cmd3="reveal refine --all --minconf=%d --nproc=%d %s.gfa &> %s.refine.out"%(minconf,nproc,run,run)
        runfile.write(cmd1+"\n")
        runfile.write(cmd2+"\n")
        runfile.write(cmd3+"\n")

    t1 = datetime.datetime.now()
    subprocess.check_output(["sh run_reveal.sh"],shell=True)
    t2 = datetime.datetime.now()
    
    gfafile=run+".realigned.gfa"
    if check:
        for f in fastas:
            subprocess.check_output(["reveal extract %s %s > tmp.fa"%(gfafile, f.replace(".fasta",""))],shell=True)
            subprocess.check_call(["diff tmp.fa %s -I \'>\'"%f],shell=True)

    return compare(gfafile,run),(t2-t1).total_seconds()

def pecan(run,fastas,minconf=0,check=False):

    with open("run_pecan.sh",'w') as runfile:
        cmd="pecan -G %s.fasta -F %s"%(run," ".join(fastas))
        runfile.write(cmd+"\n")

    t1 = datetime.datetime.now()
    subprocess.check_output(["sh run_pecan.sh"],shell=True)
    t2 = datetime.datetime.now()

    cmd="reveal convert %s.fasta --aligned"%run
    logging.debug("CMD: %s"%cmd)
    subprocess.check_output([cmd],shell=True)

    gfafile="%s.gfa"%run

    if check:
        for f in fastas:
            subprocess.check_output(["reveal extract %s %s > tmp.fa"%(gfafile, f.replace(".fasta",""))],shell=True)
            subprocess.check_call(["diff tmp.fa %s -I \'>\'"%f],shell=True)

    return compare(gfafile,run),(t2-t1).total_seconds()

def mugsy(run,fastas,check=False):

    with open("run_mugsy.sh",'w') as runfile:
        cmd="mugsy --directory %s --prefix %s %s &> %s_mugsy.out"%(os.getcwd(),run," ".join(fastas),run)
        runfile.write(cmd+"\n")
    
    for i in range(10):
        try:
            t1 = datetime.datetime.now()
            out=subprocess.check_output(["sh run_mugsy.sh"],shell=True)
            t2 = datetime.datetime.now()
        except:
            print "Mugsy run failed! (%d) retry..."%i
            continue
        break
    else:
        print "Mugsy failed after 10 retries, abort"
        sys.exit(1)
    
    cmd="reveal convert %s.maf"%run
    logging.debug(cmd)
    subprocess.check_output([cmd],shell=True)

    gfafile="%s.gfa"%run
    if check:
        for f in fastas:
            subprocess.check_output(["reveal extract %s %s > tmp.fa"%(gfafile, f.replace(".fasta",""))],shell=True)
            subprocess.check_call(["diff tmp.fa %s -I \'>\'"%f],shell=True)
    
    return compare(gfafile,run),(t2-t1).total_seconds()

def compare(gfafile,run):
    logging.debug("Checking accuracy for graph: %s"%gfafile)

    alignedg=nx.DiGraph()
    utils.read_gfa(gfafile, None, None, alignedg)
    
    genomeids=dict()
    genomes=dict()
    for pid in alignedg.graph['id2path']:
        path=alignedg.graph['id2path'][pid]
        try:
            with open(path+".seqids") as f:
                seqids=[int(v) for v in f.readline().split(',')]
            for name,seq in utils.fasta_reader(path+".fasta"):
                pass            
            genomeids[pid]=seqids
            genomes[pid]=seq
        except:
            print "Can't open file: %s.seqids"%path
            sys.exit(1)

    allpairs=dict()
    ids=sorted(alignedg.graph['id2path'].keys())
    i2p=alignedg.graph['id2path']
    for i in range(len(ids)):
        for j in range(i,len(ids)):
            p=len(set(genomeids[i]) & set(genomeids[j]))
            n=len(set(genomeids[i]) ^ set(genomeids[j]))

            key=tuple(sorted((i2p[ids[i]],i2p[ids[j]])))
            allpairs[key]={'tp':0,'fn':0,'p':p,'fp':0,'tn':0,'n':n,'_p':0,'_n':0}

    incorrectbases=0
    correctbases=0
    totbases=sum([len(genomeids[g]) for g in genomeids])
    for node,data in alignedg.nodes(data=True):
        pids=[]
        l=len(data['seq'])
        okeys=sorted(data['offsets'].keys())

        for i in range(len(okeys)):
            gi=alignedg.graph['id2path'][i]
            for j in range(i,len(okeys)):
                gj=alignedg.graph['id2path'][j]

                o=data['offsets'][okeys[i]]
                ids=genomeids[okeys[i]][o:o+l]
                seqa=genomes[okeys[i]][o:o+l]

                no=data['offsets'][okeys[j]]
                nids=genomeids[okeys[j]][no:no+l]
                seqb=genomes[okeys[j]][no:no+l]

                # tp=len(set(ids) & set(nids)) #bases correctly aligned
                # fp=len(set(ids) ^ set(nids)) #bases that are incorrectly aligned
                
                if (len(ids)!=len(nids)):
                    print seqa
                    print seqb
                    print okeys[i],okeys[j],o,no,l, len(genomeids[okeys[i]]), len(genomeids[okeys[j]])
                
                assert(len(ids)==len(nids))
                assert(len(seqa)==len(seqb))

                tp=sum([1 for x,y in zip(ids,nids) if x==y])
                fp=sum([1 for x,y in zip(ids,nids) if x!=y])

                if fp>0:
                    if fp==len(ids) or fp==len(nids): #entire node is incorrect!
                        logging.debug("Node: %s (length=%d, pair=(%s <-> %s)) is probably a spurious anchor"%(node,l,gi,gj))
                    elif ids[0]!=nids[0] or ids[-1]!=nids[-1]:
                        a=0
                        while ids[a]!=nids[a]:
                            a+=1
                        b=1
                        while ids[-b]!=nids[-b]:
                            b+=1
                        logging.debug("Node %s (length=%d, pair=(%s <-> %s)) introduces error because of edge wander (l=%d r=%d)"%(node,l,gi,gj,a,b-1))
                    else:
                        logging.debug("Node %s (length=%d, pair=(%s <-> %s)) introduces error because of convergent evolution"%(node,l,gi,gj))

                key=tuple(sorted((i2p[okeys[i]],i2p[okeys[j]])))
                allpairs[key]['tp']+=tp
                allpairs[key]['fp']+=fp
                incorrectbases+=fp
                correctbases+=tp

    #derive tn and fn
    ids=sorted(alignedg.graph['id2path'].keys())
    for i in range(len(ids)):
        for j in range(i,len(ids)):
            if i==j:
                continue
            pathi=i2p[ids[i]]
            pathj=i2p[ids[j]]
            
            key=tuple(sorted((i2p[ids[i]],i2p[ids[j]])))

            tp=allpairs[key]['tp']
            fp=allpairs[key]['fp']
            p=allpairs[key]['p']
            n=allpairs[key]['n']
            fn=p-tp
            tn=n-fp
            allpairs[key]['tn']=tn
            allpairs[key]['fn']=fn
            allpairs[key]['_p']=tp+fp
            allpairs[key]['_n']=tn+fn
            assert(allpairs[key]['_p']+allpairs[key]['_n']==allpairs[key]['p']+allpairs[key]['n'])

    return allpairs

def performance2matrices(performance,treemap,n):
    matrices=dict()
    Dbias=np.ones([n,n])/1000 #prevent divide by zero
    Dspec=np.zeros([n,n])
    Dsens=np.zeros([n,n])
    Dtp=np.zeros([n,n])
    Dfp=np.zeros([n,n])
    Dtn=np.zeros([n,n])
    Dfn=np.zeros([n,n])
    Dp=np.zeros([n,n])
    Dn=np.zeros([n,n])
    D_p=np.zeros([n,n])
    D_n=np.zeros([n,n])
    assert(len(performance)==((n*(n-1))/2)+n)
    for id1,id2 in performance:
        table=performance[(id1,id2)]
        key=tuple(sorted([treemap[id1],treemap[id2]]))
        Dtp[key]=table['tp']
        Dfp[key]=table['fp']
        Dtn[key]=table['tn']
        Dfn[key]=table['fn']
        Dp[key]=table['p']
        D_p[key]=table['_p']
        Dn[key]=table['n']
        D_n[key]=table['_n']
        Dtp[key[::-1]]=table['tp']
        Dfp[key[::-1]]=table['fp']
        Dtn[key[::-1]]=table['tn']
        Dfn[key[::-1]]=table['fn']
        Dp[key[::-1]]=table['p']
        D_p[key[::-1]]=table['_p']
        Dn[key[::-1]]=table['n']
        D_n[key[::-1]]=table['_n']

    Dspec=(Dtn+Dbias+Dbias)/((Dtn+Dbias)+(Dfp+Dbias))
    Dsens=(Dtp+Dbias+Dbias)/((Dtp+Dbias)+(Dfn+Dbias))

    matrices['tp']=Dtp
    matrices['fp']=Dfp
    matrices['tn']=Dtn
    matrices['fn']=Dfn
    matrices['p']=Dp
    matrices['n']=Dn
    matrices['_p']=Dtp
    matrices['_n']=Dtp
    matrices['spec']=Dspec
    matrices['sens']=Dsens

    return matrices
    # return Dtp,Dfp,Dtn,Dfn,Dp,Dn,D_p,D_n,Dspec,Dsens

def matrices2summary(matrices):
    #Dtp,Dfp,Dtn,Dfn,Dp,Dn,D_p,D_n,Dspec,Dsens
    summary=dict()

    if sum(sum(matrices['n']))==0:
        specificity=1
    else:
        specificity=sum(sum(matrices['tn']))/float(sum(sum(matrices['n'])))

    if sum(sum(matrices['p']))==0:
        sensitivity=1
    else:
        sensitivity=sum(sum(matrices['tp']))/float(sum(sum(matrices['p'])))

    if sum(sum(matrices['tp']))+sum(sum(matrices['fp']))==0:
        precision=1
    else:
        precision=sum(sum(matrices['tp']))/float(sum(sum(matrices['tp']))+sum(sum(matrices['fp'])))

    f1=(2*((specificity*sensitivity)/(specificity+sensitivity)))

    summary['specificity']=specificity
    summary['precision']=precision
    summary['sensitivity']=sensitivity
    summary['f1']=f1

    return summary

def simulate(seq,tree,args):
    ancestralgenomes={"root":(seq,range(len(seq)))}

    idoffset=len(seq)
    genomes=dict()

    totd=0
    for node in tree.traverse("preorder"):
        if node.up==None:
            continue
        totd+=float(node.dist)

    totmut=[]
    for node in tree.traverse("preorder"):
        if node.up==None:
            node.name="root"
            continue

        d=(node.dist/float(totd))*args.mrate

        if node.is_leaf():
            name=str(node).replace("\n","").replace("--","")
            logging.debug("Simulating genome: %s %s"%(name,d))
            mseq,mseqids,mutations,idoffset=mut(ancestralgenomes[node.up.name][0], ancestralgenomes[node.up.name][1], idoffset, rate=d, indelfrac=args.indelfrac)
            genomes[name]=(mseq,mseqids)
        else:
            name=str(uuid.uuid4().hex)
            node.name=name
            logging.debug("Simulating ancestral genome: %s %s"%(name,d))
            mseq,mseqids,mutations,idoffset=mut(ancestralgenomes[node.up.name][0], ancestralgenomes[node.up.name][1], idoffset, rate=d, indelfrac=args.indelfrac)
            ancestralgenomes[name]=(mseq,mseqids)

        totmut+=mutations

    logging.info("Total number of mutations in the simulated population: %d"%len(totmut))

    fastas=[]
    for genome in genomes:
        fn=genome+".fasta"
        fastas.append(fn)
        utils.fasta_writer(fn,[(genome,genomes[genome][0])])
        with open(genome+".seqids",'w') as seqidsf:
            for sid in genomes[genome][1][:-1]:
                seqidsf.write("%d,"%sid)
            seqidsf.write("%d"%genomes[genome][1][-1])

    return fastas

def main():
    desc="""
    Generates a synthetic dataset with a uniform distribution of variants in a specified population size,
    for which indels make up a fixed fraction of the total number of variants in the population.

    Then aligns the resulting sequence with REVEAL, MUGSY and Pecan and reports on the performance in terms
    of fp, tp, fn and tn as well as the runtime.

    Performance parameters are stored in a single pickle file.
    """
    
    parser = argparse.ArgumentParser(prog="simulate", usage="simulate.py -h for usage", description=desc)
    parser.add_argument("-l", "--log-level", type=int, dest="loglevel", default=20, help="Log level: 1=trace 10=debug 20=info 30=warn 40=error 50=fatal.")
    parser.add_argument("ancestral", type=str, help="Ancestral genome that corresponds to the root in the phylogenetic tree in fasta format.")
    parser.add_argument("-n", dest="n", type=int, required=True, help="Generate a population of this many samples (default 4).")
    parser.add_argument("-r", dest="mrate", type=int, required=True, help="Rate of variation to be simulated, a value of 100 means on average 1 variant every 100 bases.")
    parser.add_argument("-i", dest="indelfrac", type=int, required=True, help="Percentage of variants that is an indel.")
    parser.add_argument("-s", dest="seed", type=int, required=True, help="Seed value for generating random variants.")
    parser.add_argument("-c", dest="minconf", type=int, default=0, help="Confidence cut-off for REVEAL.")
    parser.add_argument("--nproc", dest="nproc", type=int, default=1, help="Use multi-processing for REVEAL.")
    parser.add_argument("--force", dest="force", action="store_true", default=False, help="Force new alignments even when the same experiment was already performed.")
    parser.add_argument("--clean", dest="clean", action="store_true", default=False, help="Remove all files expect the summary pickle file after the experiment has run.")
    parser.add_argument("--tmp", dest="tmpdir", default=None, help="Use tmp-dir to write fasta and seqid files.")
    parser.add_argument("--check", dest="check", action="store_true", default=False, help="Check that gfa files properly encode the input sequence.")
    parser.add_argument("-t", dest="tree", type=str, default=None, help="Use a predefined phylogenetic tree (in newick format).")

    args = parser.parse_args()
    logging.basicConfig(format='%(asctime)s %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p', level=args.loglevel)

    runid="%s_n%s_r%s_i%s_s%s_c%s"%(os.path.basename(args.ancestral).replace(".","_"), 
                            args.n, args.mrate, args.indelfrac, args.seed, args.minconf)

    logging.info("Starting run with id: %s"%runid)

    for name,seq in utils.fasta_reader(args.ancestral):
        break #just read one sequence

    #re-encode rates and fractions
    args.mrate= (len(seq)/float(args.mrate)) / len(seq)
    
    args.indelfrac=args.indelfrac/float(100)

    if os.path.exists(runid+"/"+runid+".pickle") and not args.force:
        logging.warn("Experiment already was already performed, quit")
        sys.exit(0)
    
    if not os.path.exists(runid):
        os.mkdir(runid)

    if args.tmpdir!=None:
        pwd=os.getcwd()
        os.chdir(args.tmpdir)
        if not os.path.exists(runid):
            os.mkdir(runid)

    os.chdir(runid)

    if args.seed!=None:
        random.seed(args.seed)
        np.random.seed(seed=args.seed)
    
    if args.tree==None:
        t=Tree()
        t.populate(args.n,[runid+"_"+str(v) for v in range(args.n)],random_branches=True)#,branch_range=(0,args.mrange))
    else:
        t=Tree(args.tree)
    
    treemap=dict() #map the leaves/genomes in the tree to ids
    for i,node in enumerate(t.get_leaves()):
        treemap[node.name]=i

    # ts = TreeStyle()
    # ts.show_leaf_name = True
    # ts.show_branch_length = False
    # ts.show_branch_support = False
    # t.convert_to_ultrametric()
    #render the (simulated) phylogenetic tree as an image
    # t.render("simtree_%d_%.2f.png"%(args.n,args.mrange),tree_style=ts)

    #write the tree to nwk format in case we want to check afterwards
    t.write(outfile=runid+".nwk")
    
    #simulate fasta and seqid files
    genomes=simulate(seq,t,args)

    #the syntethic genomes sequences are generated, now produce some graphs
    performance=dict()

    #what were the commandline arguments for this run
    performance['args']=args

    #RUN REVEAL
    allpairs,runtime=reveal(runid+"_reveal",genomes,minconf=args.minconf,nproc=args.nproc,check=args.check)
    performance['reveal_matrices']=performance2matrices(allpairs,treemap,args.n)
    performance['reveal_runtime']=runtime
    performance['reveal_summary']=matrices2summary(performance['reveal_matrices'])

    print "reveal",\
            "n=%d"%args.n,\
            "r=%d"%args.mrate,\
            "s=%d"%args.seed,\
            "i=%d"%args.indelfrac,\
            "c=%d"%args.minconf,\
            "tp=%d"%sum(sum(performance['reveal_matrices']['tp'])), \
            "fp=%d"%sum(sum(performance['reveal_matrices']['fp'])), \
            "tn=%d"%sum(sum(performance['reveal_matrices']['tn'])), \
            "fn=%d"%sum(sum(performance['reveal_matrices']['fn'])), \
            "f1=%.5f"%performance['reveal_summary']['f1'], \
            "sensitivity=%.5f"%performance['reveal_summary']['sensitivity'], \
            "specificity=%.5f"%performance['reveal_summary']['specificity'], \
            "precision=%.5f"%performance['reveal_summary']['precision'], \
            "runtime=%.2f"%runtime \


    #RUN MUGSY
    allpairs,runtime=mugsy(runid+"_mugsy",genomes,check=args.check)
    performance['mugsy_matrices']=performance2matrices(allpairs,treemap,args.n)
    performance['mugsy_runtime']=runtime
    performance['mugsy_summary']=matrices2summary(performance['mugsy_matrices'])
    print "mugsy",\
            "n=%d"%args.n,\
            "r=%d"%args.mrate,\
            "s=%d"%args.seed,\
            "i=%d"%args.indelfrac,\
            "c=%d"%args.minconf,\
            "tp=%d"%sum(sum(performance['mugsy_matrices']['tp'])), \
            "fp=%d"%sum(sum(performance['mugsy_matrices']['fp'])), \
            "tn=%d"%sum(sum(performance['mugsy_matrices']['tn'])), \
            "fn=%d"%sum(sum(performance['mugsy_matrices']['fn'])), \
            "f1=%.5f"%performance['mugsy_summary']['f1'], \
            "sensitivity=%.5f"%performance['mugsy_summary']['sensitivity'], \
            "specificity=%.5f"%performance['mugsy_summary']['specificity'], \
            "precision=%.5f"%performance['mugsy_summary']['precision'], \
            "runtime=%.2f"%runtime \


    #RUN PECAN
    allpairs,runtime=pecan(runid+"_pecan",genomes,check=args.check)#,minconf=conf)
    performance['pecan_matrices']=performance2matrices(allpairs,treemap,args.n)
    performance['pecan_runtime']=runtime
    performance['pecan_summary']=matrices2summary(performance['pecan_matrices'])
    print "pecan",\
            "n=%d"%args.n,\
            "r=%d"%args.mrate,\
            "s=%d"%args.seed,\
            "i=%d"%args.indelfrac,\
            "c=%d"%args.minconf,\
            "tp=%d"%sum(sum(performance['pecan_matrices']['tp'])),  \
            "fp=%d"%sum(sum(performance['pecan_matrices']['fp'])),  \
            "tn=%d"%sum(sum(performance['pecan_matrices']['tn'])),  \
            "fn=%d"%sum(sum(performance['pecan_matrices']['fn'])), \
            "f1=%.5f"%performance['pecan_summary']['f1'], \
            "sensitivity=%.5f"%performance['pecan_summary']['sensitivity'], \
            "specificity=%.5f"%performance['pecan_summary']['specificity'], \
            "precision=%.5f"%performance['pecan_summary']['precision'], \
            "runtime=%.2f"%runtime \

    if args.tmpdir!=None:
        for file in os.listdir(os.getcwd()):
            os.remove(file) #clean up
        os.chdir(pwd+"/"+runid)

    pickle.dump(performance,open(runid+".pickle", "wb"))

    if args.clean:
        for file in os.listdir(os.getcwd()):
            if file!=runid+".pickle":
                os.remove(file)

if __name__ == '__main__':
    main()








#         fig = plt.figure()
#         ax = plt.subplot(111)
#         plt.title("Sensitivity (conf=%s)"%conf)
#         for i,rate in enumerate(rates):
#             plt.plot([rate*.999]*niterations,psens[i],'rx',label="PECAN sens")
#             plt.plot([rate]*niterations,rmsens[i],'gx',label="REVEAL sens")
#             plt.plot([rate*1.001]*niterations,msens[i],'bx',label="MUGSY sens")
#         plt.xticks(rates)

#         box = ax.get_position()
#         ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
#         ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
#         plt.savefig("sens_n=%d_m=%d_p=%.10f_i=%.5f_conf=%d.png"%(n,m,p,indelrate,conf))

        
#         # plt.clf()
#         # plt.title("Sensitivity (conf=%s)"%conf)
#         # for i,rate in enumerate(rates):
#         #     plt.plot([rate]*niterations,rmsens[i],'gx',label="REVEAL sens")
#         #     plt.plot([rate*.999]*niterations,psens[i],'rx',label="PECAN sens")
#         #     plt.plot([rate*1.001]*niterations,msens[i],'bx',label="MUGSY sens")    
#         # fig,ax=plt.subplots()
#         # ax.set_xticks(rates)
#         # ax.set_xticklabels([str(r) for r in rate])
#         # plt.savefig("sens_n=%d_m=%d_p=%.10f_i=%.5f_conf=%d.png"%(n,m,p,indelrate,conf))








#         fig = plt.figure()
#         ax = plt.subplot(111)
#         plt.title("Specificity (conf=%s)"%conf)
#         for i,rate in enumerate(rates):
#             plt.plot([rate*.999]*niterations,pspec[i],'rx',label="PECAN spec")
#             plt.plot([rate]*niterations,rmspec[i],'gx',label="REVEAL spec")
#             plt.plot([rate*1.001]*niterations,mspec[i],'bx',label="MUGSY spec")
#         plt.xticks(rates)
#         box = ax.get_position()
#         ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
#         ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
#         plt.savefig("spec_n=%d_m=%d_p=%.10f_i=%.5f_conf=%d.png"%(n,m,p,indelrate,conf))




#         # plt.clf()
#         # plt.title("Specificity (conf=%s)"%conf)
#         # for i,rate in enumerate(rates):
#         #     plt.plot([rate]*niterations,rmspec[i],'gx',label="REVEAL spec")
#         #     plt.plot([rate*.999]*niterations,pspec[i],'rx',label="PECAN spec")
#         #     plt.plot([rate*1.001]*niterations,mspec[i],'bx',label="MUGSY spec")
#         # fig,ax=plt.subplots()
#         # ax.set_xticks(rates)
#         # ax.set_xticklabels([str(r) for r in rate])
#         # plt.savefig("spec_n=%d_m=%d_p=%.10f_i=%.5f_conf=%d.png"%(n,m,p,indelrate,conf))




#         fig = plt.figure()
#         ax = plt.subplot(111)
#         plt.title("F1 (conf=%s)"%conf)
#         for i,rate in enumerate(rates):
#             plt.plot([rate*.999]*niterations,pf1[i],'rx',label="PECAN f1")
#             plt.plot([rate]*niterations,rmf1[i],'gx',label="REVEAL f1")
#             plt.plot([rate*1.001]*niterations,mf1[i],'bx',label="MUGSY f1")
#         plt.xticks(rates)
#         box = ax.get_position()
#         ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
#         ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
#         plt.savefig("f1_n=%d_m=%d_p=%.10f_i=%.5f_conf=%d.png"%(n,m,p,indelrate,conf))


        #fig = plt.figure()
        #print [pf1,rmf1,mf1]
        #plt.boxplot([pf1[0],rmf1[0],mf1[0]],labels=["pecan","reveal","mugsy"])
        #plt.savefig("f1_boxplot_n=%d_m=%d_p=%.10f_i=%.5f_conf=%d.png"%(n,m,p,indelrate,conf))




        # plt.clf()
        # plt.title("F1 (conf=%s)"%conf)
        # for i,rate in enumerate(rates):
        #     plt.plot([rate]*niterations,rmf1[i],'gx',label="REVEAL f1")
        #     plt.plot([rate*.999]*niterations,pf1[i],'rx',label="PECAN f1")
        #     plt.plot([rate*1.001]*niterations,mf1[i],'bx',label="MUGSY f1")
        # fig,ax=plt.subplots()
        # ax.set_xticks(rates)
        # ax.set_xticklabels([str(r) for r in rate])
        # plt.savefig("f1_n=%d_m=%d_p=%.10f_i=%.5f_conf=%d.png"%(n,m,p,indelrate,conf))
        
        # plt.plot(rates,rsens,'b-.',label="REVEAL sens")
        # plt.plot(rates,rspec,'r-.',label="REVEAL spec")
        # plt.plot(rates,rf1,'g-.',label="REVEAL f1")

        # plt.plot(rates,rmsens,'g--',label="REVEAL sens")
        # plt.plot(rates,rmspec,'g:',label="REVEAL spec")
        # plt.plot(rates,rmf1,'g-',label="REVEAL f1")

        # plt.plot(rates,psens,'r--',label="PECAN sens")
        # plt.plot(rates,pspec,'r:',label="PECAN spec")
        # plt.plot(rates,pf1,'r-',label="PECAN f1")

        # plt.plot(rates,msens,'b--',label="MUGSY sens")
        # plt.plot(rates,mspec,'b:',label="MUGSY spec")
        # plt.plot(rates,mf1,'b-',label="MUGSY f1")

        # plt.legend(loc="lower left")
        # plt.show()
        # plt.savefig("perf_n=%d_m=%d_p=%.10f_i=%.5f_conf=%d.png"%(n,m,p,indelrate,conf))
        










        

        # fig = plt.figure(figsize=(12,8))
        # ax = plt.subplot(111)

        # plt.ylim(0,max([max(rtp),max(rtn),max(rfp),max(rfn),max(rmtp),max(rmtn),max(rmfp),max(rmfn),max(mtp),max(mtn),max(mfp),max(mfn),max(ptp),max(ptn),max(pfp),max(pfn)]))
        # plt.title("Alignment performance (for anchor size %sbp)"%m)

        # # ax.plot(rates,rtp,'r-',label="REVEAL tp")
        # # ax.plot(rates,rtn,'g-',label="REVEAL tn")
        # # ax.plot(rates,rfp,'b-',label="REVEAL fp")
        # # ax.plot(rates,rfn,'c-',label="REVEAL fn")

        # ax.plot(rates,rmtp,'r-',label="REVEAL + PROBCONS tp")
        # ax.plot(rates,rmtn,'g-',label="REVEAL + PROBCONS tn")
        # ax.plot(rates,rmfp,'b-',label="REVEAL + PROBCONS fp")
        # ax.plot(rates,rmfn,'c-',label="REVEAL + PROBCONS fn")

        # ax.plot(rates,ptp,'r:',label="PECAN tp")
        # ax.plot(rates,ptn,'g:',label="PECAN tn")
        # ax.plot(rates,pfp,'b:',label="PECAN fp")
        # ax.plot(rates,pfn,'c:',label="PECAN fn")

        # ax.plot(rates,mtp,'r--',label="MUGSY tp")
        # ax.plot(rates,mtn,'g--',label="MUGSY tn")
        # ax.plot(rates,mfp,'b--',label="MUGSY fp")
        # ax.plot(rates,mfn,'c--',label="MUGSY fn")

        # # Shrink current axis by 20%
        # box = ax.get_position()
        # ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
        # ax.legend(loc="center left", bbox_to_anchor=(1, 0.5))    
        # # plt.show()

        # plt.savefig("tp_tn_fp_fn_n=%d_m=%d_p=%.10f_i=%.5f_conf=%d.png"%(n,m,p,indelrate,conf))
















        # plt.clf()
        # plt.ylim(-0.1,max([max(rmfp),max(rfp),max(mfp),max(pfp)]))
        
        # plt.title("False positives (for anchor size %sbp) REVEAL"%(m))
        
        # print rfp

        # # plt.plot(rates,rfp,'b--',label="REVEAL (REM) fp")
        # plt.plot(rates,rmfp,'b-',label="REVEAL fp")
        
        # # plt.plot(rates,pfp,'g-',label="PECAN fp")
        # # plt.plot(rates,mfp,'r-',label="MUGSY fp")


        # plt.legend(loc="upper left")
        # plt.savefig("fp_mugsy_reveal_n=%d_m=%d_p=%.10f_i=%.5f_conf=%d.png"%(n,m,p,indelrate,conf))






    # print confrange, rmfp, rmfn

    # plt.clf()

    # plt.ylim(-0.1, max([max(rmfp)+1,max(rmfn)+1]))

    # plt.title("Confidence vs fp/fn REVEAL")

    # plt.plot(confrange,rmfp,'r-',label="fp")
    # plt.plot(confrange,rmfn,'b-',label="fn")

    # plt.legend(loc="upper left")
    # plt.savefig("confidence_reveal_fp_fn_n=%d_m=%d_p=%.10f_i=%.5f_conf=%d.png"%(n,m,p,indelrate,conf))





    # plt.clf()
    # plt.ylim(0.5,1.01)

    # plt.title("Confidence vs sens/spec/f1 REVEAL")

    # plt.plot(confrange,rmsens,'r-',label="sens")
    # plt.plot(confrange,rmspec,'g-',label="spec")
    # plt.plot(confrange,rmf1,'b-',label="f1")

    # plt.legend(loc="upper left")
    # plt.savefig("confidence_reveal_sens_spec_f1_n=%d_m=%d_p=%.10f_i=%.5f_conf=%d.png"%(n,m,p,indelrate,conf))














    # plt.show()
