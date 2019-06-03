from utils import *
from extract import extract,extract_path
from rem import align,prune_nodes
from random import shuffle
import bubbles
import schemes
import multiprocessing as mp
import signal
import probconslib
import math
import os
import time

def refine_bubble_cmd(args):
    if len(args.graph)<1:
        logging.fatal("Specify a gfa file for which bubbles should be realigned.")
        return
    
    # G=nx.MultiDiGraph() #TODO: make sure that refine can handle structural variant edges, so make sure we use a MultiDiGraph here!
    G=nx.DiGraph()

    logging.info("Reading graph...")
    read_gfa(args.graph[0],None,"",G)
    logging.info("Done.")

    logging.info("Paths through the graph: %s"%G.graph['paths'])

    if (args.source==None and args.sink==None) and (args.all or args.complex or args.simple):
        G=refine_all(G, **vars(args))
    else:
        if args.source==None or args.sink==None:
            logging.error("Specify source sink pair, or one of the --all --simple --complex flags.")
            sys.exit(1)

        args.source,args.sink=int(args.source),int(args.sink)
        source_idx,sink_idx=None,None

        topiter=nx.topological_sort(G)
        for i,v in enumerate(topiter):
            if v==args.source:
                source_idx=i
                break
        nodes=[args.source]
        for j,v in enumerate(topiter):
            nodes.append(v)
            if v==args.sink:
                sink_idx=i+j+1
                break
        
        if source_idx==None or sink_idx==None:
            logging.fatal("Unkown source/sink pair: %d,%d"%(args.source,args.sink))

        b=bubbles.Bubble(G,args.source,args.sink,source_idx,sink_idx,nodes)

        nn=max([node for node in G.nodes() if type(node)==int])+1

        bnodes=list(set(b.nodes)-set([b.source,b.sink]))
        sg=G.subgraph(bnodes)

        offsets=dict()
        for sid in G.node[b.source]['offsets']:
            offsets[sid]=G.node[b.source]['offsets'][sid]+len(G.node[b.source]['seq'])

        sourcesamples=set(G.node[b.source]['offsets'].keys())
        sinksamples=set(G.node[b.sink]['offsets'].keys())
        paths=sourcesamples.intersection(sinksamples)

        G.node[b.source]['aligned']=1
        G.node[b.sink]['aligned']=1

        res=refine_bubble(sg,b,offsets,paths, **vars(args))

        if res!=None:
            bubble,ng,path2start,path2end=res
            G,nn=replace_bubble(G,b,ng,path2start,path2end,nn)

    if args.outfile==None:
        fn=args.graph[0].replace(".gfa",".realigned.gfa")
    else:
        fn=args.outfile

    logging.info("Prune and contract nodes...")
    prune_nodes(G)
    contract(G,[n for n in nx.topological_sort(G) if type(n)!=str])
    logging.info("Done.")

    logging.info("Write refined graph to: %s"%fn)
    write_gfa(G,"",outputfile=fn)

def replace_bubble(G,bubble,ng,path2start,path2end,nn):
    assert(nn not in G)

    bubblenodes=bubble.nodes[1:-1] #exclude source and sink node
    
    for node in bubblenodes: #remove all bubblenodes from the original graph
        G.remove_node(node)

    mapping={}
    for node in ng.nodes(): #add nodes from newly aligned graph to original graph
        mapping[node]=nn
        nn+=1
    
    ng=nx.relabel_nodes(ng,mapping) #relabel nodes according to a unique integer range

    for node,data in ng.nodes(data=True): #add nodes from newly aligned graph to original graph
        G.add_node(node,**data)

    for edge in ng.edges(data=True):
        G.add_edge(edge[0],edge[1],**edge[2])

    for sid in path2start:
        startnode=mapping[path2start[sid][0]]
        if G.has_edge(bubble.source,startnode):
            G[bubble.source][startnode]['paths'].add(sid)
        else:
            G.add_edge(bubble.source,startnode,ofrom='+',oto='+',paths=set([sid]))

    for sid in path2end:
        endnode=mapping[path2end[sid][0]]
        if G.has_edge(endnode,bubble.sink):
            G[endnode][bubble.sink]['paths'].add(sid)
        else:
            G.add_edge(endnode,bubble.sink,ofrom='+',oto='+',paths=set([sid]))


    #SHOULD BE ABLE TO SKIP THIS HERE, AS THIS WILL BE DONE OVER THE ENTIRE GRAPH ANYWAY...
    #Just one possible path from source to start, contract nodes
    # if len(G.out_edges(bubble.source))==1 and type(bubble.source)!=str:
    #     # assert(len(set(path2start.values()))==1)
    #     startnode=mapping[path2start.values()[0][0]]
    #     G.node[bubble.source]['seq']+=G.node[startnode]['seq']
    #     for to in G[startnode]:
    #         d=G[startnode][to]
    #         G.add_edge(bubble.source,to,**d)
    #     G.remove_node(startnode)

    #Just one possible path from end to sink, contract nodes
    # if len(G.in_edges(bubble.sink))==1 and type(bubble.sink)!=str:
    #     # assert(len(set(path2end.values()))==1)
    #     endnode=mapping[path2end.values()[0][0]]
    #     G.node[bubble.sink]['seq']=G.node[endnode]['seq']+G.node[bubble.sink]['seq']
    #     G.node[bubble.sink]['offsets']=G.node[endnode]['offsets']
    #     for e0,e1,d in G.in_edges(endnode,data=True):
    #         G.add_edge(e0,bubble.sink,**d)
    #     G.remove_node(endnode)

    return G,nn

def refine_bubble(sg,bubble,offsets,paths,**kwargs):
    source=bubble.source
    sink=bubble.sink

    if len(bubble.nodes)==3:
        logging.fatal("Indel bubble, no point realigning.")
        return

    #TODO: if bubble contains structural variant edge, track these or simply refuse realignment!

    aobjs=[]

    uniqueonly=False

    t0=time.time()
    if kwargs['uniqueonly']:
        d={}
        #extract all paths
        for sid in paths:
            seq=extract(sg,sg.graph['id2path'][sid])
            if len(seq)>0:
                if seq in d:
                    d[seq].append(str(sid))
                else:
                    d[seq]=[str(sid)]

        if len(d)==1:
            logging.debug("Nothing to refine for bubble: %s - %s"%(bubble.source,bubble.sink))
            return

        aobjs=[(",".join(d[seq]),seq) for seq in d]

    else:
        for sid in paths:
            seq=extract(sg,sg.graph['id2path'][sid])
            if len(seq)>0:
                aobjs.append((str(sid),seq))

                # if seq in d:
                #     d[seq].append(str(sid))
                # else:
                #     d[seq]=[str(sid)]

        if len(aobjs)==1:
            logging.debug("Nothing to refine for bubble: %s - %s"%(bubble.source,bubble.sink))
            return
    t1=time.time()
    logging.debug("Extracting sequence for paths: %s through bubble <%s,%s> took: %.4f seconds."%(paths,bubble.source,bubble.sink,t1-t0))


    # logging.info("Realigning bubble (pid=%s) between <%s> and <%s>, %d alleles, with %s (max size %dbp, in nodes=%d)."%(os.getpid(),bubble.source,bubble.sink,len(aobjs),kwargs['method'],bubble.maxsize,len(bubble.nodes)-2))

    # for name,seq in aobjs:
    #     if len(seq)>200:
    #         logging.debug("IN %s: %s...%s (%d bp)"%(name.rjust(4,' '),seq[:100],seq[-100:],len(seq)))
    #     else:
    #         logging.debug("IN %s: %s (%d bp)"%(name.rjust(4,' '),seq,len(seq)))

    if kwargs['method']!="reveal_rem": #use custom multiple sequence aligner to refine bubble structure
        ng=msa2graph(aobjs,
                        msa=kwargs['method'],
                        minconf=kwargs['minconf'],
                        parameters=kwargs['parameters'],
                        constrans=kwargs['constrans'],
                        nrefinements=kwargs['nrefinements'],
                        consgap=kwargs['consgap']
                        )
        if ng==None:
            logging.fatal("MSA using %s for bubble: %s - %s failed."%(kwargs['method'],source,sink))
            return

    else: #use reveal with different settings
        ng,idx=align(aobjs, minlength=kwargs['minlength'],
                            minn=kwargs['minn'],
                            seedsize=kwargs['seedsize'],
                            maxmums=kwargs['maxmums'],
                            wpen=kwargs['wpen'],
                            wscore=kwargs['wscore'],
                            gcmodel=kwargs['gcmodel'],
                            sa64=kwargs['sa64'])
        T=idx.T
        seq2node(ng,T) #transfer sequence to node attributes

    #map edge atts back to original graph
    for n1,n2,data in ng.edges(data=True):
        newpaths=set()
        for p in data['paths']:
            for x in ng.graph['id2path'][p].split(','):
                newpaths.add(int(x))
        data['paths']=newpaths

    #map node atts back to original graph
    for node,data in ng.nodes(data=True):
        newoffsets={}
        for sid in data['offsets']:
            for x in ng.graph['id2path'][sid].split(','):
                newoffsets[int(x)]=data['offsets'][sid]
        data['offsets']=newoffsets

    ng.graph['paths']=sg.graph['paths']
    ng.graph['path2id']=sg.graph['path2id']
    ng.graph['id2path']=sg.graph['id2path']

    mapping={}
    
    path2start=dict()
    path2end=dict()

    #map nodes back to original offsets and idspace, and determine first/last node for every path
    for node,data in ng.nodes(data=True):
        for sid in data['offsets']:
            if sid not in path2start or data['offsets'][sid]<path2start[sid][1]:
                path2start[sid]=(node,data['offsets'][sid])

        for sid in data['offsets']:
            if sid not in path2end or data['offsets'][sid]>path2end[sid][1]:
                path2end[sid]=(node,data['offsets'][sid])

        corrected=dict()
        for sid in data['offsets']:
            corrected[sid]=data['offsets'][sid]+offsets[sid]

        ng.node[node]['offsets']=corrected

    return bubble,ng,path2start,path2end

def align_worker(G,chunk,outputq,kwargs):
    try:
        logging.info("Worker with pid=%d started on subgraph of length: %d"%(os.getpid(),len(G)))

        rchunk=[]
        for b in chunk:

            logging.debug("Start realign bubble (pid=%d) between <%s> and <%s>, max allele size %dbp (in nodes=%d)."%(os.getpid(),b.source,b.sink,b.maxsize,len(b.nodes)-2))

            G.node[b.source]['aligned']=1
            G.node[b.sink]['aligned']=1
            
            bnodes=list(set(b.nodes)-set([b.source,b.sink]))
            sg=G.subgraph(bnodes).copy()
            
            offsets=dict()
            for sid in G.node[b.source]['offsets']:
                offsets[sid]=G.node[b.source]['offsets'][sid]+len(G.node[b.source]['seq'])

            sourcesamples=set(G.node[b.source]['offsets'].keys())
            sinksamples=set(G.node[b.sink]['offsets'].keys())
            paths=sourcesamples.intersection(sinksamples)
            
            t0=time.time()
            rb=refine_bubble(sg,b,offsets,paths,**kwargs)
            t1=time.time()
            logging.debug("Realign bubble (pid=%d) between <%s> and <%s>, max allele size %dbp (in nodes=%d), took %.2f seconds."%(os.getpid(),b.source,b.sink,b.maxsize,len(b.nodes)-2,t1-t0))

            if rb==None:
                continue
            else:
                rchunk.append(rb)
                if len(rchunk)==kwargs['chunksize']:
                    t0=time.time()
                    outputq.put(rchunk)
                    t1=time.time()
                    logging.debug("Added chunk to queue in %.2f seconds."%(t1-t0))
                    rchunk=[]

        if len(rchunk)>0:
            t0=time.time()
            outputq.put(rchunk)
            t1=time.time()
            logging.debug("Added chunk to queue in %.2f seconds."%(t1-t0))

        logging.info("Worker with pid=%d is done."%(os.getpid()))
        outputq.put(-1)
    except Exception,e:
        logging.fatal("Worker with pid=%d failed at bubble <%s,%s>: %s"%(os.getpid(),b.source,b.sink,str(e)))
        exit(1)

def graph_worker(G,nn,outputq,aworkers,totbubbles):
    deadworkers=0
    nworkers=len(aworkers)
    refinedbubbles=0

    while True:

        if deadworkers==nworkers:
            break

        # try:
        t0=time.time()
        # data=outputq.get(timeout=.5)
        data=outputq.get()
        t1=time.time()
        if data!=-1:
            logging.info("Getting chunk of size %d from queue took %.4f seconds."%(len(data),t1-t0))

        # except mp.queues.Empty: #nothing to do, check if all workers are still alive...
        #     logging.info("Empty queue check on workers.")
        #     for wi in range(len(aworkers)):
        #         # p,fn,args=aworkers[wi]
        #         p=aworkers[wi]
        #         if not p.is_alive() and p.exitcode!=0: #one of the workers was killed! maybe oom...
        #             # if retry: #start a new worker and continue processing whatever is left on the queue
        #             #     logging.error("Worker %d died with exitcode: %d, start a new worker!"%(p.pid,p.exitcode))
        #             #     np=mp.Process(target=fn, args=args)
        #             #     np.start()
        #             #     aworkers[wi]=((np,fn,args)) #update it
        #             # else:
        #             raise Exception("Worker %d died with exitcode: %d. Stop refining."%(p.pid,p.exitcode))
        #     continue

        if data==-1: #worker was done
            logging.info("Worker signaled that its done.")
            deadworkers+=1
        else:
            for d in data:
                refinedbubbles+=1
                bubble,ng,path2start,path2end=d
                t0=time.time()
                G,nn=replace_bubble(G,bubble,ng,path2start,path2end,nn)
                t1=time.time()
                logging.info("Replacing bubble (%d/%d): <%s,%s> took %.4f seconds."%(refinedbubbles,totbubbles,bubble.source,bubble.sink,t1-t0))

        for wi in range(len(aworkers)):
            p=aworkers[wi]
            if not p.is_alive() and p.exitcode!=0: #one of the workers was killed! maybe oom...
                raise Exception("Worker %d died with exitcode: %d. Stop refining."%(p.pid,p.exitcode))


def refine_all(G, **kwargs):
    realignbubbles=[]
    
    if kwargs['minsize']==None:
        kwargs['minsize']=kwargs['minlength']

    #detect all bubbles
    logging.info("Extracting bubbles...")

    for b in bubbles.bubbles(G):

        if kwargs['complex']:
            if b.issimple():
                logging.debug("Skipping bubble <%s,%s>, not complex."%(b.source,b.sink))
                continue

        if kwargs['nogaps']:
            spansgap=False
            for n in b.nodes:
                if 'N' in G.nodes[n]['seq']:
                    logging.info("Skipping bubble <%s,%s>, bubble spans a gap."%(b.source,b.sink))
                    spansgap=True
                    break
            if spansgap:
                continue

        if kwargs['simple']:
            if not b.issimple():
                logging.debug("Skipping bubble <%s,%s>, not simple."%(b.source,b.sink))
                continue

        if b.maxsize-b.minsize<kwargs['mindiff']:
            logging.debug("Skipping bubble <%s,%s>, diff between smallest and largest allele (%dbp) is smaller than mindiff=%d."%(b.source,b.sink,b.maxsize-b.minsize,kwargs['mindiff']))
            continue

        if kwargs['maxdiff'] and b.maxsize-b.minsize>kwargs['maxdiff']:
            logging.debug("Skipping bubble <%s,%s>, diff between smallest and largest allele (%dbp) is larger than maxdiff=%d."%(b.source,b.sink,b.maxsize-b.minsize,kwargs['maxdiff']))
            continue

        if b.minsize<kwargs['minsize']:
            logging.debug("Skipping bubble <%s,%s>, smallest allele (%dbp) is smaller than minsize=%d."%(b.source,b.sink,b.minsize,kwargs['minsize']))
            continue

        if b.maxsize>kwargs['maxsize']:
            logging.warn("Skipping bubble <%s,%s>, largest allele (%dbp) is larger than maxsize=%d."%(b.source,b.sink,b.maxsize,kwargs['maxsize']))
            continue

        if kwargs['maxcumsize']!=None:
            if b.cumsize>kwargs['maxcumsize']:
                logging.warn("Skipping bubble <%s,%s>, cumulative size %d is larger than maxcumsize=%d."%(b.source,b.sink,b.cumsize,kwargs['maxcumsize']))
                continue

        if b.cumsize<kwargs['mincumsize']:
            logging.debug("Skipping bubble <%s,%s>, cumulative size %d is smaller than mincumsize=%d."%(b.source,b.sink,b.cumsize,kwargs['mincumsize']))
            continue

        if len(b.nodes)==3:
            logging.debug("Skipping bubble <%s,%s>, indel, no point in realigning."%(b.source,b.sink))
            continue

        b.G=None #remove reference to Graph
        realignbubbles.append(b)
    
    logging.info("Done.")

    if len(realignbubbles)==0:
        logging.info("No bubbles that qualify for realignment.")
    else:
        realignbubbles.sort(key=lambda b: b.source_idx)
        distinctbubbles=[realignbubbles[0]]
        p=0
        i=1
        for i in xrange(i,len(realignbubbles)):
            if realignbubbles[i].source_idx >= realignbubbles[p].sink_idx:
                distinctbubbles.append(realignbubbles[i])
                p=i
            else:
                logging.debug("Skipping realignment for: <%s,%s> - is contained in <%s,%s>"%(realignbubbles[i].source, realignbubbles[i].sink, realignbubbles[p].source, realignbubbles[p].sink))

        logging.info("Realigning a total of %d bubbles"%len(distinctbubbles))
        nn=max([node for node in G.nodes() if type(node)==int])+1

        if kwargs['nproc']>1:
            # inputq = mp.Queue()

            nworkers=kwargs['nproc']-1

            outputq = mp.Queue(nworkers*2)

            aworkers=[]

            shuffle(distinctbubbles) #make sure not all the big telomeric bubbles end up with one worker

            if nworkers>len(distinctbubbles):
                nworkers=len(distinctbubbles)

            # chunksize=50
            # i=0
            # while (i*chunksize)<len(distinctbubbles):
            #     chunk=distinctbubbles[(i*chunksize):((i+1)*chunksize)]
            #     print "Putting chunk:",i,len(chunk)
            #     inputq.put( (i,chunk) ) #,False
            #     i+=1

            chunksize=int(math.floor(len(distinctbubbles)/float(nworkers)))

            for i in range(nworkers):
                logging.info("Starting worker: %d"%i)
                
                if i==nworkers-1:
                    chunk=distinctbubbles[(i*chunksize):]
                else:
                    chunk=distinctbubbles[(i*chunksize):(i*chunksize)+chunksize]

                #create a subgraph for this chunk such that the worker doesnt need to load the entire graph
                t0=time.time()
                Gs=G.subgraph([node for bubble in chunk for node in bubble.nodes]).copy()
                t1=time.time()
                logging.info("Created subgraph in %.2f seconds."%(t1-t0))

                p=mp.Process(target=align_worker, args=(Gs,chunk,outputq,kwargs))
                # p=mp.Process(target=align_worker, args=(G,inputq,outputq,kwargs))
                # p=mp.Process(target=align_worker, args=(G,i,nworkers,distinctbubbles,outputq,kwargs))
                aworkers.append(p) #(p,align_worker,(G,inputq,outputq,kwargs)))

            try:
                for p in aworkers:
                    p.start()

                graph_worker(G,nn,outputq,aworkers,len(distinctbubbles))
                
            except Exception, e:
                logging.fatal("%s"%str(e))
                # for p,fn,args in aworkers:
                for p in aworkers:
                    p.terminate()
                outputq.close()
                exit(1)
            
            outputq.close()
            # inputq.close()

            logging.info("Waiting for workers to finish...")

            # for p,fn,args in aworkers:
            for p in aworkers:
                p.join()

            logging.info("Done.")

        else:
            for bubble in distinctbubbles:
                G.node[bubble.source]['aligned']=1
                G.node[bubble.sink]['aligned']=1

                bnodes=list(set(bubble.nodes)-set([bubble.source,bubble.sink]))
                
                t0=time.time()
                sg=G.subgraph(bnodes).copy()
                t1=time.time()
                logging.info("Extract subgraph for: <%s,%s> took %.4f seconds."%(bubble.source,bubble.sink,t1-t0))
                
                offsets=dict()
                for sid in G.node[bubble.source]['offsets']:
                    offsets[sid]=G.node[bubble.source]['offsets'][sid]+len(G.node[bubble.source]['seq'])

                sourcesamples=set(G.node[bubble.source]['offsets'].keys())
                sinksamples=set(G.node[bubble.sink]['offsets'].keys())
                paths=sourcesamples.intersection(sinksamples)

                # t0=time.time()
                res=refine_bubble(sg,bubble,offsets,paths, **kwargs)
                # t1=time.time()
                # logging.info("Refining bubble: <%s,%s> took %.4f seconds."%(bubble.source,bubble.sink,t1-t0))

                if res==None:
                    continue
                else:
                    bubble,ng,path2start,path2end=res
                    # G,nn=replace_bubble(G,bubble,ng,path2start,path2end,nn)
                    # t0=time.time()
                    G,nn=replace_bubble(G,bubble,ng,path2start,path2end,nn)
                    # t1=time.time()
                    # logging.info("Replacing bubble: <%s,%s> took %.4f seconds."%(bubble.source,bubble.sink,t1-t0))

    return G

def msa2graph(aobjs,idoffset=0,msa='muscle',parameters="",minconf=0,constrans=2,consgap=True,nrefinements=100):

    nn=idoffset
    ng=nx.DiGraph()
    ng.graph['paths']=[]
    ng.graph['path2id']=dict()
    ng.graph['id2path']=dict()
    ng.graph['id2end']=dict()

    maxsize=0
    for name,seq in aobjs:
        sid=len(ng.graph['paths'])
        ng.graph['path2id'][name]=sid
        ng.graph['id2path'][sid]=name
        ng.graph['id2end'][sid]=len(seq)
        ng.graph['paths'].append(name)
        if len(seq)>maxsize:
            maxsize=len(seq)

    uid=uuid.uuid4().hex
    tempfiles=[]

    if msa in {'muscle','pecan','msaprobs','probcons'}:
        logging.debug("Trying to construct MSA with %s, minconf=%d."%(msa,minconf))

        if msa=='muscle':
            cmd="muscle -in %s.fasta -quiet"%uid
            fasta_writer(uid+".fasta",aobjs)
            tempfiles.append("%s.fasta"%uid)
        elif msa=='probcons':
            # cmd="probcons %s.fasta -pre 1 -annot %s.conf"%(uid,uid)
            cmd="probcons %s.fasta -annot %s.conf %s"%(uid,uid,parameters) #-p /Users/jasperlinthorst/Documents/phd/probcons/nw.txt
            fasta_writer(uid+".fasta",aobjs)
            tempfiles.append("%s.fasta"%uid)
            tempfiles.append("%s.conf"%uid)
        elif msa=='pecan':
            cmd="pecan -G %s.fasta -F %s.*.fasta -l -p %s.conf %s && cat %s.fasta"%(uid,uid,uid,parameters,uid)
            for i,(name,seq) in enumerate(aobjs): #pecan wants sequence in separate files
                fasta_writer("%s.%d.fasta"%(uid,i),[(name,seq)])
                tempfiles.append("%s.%d.fasta"%(uid,i))
            tempfiles.append("%s.fasta"%uid)
            tempfiles.append("%s.conf"%uid)
        elif msa=='msaprobs':
            cmd="msaprobs %s.fasta -annot %s.conf %s"%(uid,uid,parameters)
            fasta_writer(uid+".fasta",aobjs)
            tempfiles.append("%s.fasta"%uid)
            tempfiles.append("%s.conf"%uid)
        else:
            logging.fatal("Unkown multiple sequence aligner: %s"%msa)
            sys.exit(1)
        
        seqs=[""]*len(aobjs)
        names=[""]*len(aobjs)

        try:
            DEVNULL = open(os.devnull, 'wb')
            for a in subprocess.check_output([cmd],shell=True,stderr=DEVNULL).split(">")[1:]:
                x=a.find('\n')
                name=a[:x]
                seq=a[x+1:].replace("\n","")
                names[ng.graph['path2id'][name]]=name
                seqs[ng.graph['path2id'][name]]=seq
        except Exception as e:
            logging.fatal("System call to %s failed: \"%s\""%(msa,e.output))
            return

        confidence=[100]*len(seq) #initialize to 100% accuracy for each column

        if os.path.exists("%s.conf"%uid): #if there's an annotation file that accompanies the msa
            with open("%s.conf"%uid) as annot:
                for i,line in enumerate(annot):
                    confidence[i]=float(line.strip()) #expected percentage of correct pairwise matches in the i'th column of the msa...
                    if confidence[i]<1: #consider it a ratio, otherwise a percentage
                        confidence[i]=confidence[i]*100

    else:
        pl=probconslib.probcons()
        t0=time.time()
        aln=pl.align(aobjs,consistency=constrans,refinement=nrefinements,pretraining=0,consgap=consgap)
        t1=time.time()
        logging.debug("ProbCons MSA took %.4f seconds for %d alleles with maxsize=%d."%((t1-t0),len(aobjs),maxsize))

        seqs=[""]*len(aobjs)
        names=[""]*len(aobjs)
        for name,seq in aln[0]:
            names[ng.graph['path2id'][name]]=name
            seqs[ng.graph['path2id'][name]]=seq
        confidence=aln[1]

    for i,seq in enumerate(seqs):
        # if len(seq)>200:
        #     logging.debug("OUT %s: %s...%s"%(str(i).rjust(4, ' '),seq[0:100],seq[-100:]))
        #     logging.debug("CONF    : %s...%s"%("".join([str(c/10) for c in confidence[:100]]),"".join([str(c/10) for c in confidence[-100:]])))
        # else:
        logging.debug("OUT %s: %s"%(str(i).rjust(4, ' '),seq))
    
    logging.debug("CONF    : %s"%"".join([str(c/10) for c in confidence]))
    


    offsets={o:-1 for o in range(len(seqs))}
    nid=nn
    for i in xrange(len(seqs[0])):
        col={}
        base2node={}
        sid2node={}
        p=confidence[i]

        for j in xrange(len(seqs)):
            if seqs[j][i] in col:
                col[seqs[j][i]].add(j)
            else:
                col[seqs[j][i]]=set([j])
            if seqs[j][i]!='-':
                offsets[j]+=1

        for base in col:
            if i==0:
                assert(len(col[base])>0)
                # if len(col[base])>0:
                if p>=minconf:
                    ng.add_node(nid,seq=base,offsets={sid:offsets[sid] for sid in offsets if sid in col[base]},p=[p])
                    base2node[base]=nid
                    for sid in col[base]:
                        assert(sid not in sid2node)
                        sid2node[sid]=nid
                    nid+=1
                else: #new node per sequence
                    for sid in col[base]:
                        ng.add_node(nid,seq=base,offsets={sid:offsets[sid]},p=[p])
                        assert(sid not in sid2node)
                        sid2node[sid]=nid
                        if base in base2node:
                            base2node[base].append(nid)
                        else:
                            base2node[base]=[nid]
                        nid+=1
            else:

                if p>=minconf and pp>=minconf:
                    for pbase in pcol:
                        overlap=pcol[pbase].intersection(col[base])
                        if len(overlap)==0:
                            continue
                        elif len(overlap)==len(col[base])==len(pcol[pbase]): #append seq
                            ng.node[pbase2node[pbase]]['seq']+=base
                            ng.node[pbase2node[pbase]]['p']+=[p]
                            
                            base2node[base]=pbase2node[pbase]
                            
                            for sid in overlap:
                                assert(sid not in sid2node)
                                sid2node[sid]=sid2pnode[sid]
                        else:
                            assert(len(overlap)>0)
                            if base not in base2node: #if not already there
                                ng.add_node(nid,seq=base,offsets=dict(),p=[p])
                                base2node[base]=nid
                                for sid in col[base]:
                                    assert(sid not in sid2node)
                                    sid2node[sid]=nid
                                nid+=1
                            for sid in overlap:
                                ng.node[base2node[base]]['offsets'][sid]=offsets[sid]

                            ng.add_edge(pbase2node[pbase],base2node[base],paths=overlap,oto='+',ofrom='+')

                elif p<minconf and pp>=minconf:

                    for sid in col[base]:
                        ng.add_node(nid,seq=base,offsets={sid:offsets[sid]},p=[p])
                        ng.add_edge(sid2pnode[sid],nid,paths={sid},oto='+',ofrom='+')
                        sid2node[sid]=nid
                        if base in base2node:
                            base2node[base].append(nid)
                        else:
                            base2node[base]=[nid]
                        nid+=1

                elif p>=minconf and pp<minconf:
                    
                    ng.add_node(nid,seq=base,offsets=dict(),p=[p])
                    for sid in col[base]:
                        ng.node[nid]['offsets'][sid]=offsets[sid]
                        if not ng.has_edge(sid2pnode[sid],nid):
                            ng.add_edge(sid2pnode[sid],nid,paths={sid},oto='+',ofrom='+')
                        else:
                            ng[sid2pnode[sid]][nid]['paths'].add(sid)
                        sid2node[sid]=nid
                        base2node[base]=nid
                    nid+=1

                elif p<minconf and pp<minconf:
                    for sid in col[base]:
                        ng.node[sid2pnode[sid]]['seq']+=base
                        ng.node[sid2pnode[sid]]['p'].append(p)
                    sid2node=sid2pnode

                else:
                    logging.error("Impossible combination!")
                    sys.exit(1)

        assert(len(sid2node)==len(seqs))
        sid2pnode=sid2node
        pbase2node=base2node
        pcol=col
        pp=p

    #remove gaps from graph
    remove=[]
    for node,data in ng.nodes(data=True):
        incroffset=False
        if data['seq'][0]=='-':
            incroffset=True

        data['seq']=data['seq'].replace("-","")
        if data['seq']=="":
            remove.append(node)
        elif incroffset:
            for sid in data['offsets']:
                data['offsets'][sid]+=1

        if len(data['offsets'])>1:
            data['aligned']=1
        else:
            data['aligned']=0

    for node in remove:
        ine=ng.in_edges(node,data=True)
        oute=ng.out_edges(node,data=True)
        for in1,in2,ind in ine:
            for out1,out2,outd in oute:
                overlap=ind['paths'].intersection(outd['paths'])
                if len(overlap)>=1:
                    if ng.has_edge(in1,out2):
                        ng[in1][out2]['paths']=ng[in1][out2]['paths'] | overlap
                    else:
                        ng.add_edge(in1,out2,paths=overlap,ofrom='+',oto='+')

    ng.remove_nodes_from(remove)

    #contract edges
    updated=True
    while updated:
        updated=False
        for v,t in ng.edges():
            if ng.out_degree(v)==ng.in_degree(t)==1:
                if ng.node[v]['offsets'].keys()==ng.node[t]['offsets'].keys():
                    ng.node[v]['seq']+=ng.node[t]['seq']
                    for suc in ng.successors(t):
                        ng.add_edge(v,suc,**ng[t][suc])
                    ng.remove_node(t)
                    updated=True
                    break

    for tmpfile in tempfiles:
        try:
            os.remove(tmpfile)
        except Exception as e:
            logging.fatal("Failed to remove tmp file: \"%s\""%tmpfile)
            return

    # logging.debug("%d nodes in refined graph:",ng.number_of_nodes())
    # for node in ng:
        # logging.debug("%s: %s"%(node,ng.node[node]['seq']))

    return ng