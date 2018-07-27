from utils import *
from extract import extract,extract_path
from rem import align,prune_nodes
from random import shuffle
import bubbles
import schemes
from multiprocessing.pool import Pool
from multiprocessing import Process,Queue
import signal
import probconslib
import math


def refine_bubble_cmd(args):
    if len(args.graph)<1:
        logging.fatal("Specify a gfa file for which bubbles should be realigned.")
        return
    
    # G=nx.MultiDiGraph() #TODO: make sure that refine can handle structural variant edges, so make sure we use a MultiDiGraph here!
    G=nx.DiGraph()

    read_gfa(args.graph[0],None,"",G)
    
    logging.info("Paths through the graph: %s"%G.graph['paths'])

    if (args.source==None and args.sink==None) and (args.all or args.complex or args.simple):
        G=refine_all(G,
                            minlength=args.minlength,
                            minn=args.minn,
                            wscore=args.wscore,
                            wpen=args.wpen,
                            maxsize=args.maxsize,
                            minmaxsize=args.minmaxsize,
                            minsize=args.minsize,
                            maxcumsize=args.maxcumsize,
                            mincumsize=args.mincumsize,
                            seedsize=args.seedsize,
                            maxmums=args.maxmums,
                            gcmodel=args.gcmodel,
                            complex=args.complex,
                            simple=args.simple,
                            nogaps=args.nogaps,
                            all=args.all,
                            method=args.method,
                            parameters=args.parameters,
                            minconf=args.minconf,
                            nproc=args.nproc,
                            sa64=args.sa64,
                            constrans=args.constrans,
                            nrefinements=args.nrefinements,
                            consgap=args.consgap
                            )
    else:
        if args.source==None or args.sink==None:
            logging.error("Specify source sink pair")
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

        res=refine_bubble(sg,b,offsets,paths,
                                minlength=args.minlength,
                                minn=args.minn,
                                wscore=args.wscore,
                                wpen=args.wpen,
                                maxsize=args.maxsize,
                                minmaxsize=args.minmaxsize,
                                seedsize=args.seedsize,
                                maxmums=args.maxmums,
                                gcmodel=args.gcmodel,
                                method=args.method,
                                parameters=args.parameters,
                                minconf=args.minconf,
                                sa64=args.sa64,
                                constrans=args.constrans,
                                nrefinements=args.nrefinements,
                                consgap=args.consgap
                                )

        if res!=None:
            bubble,ng,path2start,path2end=res
            G,nn=replace_bubble(G,b,ng,path2start,path2end,nn)

    if args.outfile==None:
        fn=args.graph[0].replace(".gfa",".realigned.gfa")
    else:
        fn=args.outfile

    prune_nodes(G)

    contract(G,[n for n in nx.topological_sort(G) if type(n)!=str])

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

    d={}
    aobjs=[]

    #extract all paths
    for sid in paths:
        seq=extract(sg,sg.graph['id2path'][sid])
        if len(seq)>0:
            if seq in d:
                d[seq].append(str(sid))
            else:
                d[seq]=[str(sid)]

    if len(d)==1:
        logging.info("Nothing to refine for bubble: %s - %s"%(bubble.source,bubble.sink))
        return

    aobjs=[(",".join(d[seq]),seq) for seq in d]

    logging.info("Realigning bubble between <%s> and <%s>, %d alleles, with %s (max size %dbp, in nodes=%d)."%(bubble.source,bubble.sink,len(aobjs),kwargs['method'],bubble.maxsize,len(bubble.nodes)-2))

    for name,seq in aobjs:
        if len(seq)>200:
            logging.debug("IN %s: %s...%s (%d bp)"%(name.rjust(4,' '),seq[:100],seq[-100:],len(seq)))
        else:
            logging.debug("IN %s: %s (%d bp)"%(name.rjust(4,' '),seq,len(seq)))

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

def align_worker(G,chunk,outputq,kwargs,chunksize=500):
    
    try:
        rchunk=[]
        for bubble in chunk:
            G.node[bubble.source]['aligned']=1
            G.node[bubble.sink]['aligned']=1

            logging.debug("Submitting realign bubble between <%s> and <%s>, max allele size %dbp (in nodes=%d)."%(bubble.source,bubble.sink,bubble.maxsize,len(bubble.nodes)-2))
            bnodes=list(set(bubble.nodes)-set([bubble.source,bubble.sink]))
            sg=G.subgraph(bnodes).copy()
            
            offsets=dict()
            for sid in G.node[bubble.source]['offsets']:
                offsets[sid]=G.node[bubble.source]['offsets'][sid]+len(G.node[bubble.source]['seq'])

            sourcesamples=set(G.node[bubble.source]['offsets'].keys())
            sinksamples=set(G.node[bubble.sink]['offsets'].keys())
            paths=sourcesamples.intersection(sinksamples)

            b=refine_bubble(sg,bubble,offsets,paths,**kwargs)
            
            if b==None:
                continue
            else:
                rchunk.append(b)
                if len(rchunk)==chunksize:
                    outputq.put(rchunk)
                    rchunk=[]

        if len(rchunk)>0:
            outputq.put(rchunk)

        outputq.put(-1)
    except:
        logging.fatal("Worker failed!")
        outputq.put(-1)
    
def graph_worker(G,nn,outputq,nworkers):
    deadworkers=0
    while True:
        data = outputq.get()
        if data==-1:
            deadworkers+=1
            if deadworkers==nworkers:
                break
        else:
            for d in data:
                bubble,ng,path2start,path2end=d
                logging.info("Replacing bubble: %s"%bubble.nodes)
                G,nn=replace_bubble(G,bubble,ng,path2start,path2end,nn)

def refine_all(G, **kwargs):
    realignbubbles=[]
    
    if kwargs['minsize']==None:
        kwargs['minsize']=kwargs['minlength']

    #detect all bubbles
    for b in bubbles.bubbles(G):

        if kwargs['complex']:
            if b.issimple():
                logging.debug("Skipping bubble %s, not complex."%str(b.nodes))
                continue

        if kwargs['nogaps']:
            for n in b.nodes:
                if 'N' in G.nodes[n]['seq']:
                    logging.info("Skipping bubble %s, bubble spans a gap."%str(b.nodes))
                    continue

        if kwargs['simple']:
            if not b.issimple():
                logging.debug("Skipping bubble %s, not simple."%str(b.nodes))
                continue

        if b.minsize<kwargs['minsize']:
            logging.debug("Skipping bubble %s, smallest allele (%dbp) is smaller than minsize=%d."%(str(b.nodes),b.minsize,kwargs['minsize']))
            continue

        if b.maxsize<kwargs['minmaxsize']:
            logging.debug("Skipping bubble %s, largest allele (%dbp) is smaller than minmaxsize=%d."%(str(b.nodes),b.maxsize,kwargs['minmaxsize']))
            continue

        if b.maxsize>kwargs['maxsize']:
            logging.warn("Skipping bubble %s, largest allele (%dbp) is larger than maxsize=%d."%(str(b.nodes),b.maxsize,kwargs['maxsize']))
            continue

        if kwargs['maxcumsize']!=None:
            if b.cumsize>kwargs['maxcumsize']:
                logging.warn("Skipping bubble %s, cumulative size %d is larger than maxcumsize=%d."%(str(b.nodes),b.cumsize,kwargs['maxcumsize']))
                continue

        if b.cumsize<kwargs['mincumsize']:
            logging.debug("Skipping bubble %s, cumulative size %d is smaller than mincumsize=%d."%(str(b.nodes),b.cumsize,kwargs['mincumsize']))
            continue

        if len(b.nodes)==3:
            logging.debug("Skipping bubble %s, indel, no point in realigning."%(str(b.nodes)))
            continue

        realignbubbles.append(b)

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
    
    logging.info("Done.")

    logging.info("Realigning a total of %d bubbles"%len(distinctbubbles))
    nn=max([node for node in G.nodes() if type(node)==int])+1

    if kwargs['nproc']>1:
        outputq = Queue()
        nworkers=kwargs['nproc']
        aworkers=[]

        shuffle(distinctbubbles) #make sure not all the big telomeric bubbles and up with one worker

        if nworkers>len(distinctbubbles):
            nworkers=len(distinctbubbles)

        chunksize=int(math.ceil(len(distinctbubbles)/float(nworkers)))
        for i in range(0, len(distinctbubbles), chunksize):
            p=Process(target=align_worker, args=(G,distinctbubbles[i:i+chunksize],outputq,kwargs))
            p.start()
            aworkers.append(p)

        try:
            graph_worker(G,nn,outputq,nworkers)
        except:
            logging.fatal("Failed to update graph with refined bubble! Signal workers to stop.")
            for p in aworkers:
                p.terminate()
            outputq.close()
            exit(1)
        
        outputq.close()
        logging.info("Waiting for workers to finish...")
        for p in aworkers:
            p.join()
        logging.info("Done.")

    else:
        for bubble in distinctbubbles:
            G.node[bubble.source]['aligned']=1
            G.node[bubble.sink]['aligned']=1

            bnodes=list(set(bubble.nodes)-set([bubble.source,bubble.sink]))
            sg=G.subgraph(bnodes)
            
            offsets=dict()
            for sid in G.node[bubble.source]['offsets']:
                offsets[sid]=G.node[bubble.source]['offsets'][sid]+len(G.node[bubble.source]['seq'])

            sourcesamples=set(G.node[bubble.source]['offsets'].keys())
            sinksamples=set(G.node[bubble.sink]['offsets'].keys())
            paths=sourcesamples.intersection(sinksamples)

            res=refine_bubble(sg,bubble,offsets,paths, **kwargs)
            if res==None:
                continue
            else:
                bubble,ng,path2start,path2end=res
                G,nn=replace_bubble(G,bubble,ng,path2start,path2end,nn)

    return G

def msa2graph(aobjs,idoffset=0,msa='muscle',parameters="",minconf=0,constrans=2,consgap=True,nrefinements=100):

    nn=idoffset
    ng=nx.DiGraph()
    ng.graph['paths']=[]
    ng.graph['path2id']=dict()
    ng.graph['id2path']=dict()
    ng.graph['id2end']=dict()

    for name,seq in aobjs:
        sid=len(ng.graph['paths'])
        ng.graph['path2id'][name]=sid
        ng.graph['id2path'][sid]=name
        ng.graph['id2end'][sid]=len(seq)
        ng.graph['paths'].append(name)

    #TODO: writing to a temporary file (for now), but this should ideally happen in memory
    uid=uuid.uuid4().hex
    tempfiles=[]
    logging.debug("Trying to construct MSA with %s, minconf=%d."%(msa,minconf))

    if msa in {'muscle','pecan','msaprobs','probcons'}:

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
            cmd="java -cp /Users/jasperlinthorst/Documents/phd/pecan bp/pecan/Pecan -G %s.fasta -F %s.*.fasta -l -p %s.conf %s && cat %s.fasta"%(uid,uid,uid,parameters,uid)
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
        logging.debug("Using probcons (in memory)")
        pl=probconslib.probcons()
        aln=pl.align(aobjs,consistency=constrans,refinement=nrefinements,pretraining=0,consgap=consgap)
        seqs=[""]*len(aobjs)
        names=[""]*len(aobjs)
        for name,seq in aln[0]:
            names[ng.graph['path2id'][name]]=name
            seqs[ng.graph['path2id'][name]]=seq
        confidence=aln[1]

        for i,seq in enumerate(seqs):
            if len(seq)>200:
                logging.debug("OUT %s: %s...%s"%(str(i).rjust(4, ' '),seq[0:100],seq[-100:]))
                logging.debug("CONF    : %s...%s"%("".join([str(c/10) for c in confidence[:100]]),"".join([str(c/10) for c in confidence[-100:]])))
            else:
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

    logging.debug("%d nodes in refined graph:",ng.number_of_nodes())
    for node in ng:
        logging.debug("%s: %s"%(node,ng.node[node]['seq']))

    return ng