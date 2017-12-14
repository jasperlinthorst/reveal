from utils import *
from extract import extract
from align import align,prune_nodes
import bubbles
import schemes

def seq2node(G,T,toupper=True):
    for node in G:
        if isinstance(node,Interval):
            if toupper:
                G.node[node]['seq']=T[node.begin:node.end].upper()
            else:
                G.node[node]['seq']=T[node.begin:node.end]

def realign_bubble_cmd(args):
    if len(args.graph)<1:
        logging.fatal("Specify a gfa file for which bubbles should be realigned.")
        return
    
    G=nx.DiGraph()
    read_gfa(args.graph[0],None,"",G)
    
    if (args.source==None and args.sink==None) and (args.all or args.complex or args.simple):
        G=realign_all(G,
                            minlength=args.minlength,
                            minn=args.minn,
                            wscore=args.wscore,
                            wpen=args.wpen,
                            maxlen=args.maxlen,
                            minlen=args.minlen,
                            seedsize=args.seedsize,
                            maxmums=args.maxmums,
                            gcmodel=args.gcmodel,
                            complex=args.complex,
                            simple=args.simple,
                            all=args.all,
                            minsize=args.minsize,
                            muscle=args.muscle,
                            sa64=args.sa64)
    else:
        if args.source==None or args.sink==None:
            logging.error("Specify source sink pair")
            sys.exit(1)

        b=bubbles.Bubble(G,args.source,args.sink)

        G=realign_bubble(G,b,
                                minlength=args.minlength,
                                minn=args.minn,
                                wscore=args.wscore,
                                wpen=args.wpen,
                                maxlen=args.maxlen,
                                seedsize=args.seedsize,
                                maxmums=args.maxmums,
                                gcmodel=args.gcmodel,
                                muscle=args.muscle,
                                sa64=args.sa64)
    
    if args.outfile==None:
        fn=args.graph[0].replace(".gfa",".realigned.gfa")
    else:
        fn=args.outfile
    
    write_gfa(G,"",outputfile=fn)

def realign_bubble(G,bubble,**kwargs):
    
    nn=max(G.nodes())+1

    source=bubble.source
    sink=bubble.sink
    
    assert(source in G)
    assert(sink in G)

    sourcesamples=set(G.node[source]['offsets'].keys())
    sinksamples=set(G.node[sink]['offsets'].keys())
    
    if sourcesamples!=sinksamples:
        logging.error("Specify proper source/sink pair.")
        return G
    
    if bubble.seqsize>kwargs['maxlen']:
        logging.fatal("Bubble (%s,%s) is too big. Increase --maxlen."%(source,sink))
        return G

    if len(bubble.nodes)==3:
        logging.fatal("Indel bubble, no point realigning.")
        return G
    
    bubblenodes=bubble.nodes[1:-1]

    sg=G.subgraph(bubblenodes)
    d={}
    aobjs=[]
    
    #extract all paths
    for sid in sourcesamples.intersection(sinksamples): #should be equal..
        seq=extract(sg,G.graph['id2path'][sid])
        if len(seq)>0:
            aobjs.append((G.graph['id2path'][sid],seq))

    for name,seq in aobjs:
        logging.debug("IN %s: %s%s"%(name,seq[:200],'...'if len(seq)>200 else ''))

    if kwargs['muscle']: #use muscle multiple sequence aligner to align bubble
        ng=msa2graph(aobjs,idoffset=nn,msa="muscle")

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
        prune_nodes(ng,T)
        seq2node(ng,T) #transfer sequence to node attributes

    #map edge atts back to original graph
    for n1,n2,data in ng.edges(data=True):
        old=data['paths']
        new=set()
        for sid in old:
            new.add( G.graph['path2id'][ng.graph['id2path'][sid]] )
        data['paths']=new

    #map node atts back to original graph
    for node,data in ng.nodes(data=True):
        old=data['offsets']
        new=dict()
        for sid in old:
            new[G.graph['path2id'][ng.graph['id2path'][sid]]]=old[sid]
        data['offsets']=new

    ng.graph['paths']=G.graph['paths']
    ng.graph['path2id']=G.graph['path2id']
    ng.graph['id2path']=G.graph['id2path']

    for sample in G.graph['paths']:
        assert(G.graph['path2id'][sample]==ng.graph['path2id'][sample])

    mapping={}
    
    path2start=dict()
    path2end=dict()

    #map nodes back to original offsets and idspace, and determine first/last node for every path
    for node,data in ng.nodes(data=True):
        for sid in data['offsets']:
            if sid not in path2start or data['offsets'][sid]<path2start[sid][1]:
                path2start[sid]=(nn,data['offsets'][sid])

        for sid in data['offsets']:
            if sid not in path2end or data['offsets'][sid]>path2end[sid][1]:
                path2end[sid]=(nn,data['offsets'][sid])

        corrected=dict()
        for sid in data['offsets']:
            # corrected[sid]=data['offsets'][sid]+G.node[source]['offsets'][sid]
            corrected[sid]=data['offsets'][sid]+(G.node[source]['offsets'][sid]+len(G.node[source]['seq']))

        ng.node[node]['offsets']=corrected
        mapping[node]=nn
        nn+=1

    startnodes=set([p[0] for p in path2start.values()])
    endnodes=set([p[0] for p in path2end.values()])
    ng=nx.relabel_nodes(ng,mapping)
    
    for node in bubblenodes: #remove all bubblenodes from the original graph
        G.remove_node(node)

    for node,data in ng.nodes(data=True): #add nodes from newly aligned graph to original graph
        G.add_node(node,data)
    
    for edge in ng.edges(data=True):
        G.add_edge(*edge)

    for start in startnodes:
        x=set(G.node[source]['offsets'].keys()).intersection(set(G.node[start]['offsets'].keys()))
        if len(x)>0:
            G.add_edge(source,start,ofrom='+',oto='+',paths=x)

    for end in endnodes:
        x=set(G.node[sink]['offsets'].keys()).intersection(set(G.node[end]['offsets'].keys()))
        if len(x)>0:
            G.add_edge(end,sink,ofrom='+',oto='+',paths=x)
    
    return G

def realign_all(G,  **kwargs):
    realignbubbles=[]
    
    if kwargs['minsize']==None:
        kwargs['minsize']=kwargs['minlength']

    #detect all bubbles
    for b in bubbles.bubbles(G):

        if kwargs['complex']:
            if b.issimple():
                logging.debug("Skipping bubble %s, not complex."%str(b.nodes))
                continue

        if kwargs['simple']:
            if not b.issimple():
                logging.debug("Skipping bubble %s, not simple."%str(b.nodes))
                continue

        if b.minsize<kwargs['minsize']:
            logging.debug("Skipping bubble %s, smaller %d than minsize=%d."%(str(b.nodes),b.minsize,kwargs['minsize']))
            continue

        if b.seqsize>kwargs['maxlen']:
            logging.debug("Skipping bubble %s, larger %d than maxlen=%d."%(str(b.nodes),b.seqsize,kwargs['maxlen']))
            continue

        if b.seqsize<kwargs['minlen']:
            logging.debug("Skipping bubble %s, smaller %d than minlen=%d."%(str(b.nodes),b.seqsize,kwargs['minlen']))
            continue

        if len(b.nodes)==3:
            logging.debug("Skipping bubble %s, indel, no point in realigning."%(str(b.nodes)))
            continue

        sourcesamples=set(G.node[b.source]['offsets'].keys())
        sinksamples=set(G.node[b.sink]['offsets'].keys())

        if sourcesamples!=sinksamples:
            logging.warn("Invalid bubble between %s and %s"%(str(b.source),str(b.sink)))
            continue
        
        realignbubbles.append(b)

    distinctbubbles=[]
    for b1 in realignbubbles:
        # nodes=realignbubbles[(source,sink)]
        for b2 in realignbubbles:
            if set(b2.nodes).issuperset(set(b1.nodes)) and not set(b1.nodes)==set(b2.nodes):
                logging.debug("Skipping bubble %s, because its nested in %s."%(str(b1.nodes),str(b2.nodes)))
                break
        else:
            distinctbubbles.append(b1)

    logging.info("Realigning a total of %d bubbles"%len(distinctbubbles))
    i=1
    for bubble in distinctbubbles:
        logging.info("Realigning bubble between <%s> and <%s>, cumulative size %dbp (in nodes=%d)."%(bubble.source,bubble.sink,bubble.seqsize,len(bubble.nodes)-2))
        G=realign_bubble(G,bubble, **kwargs)
        i+=1
    
    return G