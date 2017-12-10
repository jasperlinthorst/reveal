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
    
    if args.all or args.complex:
        G=realign_all(G,    minlength=args.minlength,
                            minn=args.minn,
                            wscore=args.wscore,
                            wpen=args.wpen,
                            maxlen=args.maxlen,
                            seedsize=args.seedsize,
                            maxmums=args.maxmums,
                            gcmodel=args.gcmodel,
                            complex=args.complex,
                            minsize=args.minsize,
                            sa64=args.sa64)
    else:
        if args.source==None or args.sink==None:
            logging.error("Specify source sink pair")
            sys.exit(1)

        G=realign_bubble(G,args.source,args.sink,minlength=args.minlength,
                                                 minn=args.minn,
                                                 wscore=args.wscore,
                                                 wpen=args.wpen,
                                                 maxlen=args.maxlen,
                                                 seedsize=args.seedsize,
                                                 maxmums=args.maxmums,
                                                 gcmodel=args.gcmodel,
                                                 sa64=args.sa64)
    
    if args.outfile==None:
        fn=args.graph[0].replace(".gfa",".realigned.gfa")
    else:
        fn=args.outfile
    
    write_gfa(G,"",outputfile=fn)

def realign_bubble(G,source,sink,minlength=20,
                                 minn=2,
                                 maxlen=10000000,
                                 wscore=1,
                                 wpen=1,
                                 seedsize=None,
                                 maxmums=None,
                                 gcmodel="sumofpairs",
                                 sa64=False):

    nn=max(G.nodes())+1
    bubblenodes=[]
    
    assert(source in G)
    assert(sink in G)
    sourcesamples=set(G.node[source]['offsets'].keys())
    sinksamples=set(G.node[sink]['offsets'].keys())
    
    if sourcesamples!=sinksamples:
        logging.error("Specify proper source/sink pair.")
        sys.exit(1)
    
    bubblesnodes=[]
    add=False
    for node in nx.topological_sort(G):
        if node==source:
            add=True
        if add:
            bubblenodes.append(node) #TOOD: this comes from bubble class now..
        if node==sink:
            add=False
    
    sg=G.subgraph(bubblenodes)
    d={}
    aobjs=[]
    cumsum=0
    
    #extract all paths
    for sid in sourcesamples.intersection(sinksamples): #should be equal..
        seq=extract(sg,G.graph['id2path'][sid])

        cumsum+=len(seq)
        if len(seq)>0:
            aobjs.append((G.graph['id2path'][sid],seq))

        if cumsum>maxlen:
            logging.fatal("Bubble (%s,%s) is too big. Increase --maxlen."%(source,sink))
            sys.exit(1)

    ng,idx=align(aobjs, minlength=minlength,
                        minn=minn,
                        seedsize=seedsize,
                        maxmums=maxmums,
                        wpen=wpen,
                        wscore=wscore,
                        gcmodel=gcmodel,
                        sa64=sa64)
    T=idx.T
    
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

    prune_nodes(ng,T)
    
    seq2node(ng,T) #transfer sequence to node attributes

    mapping={}
    startnodes=set()
    endnodes=set()
    
    #map nodes back to original offsets and ids
    for node,data in ng.nodes(data=True):

        if len(ng.predecessors(node))==0:
            startnodes.add(nn)

        if len(ng.successors(node))==0:
            endnodes.add(nn)

        corrected=dict()
        for sid in data['offsets']:
            corrected[sid]=data['offsets'][sid]+G.node[source]['offsets'][sid]

        ng.node[node]['offsets']=corrected
        mapping[node]=nn
        nn+=1

    ng=nx.relabel_nodes(ng,mapping) #relabel nodes to get rid of Intervals
    
    #store edges that have to be reconnected after removing
    pre=G.predecessors(source)
    suc=G.successors(sink)
    
    for node in bubblenodes: #remove all bubblenodes from the original graph
        G.remove_node(node)

    for node,data in ng.nodes(data=True): #add nodes from newly aligned graph to original graph, except for source and sink, update those with seq attribute
        assert(node not in G.node)
        G.add_node(node,data)
    
    for edge in ng.edges(data=True):
        assert(edge[0] in G)
        assert(edge[1] in G)
        assert(edge[1] not in G[edge[0]])
        G.add_edge(*edge)
    
    #reconnect nodes
    for p in pre:
        for start in startnodes:
            x=set(G.node[p]['offsets'].keys()).intersection(set(G.node[start]['offsets'].keys()))
            if len(x)!=0:
                G.add_edge(p,start,ofrom='+',oto='+',paths=x )
            else:
                print "!",source, sink, p, start, G.node[p]['offsets'], G.node[start]['offsets']
    for s in suc:
        for end in endnodes:
            x=set(G.node[s]['offsets'].keys()).intersection(set(G.node[end]['offsets'].keys()))
            if len(x)!=0:
                G.add_edge(end,s,ofrom='+',oto='+',paths=x )
            else:
                print "!",source, sink, end, s, G.node[end]['offsets'], G.node[s]['offsets']
    
    return G

def realign_all(G,  minlength=20,
                    minn=2,
                    maxlen=10000000,
                    wscore=1,
                    wpen=1,
                    seedsize=None,
                    maxmums=None,
                    gcmodel="sumofpairs",
                    sa64=False,
                    minsize=None,
                    complex=False):

    complexbubbles=dict()
    source2sink=dict()
    sink2source=dict()
    
    if minsize==None:
        minsize=minlength

    #detect all complex bubbles
    for b in bubbles.bubbles(G):

        if complex:
            if b.issimple():
                continue

        if b.minsize<minsize:
            continue

        pair=(b.source,b.sink)
        bubblenodes=b.nodes
        ordD=b.ordD

        sourcesamples=set(G.node[pair[0]]['offsets'].keys())
        sinksamples=set(G.node[pair[1]]['offsets'].keys())

        if sourcesamples!=sinksamples:
            logging.warn("Invalid bubble between %s and %s"%(str(pair[0]),str(pair[1])))
            continue
        
        sucs=set(G.successors(pair[0]))
        pres=set(G.predecessors(pair[1]))
        sucs.discard(pair[1])
        pres.discard(pair[0])
        
        complexbubbles[pair]=bubblenodes
        source2sink[pair[0]]=pair[1]
        sink2source[pair[1]]=pair[0]
     
    #join complex bubbles that share a source/sink node
    converged=False
    while not converged:
        converged=True
        for source,sink in complexbubbles:
            if sink in source2sink:
                #sink is also a source, combine the two
                #print "sink is also source, merging",(source,sink),"with",(sink,source2sink[sink])
                complexbubbles[(source,source2sink[sink])]=complexbubbles[(source,sink)]+complexbubbles[(sink,source2sink[sink])]
                del complexbubbles[(sink, source2sink[sink])]
                del complexbubbles[(source,sink)]
                source2sink[source]=source2sink[sink]
                sink2source[source2sink[sink]]=source
                del sink2source[sink]
                del source2sink[sink]
                converged=False
                break
            elif source in sink2source:
                #source is also sink, combine
                #print "source is also sink, merging",(sink2source[source],source),"with",(source,sink)
                complexbubbles[(sink2source[source],sink)]=complexbubbles[(source,sink)]+complexbubbles[(sink2source[source],source)]
                del complexbubbles[(sink2source[source],source)]
                del complexbubbles[(source,sink)]
                source2sink[sink2source[source]]=sink
                sink2source[sink]=sink2source[source]
                del sink2source[source]
                del source2sink[source]
                converged=False
                break
    
    #filter out any bubbles that are a subset of another complex bubble
    distinctbubbles=dict()
    values=complexbubbles.values()
    for source,sink in complexbubbles:
        nodes=complexbubbles[(source,sink)]
        for bubble in values:
            if set(bubble).issuperset(set(nodes)) and not set(bubble)==set(nodes):
                #contained, dont add
                break
        else:
            distinctbubbles[(source,sink)]=nodes

    logging.info("Realigning a total of %d bubbles"%len(distinctbubbles))
    i=1
    for source,sink in distinctbubbles:
        logging.info("Realigning bubble between <%s> and <%s> of size (in nodes) %d."%(source,sink,len(distinctbubbles[(source,sink)])))
        G=realign_bubble(G,source,sink, maxmums=maxmums,
                                        seedsize=seedsize,
                                        gcmodel=gcmodel,
                                        minlength=minlength,
                                        minn=minn,
                                        wscore=wscore,
                                        maxlen=maxlen,
                                        wpen=wpen,
                                        sa64=sa64)
        
        i+=1
    
    return G
