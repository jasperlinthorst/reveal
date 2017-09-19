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

    if args.exp==None:
        args.exp=len(G.graph['samples'])
    
    if args.all:
        G=realign_all(G,minscore=args.minscore,minlength=args.minlength,minn=args.minn,exp=args.exp,wscore=args.wscore,wpen=args.wpen,maxsize=args.maxsize,maxlen=args.maxlen,sa64=args.sa64)
    else:
        G=realign_bubble(G,args.source,args.sink,minscore=args.minscore,minlength=args.minlength,minn=args.minn,exp=args.exp,wscore=args.wscore,wpen=args.wpen,maxsize=args.maxsize,maxlen=args.maxlen,sa64=args.sa64,pcutoff=args.pcutoff)
    
    if args.outfile==None:
        fn=args.graph[0].replace(".gfa",".realigned.gfa")
    else:
        fn=args.outfile
    
    write_gfa(G,"",outputfile=fn)

def realign_bubble(G,source,sink,minscore=0,minlength=20,minn=2,maxsize=100,maxlen=10000000,exp=2,wscore=3,wpen=1,pcutoff=None,sa64=False):
    #print "Realigning graph between %s and %s"%(source,sink)
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
            bubblenodes.append(node)
        if node==sink:
            add=False
    
    sg=G.subgraph(bubblenodes)
    d={}
    aobjs=[]
    cumsum=0
    
    #extract all paths
    for sample in sourcesamples.intersection(sinksamples): #should be equal..
        #TODO: can also be that input was an assembly graph! What to do...
        for i,seq in enumerate(extract(sg,sample)): #has to be just one component
            pass
        assert(i==0)
        cumsum+=len(seq)
        if len(seq)>0:
            aobjs.append((sample,seq))
        
        if cumsum>maxlen:
            print "Bubble (%s,%s) is too big. Increase --maxlen."%(source,sink)
            return G
    
    schemes.wpen=wpen
    schemes.wscore=wscore
    schemes.exp=exp
    schemes.pcutoff=pcutoff
    ng,idx=align(aobjs,global_align=False,minscore=minscore,minlength=minlength,minn=minn,sa64=sa64)
    T=idx.T
    prune_nodes(ng,T)
    
    seq2node(ng,T)
    mapping={}
    startnodes=set()
    endnodes=set()
    
    for node,data in ng.nodes(data=True): 
        if len(ng.predecessors(node))==0:
            startnodes.add(nn)
        if len(ng.successors(node))==0:
            endnodes.add(nn)
        corrected=dict()
        for sample in data['offsets']:
            corrected[sample]=data['offsets'][sample]+G.node[source]['offsets'][sample]
        ng.node[node]['offsets']=corrected
        mapping[node]=nn
        nn+=1
    
    ng=nx.relabel_nodes(ng,mapping) #relabel nodes to get rid of Intervals
    
    #store edges that have to be reconnected after removing
    pre=G.predecessors(source)
    suc=G.successors(sink)
    
    for node in bubblenodes: #remove all bubblenodes from the original graph except for source and sink
        G.remove_node(node)
    for node,data in ng.nodes(data=True): #add nodes from newly aligned graph to original graph, except for source and sink, update those with seq attribute
        G.add_node(node,data)
    for edge in ng.edges():
        G.add_edge(*edge)
    
    #reconnect nodes
    for p in pre:
        for start in startnodes:
            G.add_edge(p,start)
    for s in suc:
        for end in endnodes:
            G.add_edge(end,s)
    
    return G

def realign_all(G,minscore=0,minlength=20,minn=2,maxlen=10000000,exp=2,wscore=3,wpen=1,maxsize=100,sa64=False):
    complexbubbles=dict()
    source2sink=dict()
    sink2source=dict()
    
    #detect all complex bubbles
    for pair,bubblenodes,size,ordD in bubbles.bubbles(G):
        sourcesamples=set(G.node[pair[0]]['offsets'].keys())
        sinksamples=set(G.node[pair[1]]['offsets'].keys())

        if not isinstance(sourcesamples,set):
            print pair,"skipping, not set"
            continue
        if sourcesamples!=sinksamples:
            print pair,sourcesample,sinksamples
            continue
        
        sucs=set(G.successors(pair[0]))
        pres=set(G.predecessors(pair[1]))
        sucs.discard(pair[1])
        pres.discard(pair[0])
        
        simple=True
        if len(sucs)!=len(pres) or len(sucs)<size or len(pres)<size:
            simple=False
        else:
            for suc in sucs:
                if len(G.successors(suc))!=1:
                    simple=False
            for pre in pres:
                if len(G.predecessors(pre))!=1:
                    simple=False

        if not simple and size<maxsize:
            assert(len(bubblenodes)<maxsize)
            complexbubbles[pair]=bubblenodes
            source2sink[pair[0]]=pair[1]
            sink2source[pair[1]]=pair[0]
     
    #join complex bubbles that share a source/sink nodes
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
                #print source,sink,"is contained"
                break
        else:
            distinctbubbles[(source,sink)]=nodes
    
    for bubble in distinctbubbles:
        assert(len(distinctbubbles[bubble])<500)

    #print "Total of %d bubbles to realign"%len(distinctbubbles)
    i=1
    for source,sink in distinctbubbles:
        #print i,"realigning",source,sink,len(distinctbubbles[(source,sink)])
        G=realign_bubble(G,source,sink,minscore=minscore,minlength=minlength,minn=minn,exp=exp,wscore=wscore,wpen=wpen)
        i+=1
    
    return G

