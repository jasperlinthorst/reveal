import networkx as nx
from utils import *
import sys
import bubbles

def stats_cmd(args):
    if len(args.gfa)!=1:
        logging.fatal("Specify 1 gfa file.")
        return
    
    if not args.gfa[0].endswith('.gfa'):
        logging.fatal("Use .gfa as extension of the gfa file.")
        return
    
    stats(args.gfa[0])

def stats(gfafile):
    G=nx.DiGraph()
    read_gfa(gfafile,None,"",G)
    nsamples=len(G.graph['samples'])

    stats=dict()
    seqperngenomes=dict()
    
    i=1
    for sample in G.graph['samples']:
        seqperngenomes[i]=0
        i+=1
    
    for node,data in G.nodes(data=True):
        seqperngenomes[len(data['offsets'])]+=len(data['seq'])

    for n in seqperngenomes:
        stats["Sequence observed in %d genomes"%n]=seqperngenomes[n]
    
    #count bubble stats
    complexbubbles=0
    simplebubbles=0
    for bubble in bubbles.bubbles(G):
        pair,bubblenodes,size,ordD=bubble
        if bubbles.issimple(G,pair[0],pair[1]):
            simplebubbles+=1
        else:
            complexbubbles+=1
    
    #chain stats
    chain=[]
    chainweight=0
    chainpenalty=0
    chainlength=0
    chainlengthbp=0
    for node,data in G.nodes(data=True):
        offsets=data['offsets']
        l=len(data['seq'])
        if len(offsets)==nsamples:
            coords=tuple([offsets[k] for k in sorted(offsets.keys())])
            chain.append((coords,l))
            chainweight+=l*len(offsets)
            chainlengthbp+=l
            chainlength+=1
    
    chain.sort(key=lambda l: l[0])
    
    ppoint=tuple([c+chain[0][1] for c in chain[0][0]])
    for point,length in chain[1:]:
        for i in range(len(point)):
            assert(point[i]>=ppoint[i])
        p=sumofpairs(ppoint,point)
        ppoint=tuple([c+length for c in point])
        chainpenalty+=p
    
    stats["Chain length"]=chainlength
    stats["Chain length basepairs"]=chainlengthbp
    stats["Chain weight"]=chainweight
    stats["Chain penalty"]=chainpenalty
    stats["Chain score"]=chainweight-chainpenalty
    stats["Number of bubbles (total)"]=complexbubbles+simplebubbles
    stats["Number of bubbles (simple)"]=simplebubbles
    stats["Number of bubbles (complex)"]=complexbubbles
    stats["Number of samples"]=len(G.graph['samples'])
    #stats["Samples"]=",".join(G.graph['samples'])
    stats["Number of nodes"]=G.number_of_nodes()
    stats["Number of edges"]=G.number_of_edges()
    
    for label in sorted(stats.keys()):
        sys.stdout.write("%s:\t%s\n"%(label.ljust(35),str(stats[label]).rjust(15)))
