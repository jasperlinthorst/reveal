import networkx as nx
from utils import *
import sys
from bubbles import bubbles,Bubble,Variant

def stats_cmd(args):
    if len(args.gfa)!=1:
        logging.fatal("Specify 1 gfa file.")
        return
    
    if not args.gfa[0].endswith('.gfa'):
        logging.fatal("Use .gfa as extension of the gfa file.")
        return
    
    stats(args.gfa[0])

def stats(gfafile):
    stats=dict()
    G=nx.DiGraph()
    read_gfa(gfafile,None,"",G)
    samples=G.graph['paths']
    nsamples=len(samples)
    
    stats["Graph"]=os.path.basename(gfafile)
    stats["Number of samples"]=nsamples
    for i,sample in enumerate(samples):
        stats["Sample %d"%i]=sample

    stats["Number of connected components"]=0
    stats["A count"]=0
    stats["C count"]=0
    stats["G count"]=0
    stats["T count"]=0
    stats["N count"]=0
    for node,data in G.nodes(data=True):
        stats["A count"]+=data['seq'].count('A')
        stats["C count"]+=data['seq'].count('C')
        stats["G count"]+=data['seq'].count('G')
        stats["T count"]+=data['seq'].count('T')
        stats["N count"]+=data['seq'].count('N')

    #for each connected component
    for sgi,sub in enumerate(nx.connected_components(G.to_undirected())):
        stats["Number of connected components"]+=1
        sg=G.subgraph(sub)
        
        #determine samples in subgraph
        nsgsamples=1
        sgsamples=set()
        for node,data in sg.nodes(data=True):
            if len(data['offsets'])>nsgsamples:
                nsgsamples=len(data['offsets'])
            for sid in data['offsets']:
                sgsamples.add(G.graph['id2path'][sid])
        
        stats["Composition of component %d"%sgi]=",".join(sgsamples)

        seqperngenomes=dict()
        
        i=1
        for sample in sg.graph['paths']:
            seqperngenomes[i]=0
            i+=1
        
        for node,data in sg.nodes(data=True):
            seqperngenomes[len(data['offsets'])]+=len(data['seq'])

        for n in seqperngenomes:
            stats["Sequence observed in %d genomes"%n]=seqperngenomes[n]
        
        #count bubble stats
        complexbubbles=0
        simplebubbles=0
        snpcount=0
        indelcount=0
        multicount=0
        regioncount=0
        unknowncount=0
        
        for bubble in bubbles(sg):
            if bubble.issimple():
                simplebubbles+=1
            else:
                complexbubbles+=1
            
            v=Variant(bubble)
            
            if v.vtype=='snp':
                snpcount+=1
            elif v.vtype=='indel':
                indelcount+=1
            elif v.vtype=='multi-allelic':
                multicount+=1
            elif v.vtype=='region':
                regioncount+=1
            else:
                unknowncount+=1
        
        stats["Number of bubbles in component %d (total)"%sgi]=complexbubbles+simplebubbles
        stats["Number of bubbles in component %d (simple)"%sgi]=simplebubbles
        stats["Number of bubbles in component %d (complex)"%sgi]=complexbubbles
        stats["Number of variants in component %d (snps)"%sgi]=snpcount
        stats["Number of variants in component %d (indels)"%sgi]=indelcount
        stats["Number of variants in component %d (multi-allelic)"%sgi]=multicount
        stats["Number of variants in component %d (complex)"%sgi]=unknowncount
        stats["Number of variants in component %d (regions)"%sgi]=regioncount

        #chain stats
        chain=[]
        chainweight=0
        chainpenalty=0
        chainlength=0
        chainlengthbp=0
        
        for node,data in sg.nodes(data=True):
            if type(node)==str: #skip start and end nodes
                continue
            offsets=data['offsets']
            l=len(data['seq'])
            if len(offsets)==nsgsamples:
                coords=tuple([offsets[k] for k in sorted(offsets.keys())])
                chain.append((coords,l))
                chainweight+=l*((len(offsets)*(len(offsets)-1))/2) #sumofpairs score!
                chainlengthbp+=l
                chainlength+=1
        
        chain.sort(key=lambda l: l[0])
        
        ppoint=tuple([c+chain[0][1] for c in chain[0][0]])
        for point,length in chain[1:]:
            for i in range(len(point)):
                assert(point[i]>=ppoint[i])
            p=gapcost(ppoint,point) #sumofpairs penalty!
            ppoint=tuple([c+length for c in point])
            chainpenalty+=p
        
        stats["Chain length in component %d"%sgi]=chainlength
        stats["Chain length basepairs in component %d"%sgi]=chainlengthbp
        stats["Chain weight (sum-of-pairs) in component %d"%sgi]=chainweight
        stats["Chain penalty (sum-of-pairs) in component %d"%sgi]=chainpenalty
        stats["Chain score in component %d"%sgi]=chainweight-chainpenalty

    remove=[]
    for node in G.nodes():
        if type(node)==str:
            remove.append(node)
    G.remove_nodes_from(remove)
    
    stats["Number of nodes"]=G.number_of_nodes()
    stats["Number of edges"]=G.number_of_edges()

    for label in sorted(stats.keys()):
        sys.stdout.write("%s:\t%s\n"%(label.ljust(50),str(stats[label]).rjust(50)))

