
from utils import *
import sys

def bubbles_cmd(args):
    if len(args.graph)<1:
        logging.fatal("Specify a gfa file to extract bubbles.")
        return
    
    reference=args.reference
    G=nx.DiGraph()
    read_gfa(args.graph[0],None,"",G)
    complexbubblenodes=[]
    
    if 'samples' in G.graph:
        gori=sorted(G.graph['samples'])
    else:
        gori=[]
    
    sys.stdout.write("#source\tsink\tsubgraph\tref\tpos\tvariant")
    for sample in gori:
        sys.stdout.write("\t%s"%sample)
    sys.stdout.write("\n")
    
    for pair,bubblenodes,size,ordD in bubbles(G):
        gt=G.successors(pair[0])
        gt.sort(key=lambda l: ordD[l]) #topological sort neighbors
        sourcesamples=set(G.node[pair[0]]['offsets'].keys())
        simple=issimple(G,pair[0],pair[1])

        #determine genotypes, associate genotype with numerical value
        if len(gt)<size or not simple:
            #complex bubble; not possible to output the actual alleles, instead we can make a call based on branching at source node and output the nodes that make up the complex bubble so we can export and vizualize the subgraph
            genotypes=['N']*len(gt)
        else:
            genotypes=[None]*len(gt)
            for i,node in enumerate(gt):
                genotypes[i]=G.node[node]['seq']
                if node==pair[1]:
                    genotypes[i]="-"
        
        #make call, associate sample with the numerical value of the genotype
        tmpsource=sourcesamples.copy()
        calls={}
        for i,node in enumerate(gt):
            n=G.node[node]
            for sample in n['offsets'].keys():
                if sample not in tmpsource:
                    continue
                calls[sample]=i
                tmpsource.discard(sample)
        
        d=G.node[pair[0]]
        
        if reference==None or reference not in d['offsets'].keys():
            ref=d['offsets'].keys()[0]
            pos=d['offsets'][ref]+len(d['seq'])+1
        else:
            ref=reference
            pos=d['offsets'][reference]+len(d['seq'])+1
        
        sys.stdout.write("%d\t%d\t%s\t%s\t%s\t%s"%(pair[0],pair[1],",".join([str(x) for x in bubblenodes]),ref,pos,",".join(genotypes)))
        
        if args.exportcomplex and not simple:
            if args.separate:
                sg=G.subgraph(bubblenodes)
                if args.gml:
                    write_gml(sg,None,outputfile="%s_%s.gml"%(pair[0],pair[1]),partition=False)
                else:
                    write_gfa(sg,None,remap=False,outputfile="%s_%s.gfa"%(pair[0],pair[1]))
            else:
                complexbubblenodes+=bubblenodes
        
        for sample in gori:
            if sample in calls:
                sys.stdout.write("\t%s"%calls[sample])
            else:
                sys.stdout.write("\t-")
        sys.stdout.write("\n")
    
    if args.exportcomplex and not args.separate:
        sg=G.subgraph(set(complexbubblenodes))
        if args.gml:
            write_gml(sg,None,outputfile=args.graph[0].replace(".gfa",".complex.gml"),partition=False)
        else:
            write_gfa(sg,None,remap=False,outputfile=args.graph[0].replace(".gfa",".complex.gfa"))

def issimple(G,source,sink):
    gt=G.successors(source)
    sucs=set(gt)
    pres=set(G.predecessors(sink))
    sucs.discard(sink)
    pres.discard(source)
    simple=True
    for suc in sucs:
        if len(G.successors(suc))!=1:
            simple=False
    for pre in pres:
        if len(G.predecessors(pre))!=1:
            simple=False
    return simple

def bubbles(G):
    def entrance(G,v):
        for c in G.successors(v):
            if len(G.predecessors(c))==1:
                return True
        return False 

    def exit(G,v):
        for p in G.predecessors(v):
            if len(G.successors(p))==1:
                return True
        return False
    
    def nextentrance(candidates,v):
        #TODO: rewrite this
        for candidate in candidates[candidates.index((v,0))+1:]:
            if candidate[1]==0:
                return candidate
    
    def superbubble(G):
        candidates=[]
        sspairs=[]
        #prevEnt=None
        prevEnti=None
        alternativeEntrance={}
        previousEntrance={}
        ordD={}
        ordD_=nx.topological_sort(G)
        
        #construct candidates array
        for i,v in enumerate(ordD_):
            ordD[v]=i
            alternativeEntrance[v]=None
            previousEntrance[v]=prevEnti
            if exit(G,v):
                candidates.append((v,1))
            if entrance(G,v):
                candidates.append((v,0))
                #prevEnt=v
                prevEnti=i
        
        #construct outparent
        outparent=[None]*(len(ordD))
        for i,c in enumerate(ordD):
            tmp=[]
            for p in G.predecessors(c):
                tmp.append(ordD[p])
            if len(tmp)>0:
                outparent[ordD[c]]=min(tmp)
        
        #construct outchild
        outchild=[None]*(len(ordD))
        for i,c in enumerate(ordD):
            tmp=[]
            for p in G.successors(c):
                tmp.append(ordD[p])
            if len(tmp)>0:
                outchild[ordD[c]]=max(tmp)
        
        #loop
        while len(candidates)!=0:
            if candidates[-1][1]==0:
                del candidates[-1]
            else:
                reportsuperbubble(candidates[0],candidates[-1],candidates,previousEntrance,alternativeEntrance,G,ordD,ordD_,outchild,outparent,sspairs)
        
        return ordD,ordD_,sspairs
    
    def reportsuperbubble(vstart,vexit,candidates,previousEntrance,alternativeEntrance,G,ordD,ordD_,outchild,outparent,sspairs):
        if (vstart[0] == None) or (vexit[0] == None) or (ordD[vstart[0]] >= ordD[vexit[0]]):
            del candidates[-1]
            return
        si=previousEntrance[vexit[0]]
        s=ordD_[si]
        while ordD[s] >= ordD[vstart[0]]:
            valid = validatesuperbubble(s, vexit[0], ordD, ordD_, outchild, outparent, previousEntrance, G)
            if valid==s or valid==alternativeEntrance[s] or valid==-1:
                break
            alternativeEntrance[s] = valid
            s = valid

        del candidates[-1]
        if (valid == s):
            sspairs.append((s, vexit[0]))
            while (candidates[-1][0] is not s):
                if candidates[-1][1]==1:
                    ne=nextentrance(candidates,s)
                    if ne!=None:
                        reportsuperbubble(ne, candidates[-1], candidates, previousEntrance, alternativeEntrance, G, ordD, ordD_, outchild, outparent,sspairs)
                    else:
                        del candidates[-1]
                else:
                    del candidates[-1]
    
    def validatesuperbubble(startVertex, endVertex, ordD, ordD_, outchild, outparent, previousEntrance, G):
        start=ordD[startVertex]
        end=ordD[endVertex]
        if start+1!=end:
            oc=max(outchild[start:end])
            op=min(outparent[start+1:end+1])
        else:
            oc=outchild[start]
            op=outparent[end]
        if oc!=end:
            return -1
        if op==start:
            return startVertex
        if op==None:
            return -1
        elif entrance(G, ordD_[op]):
            return ordD_[op]
        elif previousEntrance[ordD_[op]]==None: #
            return -1
        else:
            return ordD_[previousEntrance[ordD_[op]]]
        return startVertex
    
    ordD,ordD_,sspairs=superbubble(G)
    
    sg=set()
    for pair in sspairs:
        size=(ordD[pair[1]]-ordD[pair[0]])-1
        bubblenodes=ordD_[ordD[pair[0]]:ordD[pair[1]]+1]
        sourcenode=G.node[pair[0]]
        sourcesamples=set(sourcenode['offsets'].keys())
        sinknode=G.node[pair[1]]
        sinksamples=set(sinknode['offsets'].keys())
        
        #TODO: supbub algorithm detects invalid bubbles; filter them out here
        if sourcesamples!=sinksamples:
            continue
        yield pair,bubblenodes,size,ordD


