
from utils import *
import sys

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
        
        logging.debug("Topologically sort the graph.")
        ordD_=nx.topological_sort(G)
        logging.debug("Done.")
        
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
        
        #print "reporting supbub",vstart,vexit

        si=previousEntrance[vexit[0]]
        
        if si==None:
            del candidates[-1]
            return

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
        
        # if sourcesamples!=sinksamples:
            #logging.warn("Invalid bubble detected between node %s and node %s."%pair)
            # continue

        if len(bubblenodes)==2: #only source sink, no variation
            continue

        yield Bubble(G,pair[0],pair[1],nodes=bubblenodes,ordD=ordD)

def bubbles_cmd(args):
    if len(args.graph)<1:
        logging.fatal("Specify a gfa file to extract bubbles.")
        return
    
    G=nx.DiGraph()
    read_gfa(args.graph[0],None,"",G)
    
    allcomplexnodes=[]
    sys.stdout.write("#source\tsink\tsubgraph\ttype\n")
    for b in bubbles(G):
        t=b.issimple()
        sys.stdout.write("%d\t%d\t%s\t%s\n"%(b.source,b.sink,",".join([str(x) for x in b.nodes]), 'simple' if t else 'complex'))

        if not t:
            if args.exportcomplex:
                if args.separate:
                    sg=G.subgraph(set(b.nodes))
                    if args.gml:
                        write_gml(sg,None,outputfile=args.graph[0].replace(".gfa",".%d.%d.complex.gml"%(b.source,b.sink)),partition=False)
                    else:
                        write_gfa(sg,None,remap=False,outputfile=args.graph[0].replace(".gfa","%d.%d.complex.gfa"%(b.source,b.sink)))
                else:
                    allcomplexnodes+=b.nodes

    if args.exportcomplex and not args.separate:
        sg=G.subgraph(allcomplexnodes)
        if args.gml:
            write_gml(sg,None,outputfile=args.graph[0].replace(".gfa",".complex.gml"),partition=False)
        else:
            write_gfa(sg,None,remap=False,outputfile=args.graph[0].replace(".gfa",".complex.gfa"))

def variants_cmd(args):
    if len(args.graph)<1:
        logging.fatal("Specify a gfa file to extract bubbles.")
        return
    
    reference=args.reference
    G=nx.DiGraph() #if we parse a DiGraph, the edges introduced by structural variants will be ignored
    read_gfa(args.graph[0],None,"",G)
    complexbubblenodes=[]
    
    if 'paths' in G.graph:
        gori=sorted(G.graph['paths'])
    else:
        gori=[]
    
    if args.reference==None:
        args.reference=gori[0]
        logging.warn("No reference specified as a coordinate system, use %s where possible."%args.reference)
    else:
        if args.reference in G.graph['path2id']:
            args.reference=G.graph['path2id'][args.reference]
        else:
            logging.fatal("Specified reference (%s) not available in graph, graph knows of: %s."%(args.reference,str(G.graph['paths'])))
            sys.exit(1)
    
    if not args.fastaout:
        sys.stdout.write("#reference\tpos\tminflanksize\tsource\tsink\ttype\tgenotypes")
        for sample in gori:
            sys.stdout.write("\t%s"%sample)
        sys.stdout.write("\n")
    
    for b in bubbles(G):
        v=Variant(b)
        
        if v.size<args.minsize:
            continue

        if v.vtype!=args.type and args.type!='all':
            continue
        
        if args.fastaout:
            if args.split:
                with open("%s_%s.fasta"%(v.source,v.sink),'w') as of:
                    for i,seq in enumerate(v.genotypes):
                        if seq!='-':
                            of.write(">%s_%d\n"%(v.source,i))
                            of.write("%s\n"%seq)
            else:
                for i,seq in enumerate(v.genotypes):
                    if seq!='-':
                        sys.stdout.write(">%s_%d\n"%(v.source,i))
                        sys.stdout.write("%s\n"%seq)
            
            continue
        
        if args.reference in v.vpos:
            cds=args.reference
        else:
            cds=v.vpos.keys()[0]
        
        minflank=min([len(G.node[v.source]['seq']),len(G.node[v.sink]['seq'])])
        if minflank<args.minflank:
            continue

        sys.stdout.write("%s\t%d\t%d\t%s\t%s\t%s\t%s"%(G.graph['id2path'][cds],v.vpos[cds],minflank,v.source,v.sink, v.vtype, ",".join(v.genotypes) ))
        for sample in gori:
            if sample in v.calls:
                sys.stdout.write("\t%s"%v.calls[sample])
            else:
                sys.stdout.write("\t-")
        sys.stdout.write("\n")

class Bubble:
    def __init__(self,G,source,sink,nodes=None,ordD=None):
        self.source=source
        self.sink=sink
        self.G=G

        if ordD==None or nodes==None:
            ordD_=nx.topological_sort(G)

        if ordD!=None:
            self.ordD=ordD
        else:
            self.ordD={}
            for i,v in enumerate(ordD_):
                self.ordD[v]=i
        
        if nodes!=None:
            self.nodes=nodes
        else:
            self.nodes=ordD_[self.ordD[source]:self.ordD[sink]+1]

        if len(self.nodes)==2:
            raise InvalidBubble("Not a valid source sink pair as bubble")

        self.simple=None

        self.minsize=min([len(G.node[node]['seq']) for node in self.nodes[1:-1]])
        self.seqsize=sum([len(G.node[node]['seq'])*len(G.node[node]['offsets']) for node in self.nodes[1:-1]])
    
    def issimple(self):
        if self.simple==None:
            
            sucs=set(self.G.successors(self.source))
            pres=set(self.G.predecessors(self.sink))
            
            sucs.discard(self.sink)
            pres.discard(self.source)
            
            for suc in sucs:
                if len(self.G.successors(suc))!=1 or self.G.successors(suc)[0]!=self.sink:
                    self.simple=False
                    return self.simple
            
            for pre in pres:
                if len(self.G.predecessors(pre))!=1 or self.G.predecessors(pre)[0]!=self.source:
                    self.simple=False
                    return self.simple
            
            self.simple=True
            
            return self.simple
        else:
            return self.simple 

class Variant(Bubble):
    def __init__(self,bubble):

        Bubble.__init__(self,bubble.G,bubble.source,bubble.sink,bubble.nodes,bubble.ordD)
        
        self.genotypes=[] #list of variant sequence
        self.vtype='undefined' #type definition of the variant
        self.calls=dict() #key is sample, value is index within genotypes
        self.vpos=dict() #key is sample, value is position within sample
        self.size=1 #length of the largest allele
        
        gt=self.G.successors(self.source)
        gt.sort(key=lambda l: self.ordD[l])

        if self.issimple():
            for v in gt:
                if v==self.sink:
                    self.genotypes.append('-')
                else:
                    s=self.G.node[v]['seq']
                    self.genotypes.append(s)
                    if len(s)>self.size:
                        self.size=len(s)
            #self.genotypes.sort()
        else:
            for v in gt:
                self.genotypes.append('N')
        
        if self.issimple():
            if self.G.has_edge(self.source,self.sink):
                self.vtype='indel'
            elif len(self.genotypes)==2:
                if len(self.genotypes[0])==1 and len(self.genotypes[1])==1:
                    self.vtype='snp'
                else:
                    self.vtype='region'
            else:
                self.vtype='multi-allelic'
        
        n=self.G.node[self.source]
        
        for s in n['offsets']:
            self.vpos[s]=n['offsets'][s]+len(n['seq'])+1

        #make call, associate sample with the numerical value of the genotype
        sourcesamples=set(self.G.node[self.source]['offsets'].keys())
        tmpsource=sourcesamples.copy()
        
        for i,node in enumerate(gt):
            n=self.G.node[node]
            #for sample in n['offsets'].keys():
            for sampleid in n['offsets'].keys():
                #if sample not in tmpsource:
                if sampleid not in tmpsource:
                    continue
                #self.calls[sample]=i
                self.calls[bubble.G.graph['id2path'][sampleid]]=i
                #tmpsource.discard(sample)
                tmpsource.discard(sampleid)

    def isinversion(self):
        if self.vtype=='undefined':
            return 'maybe'
            pass
        else:
            return False
    
