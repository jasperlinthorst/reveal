
from utils import *
import sys

def bubbles(G):
    def entrance(G,v):
        for c in G.successors(v):
            if len(list(G.predecessors(c)))==1:
                return True
        return False 

    def exit(G,v):
        for p in G.predecessors(v):
            if len(list(G.successors(p)))==1:
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
        structural_variants=[]

        if type(G)==nx.MultiDiGraph: #convert to DiGraph first so we can actually toposort it
            logging.debug("Converting MultiDigraph to DiGraph, before bubble extraction.")
            orgpaths=set([G.graph['path2id'][p] for p in G.graph['paths'] if p.startswith('*')])
            refpaths=set([G.graph['path2id'][p] for p in G.graph['paths'] if not p.startswith('*')])
            refpathnames=[p for p in G.graph['paths'] if not p.startswith('*')]
            for e0,e1,k,d in G.edges(keys=True,data=True):
                if len(d['paths'] & refpaths)==0: #edge that exclusively represents a structural event 
                    # logging.debug("STRUCTURAL VARIANT: %s %s %s"%(e0,e1,d))
                    structural_variants.append((e0,e1,k,d))

            G.remove_edges_from(structural_variants)
            G.graph['paths']=refpathnames

        logging.debug("Topologically sort the graph.")
        ordD_=list(nx.topological_sort(G))
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
        
        return ordD,ordD_,sspairs,structural_variants
    
    def reportsuperbubble(vstart,vexit,candidates,previousEntrance,alternativeEntrance,G,ordD,ordD_,outchild,outparent,sspairs):
        
        if (vstart[0] == None) or (vexit[0] == None) or (ordD[vstart[0]] >= ordD[vexit[0]]):
            del candidates[-1]
            return

        si=previousEntrance[vexit[0]]
        
        if si==None:
            del candidates[-1]
            return

        s=ordD_[si]
        while ordD[s] >= ordD[vstart[0]]:
            valid = validatesuperbubble(s, vexit[0], ordD, ordD_, outchild, outparent, previousEntrance, G)
            
            if (valid==s):
                break

            if (valid==alternativeEntrance[s]):
                break

            if valid==-1:
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
    
    ordD,ordD_,sspairs,structural_variants=superbubble(G)
    
    allpairs=[]
    for v,u,k,d in structural_variants:
        if d['ofrom']==d['oto'] and (G.has_edge(v,u) or G.has_edge(u,v)):
            continue
        allpairs.append((v,u,d))

    for pair in sspairs:
        allpairs.append((pair[0],pair[1],None))

    allpairs.sort(key=lambda a: ordD[a[0]],reverse=True) #sort by topological order of the source

    for v,u,d in allpairs:
        if d==None:
            bubblenodes=ordD_[ordD[v]:ordD[u]+1]
            sourcenode=G.node[v]
            sourcesamples=set(sourcenode['offsets'].keys())
            sinknode=G.node[u]
            sinksamples=set(sinknode['offsets'].keys())
            if len(bubblenodes)==2: #only source sink, no variation
                continue
            yield Bubble(G,v,u,nodes=bubblenodes,ordD=ordD)
        else:
            yield v,u,d

def bubbles_cmd(args):
    if len(args.graph)<1:
        logging.fatal("Specify a gfa file to extract bubbles.")
        return
    
    G=nx.MultiDiGraph()
    read_gfa(args.graph[0],None,"",G)
    
    allcomplexnodes=[]
    sys.stdout.write("#source\tsink\tsubgraph\ttype\n")
    for b in bubbles(G):
        if type(b)!=tuple:
            t=b.issimple()
            sys.stdout.write("%s\t%s\t%s\t%s\n"%(b.source if type(b.source)!=str else '<start>',
                                                b.sink if type(b.sink)!=str else '<end>',
                                                ",".join([str(x) for x in b.nodes if type(x)!=str]),
                                                'simple' if t else 'complex'))

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
    G=nx.MultiDiGraph() #if we parse a DiGraph, the edges introduced by structural variants will be ignored
    
    logging.debug("Reading graph...")
    read_gfa(args.graph[0],None,"",G)
    logging.debug("Done.")
    complexbubblenodes=[]
    
    if 'paths' in G.graph:
        gori=sorted([p for p in G.graph['paths'] if not p.startswith('*')])
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
        sys.stdout.write("#reference\tpos\tminflanksize\tsource\tsink\tsource_seq\tsink_seq\ttype\tgenotypes")
        for sample in gori:
            sys.stdout.write("\t%s"%sample)
        sys.stdout.write("\n")
    
    for b in bubbles(G):
        if type(b)!=tuple:
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
                for cds in v.vpos.keys():
                    if not G.graph['id2path'][cds].startswith('*'): #use ref layout if its there
                        break
            
            minflank=min([len(G.node[v.source]['seq']),len(G.node[v.sink]['seq'])])
            if minflank<args.minflank:
                continue

            genotypestr=",".join(v.genotypes)

            if args.nogaps:
                if 'N' in genotypestr:
                    continue

            sys.stdout.write("%s\t%s\t%d\t%s\t%s\t%s\t%s\t%s\t%s"% (G.graph['id2path'][cds],v.vpos[cds],minflank,
                                                                    v.source if type(v.source)!=str else '<start>',
                                                                    v.sink if type(v.sink)!=str else '<end>',
                                                                    G.node[v.source]['seq'][-20:] if v.source in G else '-',
                                                                    G.node[v.sink]['seq'][:20] if v.sink in G else '-',
                                                                    v.vtype,
                                                                    genotypestr))
            for sample in gori:
                if sample in v.calls:
                    sys.stdout.write("\t%s"%v.calls[sample])
                else:
                    sys.stdout.write("\t-")
            
            sys.stdout.write("\n")

        else: #structural variant, handle output differently
            v,u,d=b
            if type(v)==str or type(u)==str:
                continue #just start/end
            else:
                if args.reference in G.node[v]['offsets']:
                    cds=args.reference
                else:
                    for p in G.node[v]['offsets']:
                        if not G.graph['id2path'][p].startswith("*"): #make sure we output something useful instead of a location on a contig
                            cds=p
                            break
                    for p in G.node[u]['offsets']:
                        if not G.graph['id2path'][p].startswith("*"): #make sure we output something useful instead of a location on a contig
                            cdsu=p
                            break

                minflank=min([len(G.node[v]['seq']),len(G.node[u]['seq'])])
                sys.stdout.write("%s\t%s\t%d\t%s\t%s\t%s\t%s\t%s\t%s"% (G.graph['id2path'][cds],G.node[v]['offsets'][cds]+len(G.node[v]['seq']),minflank,
                                                                    v,
                                                                    u,
                                                                    G.node[v]['seq'][-20:] if v in G else '-',
                                                                    G.node[u]['seq'][:20] if u in G else '-',
                                                                    'struct_inv' if d['ofrom']!=d['oto'] else 'struct',
                                                                    G.graph['id2path'][cds]+" <--> "+G.graph['id2path'][cdsu]))
                for sample in gori:
                    if G.graph['path2id'][sample] in d['paths']:
                        sys.stdout.write("\t1")
                    else:
                        sys.stdout.write("\t0")
                
                sys.stdout.write("\n")

class InvalidBubble(Exception):
    pass

class Bubble:
    def __init__(self,G,source,sink,nodes=None,ordD=None):
        self.source=source
        self.sink=sink
        self.G=G

        if ordD==None or nodes==None:
            ordD_=list(nx.topological_sort(G))

        if ordD!=None:
            self.ordD=ordD
        else:
            self.ordD={}
            for i,v in enumerate(ordD_):
                self.ordD[v]=i
        
        if nodes!=None:
            self.nodes=nodes
        else:
            if self.ordD[source]>=self.ordD[sink]:
                raise InvalidBubble("Sink comes before Source, flip source and sink to form valid bubble.")
            self.nodes=ordD_[self.ordD[source]:self.ordD[sink]+1]

        if len(self.nodes)<=2:
            raise InvalidBubble("Not a valid source sink pair as bubble")

        self.simple=None
        
        #TODO: use source/sink offsets to calculate this!!

        self.paths=set(G.node[self.source]['offsets'].keys()) & set(G.node[self.sink]['offsets'].keys())

        if 'seq' in G.node[self.source]:
            l=len(G.node[self.source]['seq'])
        else:
            l=(self.source[1]-self.source[0])

        self.allelesizes=[G.node[self.sink]['offsets'][p]-(G.node[self.source]['offsets'][p]+l) for p in self.paths]

        self.minsize=min(self.allelesizes)
        # self.minsize=min([len(G.node[node]['seq']) for node in self.nodes[1:-1]])

        # self.maxsize=max([len(G.node[node]['seq']) for node in self.nodes[1:-1]])
        self.maxsize=max(self.allelesizes)

        self.cumsize=sum(self.allelesizes)
        # self.cumsize=sum([len(G.node[node]['seq'])*len(G.node[node]['offsets']) for node in self.nodes[1:-1]])
    
    def issimple(self):
        if self.simple==None:
            
            sucs=set(self.G.successors(self.source))
            pres=set(self.G.predecessors(self.sink))
            
            sucs.discard(self.sink)
            pres.discard(self.source)
            
            for suc in sucs:
                if len(list(self.G.successors(suc)))!=1 or list(self.G.successors(suc))[0]!=self.sink:
                    self.simple=False
                    return self.simple
            
            for pre in pres:
                if len(list(self.G.predecessors(pre)))!=1 or list(self.G.predecessors(pre))[0]!=self.source:
                    self.simple=False
                    return self.simple
            
            self.simple=True
            
            return self.simple
        else:
            return self.simple

    #returns the amount the left and right margin for indel positioning
    def getwiggle(self,minwiggle=0):
        if self.issimple():

            if self.G.has_edge(self.source,self.sink):
                #how far can we move this bubble to the right?

                if 'seq' in self.G.node[self.sink]:
                    sink=self.G.node[self.sink]['seq']
                else:
                    sink=""

                if 'seq' in self.G.node[self.source]:
                    source=self.G.node[self.source]['seq']
                else:
                    source=""

                vs=[self.G.node[n]['seq']+sink for n in self.nodes[1:-1]]
                lvs=[len(s) for s in vs]+[len(sink)]
                i=0
                while i<min(lvs) and sink[i]==vs[0][i]:
                    for v in vs[1:]:
                        if not v[i]==sink[i]:
                            break
                    i+=1
                
                vs=[source+self.G.node[n]['seq'] for n in self.nodes[1:-1]]
                lvs=[len(s) for s in vs]+[len(source)]
                j=1
                while j<=min(lvs) and source[-j]==vs[0][-j]:
                    for v in vs[1:]:
                        if not v[-j]==source[-j]:
                            break
                    j+=1

                return (minwiggle+j-1,minwiggle+i) #tuple with margin on the left and margin on the right
        
        return (minwiggle,minwiggle)

        # return (0,0)

class Variant(Bubble):
    def __init__(self,bubble):

        Bubble.__init__(self,bubble.G,bubble.source,bubble.sink,bubble.nodes,bubble.ordD)
        
        self.genotypes=[] #list of variant sequence
        self.vtype='undefined' #type definition of the variant
        self.calls=dict() #key is sample, value is index within genotypes
        self.vpos=dict() #key is sample, value is position within sample
        self.size=1 #length of the largest allele
        
        gt=list(self.G.successors(self.source))
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
