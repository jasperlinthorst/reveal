
from utils import *
import sys


def bubbles__(G):
    stack={}
    bubblenodestack={}
    for i,node in enumerate(nx.topological_sort(G)):
        
        for (source,source_idx) in stack.keys():
            bubblenodestack[(source,source_idx)].append(node)

        nei=set(G[node].keys())
        
        if len(nei)>1: #potential source
            stack[(node,i)]=nei
            bubblenodestack[(node,i)]=[node]

        ine=[v for v,t in G.in_edges(node)]
        if len(ine)>1: #potential sink
            for (source,source_idx) in stack.keys():
                
                stack[(source,source_idx)].discard(node)
                for v in ine:
                    stack[(source,source_idx)].discard(v)
                
                # for (source,source_idx) in stack.keys():
                if stack[(source,source_idx)]==set():
                    yield Bubble(G,source,node,source_idx,i,bubblenodestack[(source,source_idx)])
                    del stack[(source,source_idx)]
                    del bubblenodestack[(source,source_idx)]
        else:
            for (source,source_idx) in stack.keys():
                for v in ine:
                    stack[(source,source_idx)].discard(v)

def bubbles_(G):
    outstack={}
    instack={}
    bubblenodestack={}
    for i,node in enumerate(nx.topological_sort(G)):
        
        for (source,source_idx) in outstack.keys():
            bubblenodestack[(source,source_idx)].append(node)

        incoming=set([v for v,t in G.in_edges(node)])
        outgoing=set([t for v,t in G.out_edges(node)])

        if len(outgoing)>1: #potential source, open a bubble
            outstack[(node,i)]=outgoing
            bubblenodestack[(node,i)]=[node]

        for source,source_idx in outstack.keys():

            outstack[(source,source_idx)].discard(node)

            if outstack[(source,source_idx)]==set():
                # print "bubble: %s <-> %s"%(source,node)
                yield Bubble(G,source,node,source_idx,i,bubblenodestack[(source,source_idx)])

                del outstack[(source,source_idx)]
                del bubblenodestack[(source,source_idx)]
            else:
                for vo in outgoing:
                    outstack[(source,source_idx)].add(vo)

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
        entrance2candidateidx=dict()
        sspairs=[]
        #prevEnt=None
        prevEnti=None
        alternativeEntrance={}
        previousEntrance={}

        ordD={}

        # assert(type(G)==nx.DiGraph)

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
                entrance2candidateidx[(v,0)]=len(candidates)-1
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
                reportsuperbubble(candidates[0],candidates[-1],candidates,previousEntrance,alternativeEntrance,G,ordD,ordD_,outchild,outparent,sspairs,entrance2candidateidx)

        return ordD,ordD_,sspairs
    
    def reportsuperbubble(vstart,vexit,candidates,previousEntrance,alternativeEntrance,G,ordD,ordD_,outchild,outparent,sspairs,entrance2candidateidx):
        
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
                    
                    # ne=None
                    # for candidate in candidates[entrance2candidateidx[(s,0)]+1:]:
                    #     if candidate[1]==0:
                    #         ne=candidate
                    ne=nextentrance(candidates,s)

                    if ne!=None:
                        reportsuperbubble(ne, candidates[-1], candidates, previousEntrance, alternativeEntrance, G, ordD, ordD_, outchild, outparent,sspairs,entrance2candidateidx)
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

    allpairs=[]
    for pair in sspairs:
        allpairs.append((pair[0],pair[1]))

    allpairs.sort(key=lambda a: ordD[a[0]])#,reverse=True) #sort by topological order of the source

    for v,u in allpairs:
        bubblenodes=ordD_[ordD[v]:ordD[u]+1]
        sourcenode=G.node[v]
        sourcesamples=set(sourcenode['offsets'].keys())
        sinknode=G.node[u]
        sinksamples=set(sinknode['offsets'].keys())

        if sinksamples!=sourcesamples:
            logging.debug("Invalid bubble, between %s and %s"%(v,u))
            continue

        if len(bubblenodes)==2: #only source sink, no variation
            continue

        yield Bubble(G,v,u,ordD[v],ordD[u],bubblenodes)

def bubbles_cmd(args):
    if len(args.graph)<1:
        logging.fatal("Specify a gfa file to extract bubbles.")
        return
    
    G=nx.DiGraph()
    read_gfa(args.graph[0],None,"",G,remap=False)

    # bubbles(G)
    # sys.exit(0)

    sys.stdout.write("#source\tsink\tsubgraph\ttype\n")
    for i,g in enumerate(nx.weakly_connected_component_subgraphs(G)):
        logging.info("Reporting bubbles for subgraph: %d"%i)
        allcomplexnodes=[]
        for b in bubbles(g):
            if type(b)!=tuple:
                t=b.issimple()
                sys.stdout.write("%s\t%s\t%s\t%s\n"%(b.source if type(b.source)!=str else '<start>',
                                                    b.sink if type(b.sink)!=str else '<end>',
                                                    ",".join([str(x) for x in b.nodes if type(x)!=str]),
                                                    'simple' if t else 'complex'))

                if not t:
                    if args.exportcomplex:
                        if args.separate:
                            sg=g.subgraph(set(b.nodes))
                            if args.gml:
                                write_gml(sg,None,outputfile=args.graph[0].replace(".gfa",".%d.%d.complex.gml"%(b.source,b.sink)),partition=False)
                            else:
                                write_gfa(sg,None,remap=False,outputfile=args.graph[0].replace(".gfa","%d.%d.complex.gfa"%(b.source,b.sink)))
                        else:
                            allcomplexnodes+=b.nodes

        if args.exportcomplex and not args.separate:
            sg=g.subgraph(allcomplexnodes)
            if args.gml:
                write_gml(sg,None,outputfile=args.graph[0].replace(".gfa",".complex.gml"),partition=False)
            else:
                write_gfa(sg,None,remap=False,outputfile=args.graph[0].replace(".gfa",".complex.gfa"))

def rearrangements_cmd(args):

    G=nx.MultiDiGraph() #if we parse a DiGraph, the edges introduced by structural variants will be ignored
    
    logging.debug("Reading graph...")
    read_gfa(args.graph[0],None,"",G)
    logging.debug("Done.")

    logging.info("Determine rearrangement edges...")
    if type(G)==nx.MultiDiGraph: #convert to DiGraph first so we can actually toposort it
        rearrangements=MultiGraphToDiGraph(G)
    logging.info("Done (%d)."%len(rearrangements))

    gori=sorted([p for p in G.graph['paths'] if not p.startswith('*')])

    if args.reference==None:
        args.reference=gori[0]

    cds=G.graph['path2id'][args.reference] if args.reference in G.graph['path2id'] else G.graph['path2id'][gori[0]]

    sys.stdout.write("#reference\tapproximate_pos\tcontigs\tsource\tsink\tinvert\tpaths\n")

    for b in rearrangements:
        v,u,k,d=b

        if type(v)==str or type(u)==str:
            continue #just start/end
        else:
            
            paths=[G.graph['id2path'][sid] for sid in d['paths']] #all paths that go through the rearrangement edge

            for p in G.node[u]['offsets'].keys():
                if G.graph['id2path'][p].startswith(args.reference):
                    vcds=p
                    break
            else:
                logging.warn("Edge %s could not be located on reference: %s."% (str((v,u)), args.reference))
                vcds=G.node[u]['offsets'].keys()[0]

            vpos=G.node[u]['offsets'][vcds]



            contigs=[]
            for p in d['paths']:
                path=G.graph['id2path'][p]
                if path.startswith("*"):
                    contigs.append(path)

            sys.stdout.write("%s\t"*7 % (G.graph['id2path'][vcds], vpos, contigs, v, u, d['oto']==d['ofrom'], ",".join(paths)))
            
            sys.stdout.write("\n")
            sys.stdout.flush()

    logging.info("Done")

def variants_cmd(args):
    if len(args.graph)<1:
        logging.fatal("Specify a gfa file to extract bubbles.")
        return
    
    reference=args.reference
    g=nx.DiGraph() #if we parse a DiGraph, the edges introduced by structural variants will be ignored
    
    logging.debug("Reading graph...")
    read_gfa(args.graph[0],None,"",g)
    logging.debug("Done.")

    complexbubblenodes=[]
    
    if 'paths' in g.graph:
        gori=sorted([p for p in g.graph['paths'] if not p.startswith('*')])
    else:
        gori=[]

    if args.reference==None:
        args.reference=gori[0]
        logging.warn("No reference specified as a coordinate system, use %s where possible."%args.reference)
        args.reference=g.graph['path2id'][args.reference]
    else:
        if args.reference in g.graph['path2id']:
            args.reference=g.graph['path2id'][args.reference]
        else:
            logging.fatal("Specified reference (%s) not available in graph, graph knows of: %s."%(args.reference,str(g.graph['paths'])))
            sys.exit(1)
    
    try:
        if not args.fastaout and not args.bedout and not args.vcfout:
            sys.stdout.write("#reference\tpos_start\tpos_end\tsource_size\tsink_size\tmax_allele_size\tmin_allele_size\tdiff_allele_size\tsource\tsink\tsource_seq\tsink_seq\ttype\tgenotypes")
            for sample in gori:
                sys.stdout.write("\t%s"%sample)
            sys.stdout.write("\n")
        elif args.vcfout:
            sys.stdout.write("##fileformat=VCFv4.0\n")#?
            sys.stdout.write("##source=REVEAL\n")
            for sid in g.graph['id2path']:
                size=g.graph['id2end'][sid]
                sys.stdout.write("##contig=<ID=%s,LENGTH=%d>\n"%(g.graph['id2path'][sid],size))
            sys.stdout.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")
            sys.stdout.write("##FORMAT=<ID=SZ,Number=1,Type=Integer,Description=\"Allele size\">\n")
            sys.stdout.write("##INFO=<ID=reveal_diffsize,Number=1,Type=Integer,Description=\"Difference between the shortest and longest allele.\">\n")
            sys.stdout.write("##INFO=<ID=reveal_source,Number=1,Type=String,Description=\"Source of the node pair.\">\n")
            sys.stdout.write("##INFO=<ID=reveal_sink,Number=1,Type=String,Description=\"Sink of the node pair.\">\n")
            sys.stdout.write("##INFO=<ID=reveal_bubbletype,Number=1,Type=String,Description=\"Simplistic interpretation of the variant.\">\n")
            sys.stdout.write("##INFO=<ID=reveal_start,Number=1,Type=String,Description=\"Start position on the specified reference.\">\n")
            sys.stdout.write("##INFO=<ID=reveal_end,Number=1,Type=String,Description=\"End position on the specified reference.\">\n")
            sys.stdout.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT")
            for sample in gori:
                sys.stdout.write("\t%s"%sample)
            sys.stdout.write("\n")

        for bi,b in enumerate(bubbles(g)):

            v=Variant(b)
            
            if v.maxsize<args.minsize:
                continue

            if v.maxsize-v.minsize<args.mindiff:
                continue

            if args.maxdiff!=None and v.maxsize-v.minsize>args.maxdiff:
                continue

            if v.vtype!=args.type and args.type!='all':
                continue
            
            genotypestr=",".join(v.genotypes)
            
            if args.nogaps:
                if v.spans_gap:
                    continue

            minflank=min([len(g.node[v.source]['seq']),len(g.node[v.sink]['seq'])])
            
            if minflank<args.minflank:
                continue

            if args.reference in v.vpos:
                cds=args.reference
            else: #source does not occur on specified reference, pick any other path that does have a location for this variant
                if args.refonly: #skip the variant if its not positionable on the reference
                    continue
                for cds in v.vpos.keys():
                    if not g.graph['id2path'][cds].startswith('*'): #use ref layout if its there
                        break

            sourcelen=len(g.node[v.source]['seq'])
            sinklen=len(g.node[v.sink]['seq'])
            
            startpos=g.node[v.source]['offsets'][cds]+sourcelen
            endpos=g.node[v.sink]['offsets'][cds]

            if args.fastaout:
                if args.split:
                    with open("%s_%s.fasta"%(v.source,v.sink),'w') as of:
                        for i,seq in enumerate(v.genotypes):
                            if seq!='-':
                                of.write(">%s:%d-%d_%d\n"%(g.graph['id2path'][cds],startpos,endpos,i))
                                of.write("%s\n"%seq)
                else:
                    for i,seq in enumerate(v.genotypes):
                        if seq!='-':
                            sys.stdout.write(">%s:%d-%d_%d\n"%(g.graph['id2path'][cds],startpos,endpos,i))
                            sys.stdout.write("%s\n"%seq)
                continue

            if args.bedout:
                sys.stdout.write("%s\t%d\t%s\t%s\n"%(g.graph['id2path'][cds],startpos,endpos,v.vtype))
                continue

            allelesizes=[]

            for gt in v.genotypes:
                if gt=='-':
                    allelesizes.append(0)
                else:
                    allelesizes.append(len(gt))
            
            maxa=max(allelesizes)
            mina=min(allelesizes)

            if args.vcfout:
                startpos+=1
                if maxa-mina>0:
                    startpos-=1
                    genotypes=[]
                    for gt in v.genotypes:
                        if gt=='-':
                            gt=""
                        genotypes.append(g.node[v.source]['seq'][-1:]+gt)
                    v.genotypes=genotypes

                if v.calls[g.graph['id2path'][cds]]!=0: #for vcf output flip alleles to make reference allele 0
                    v.genotypes[0],v.genotypes[v.calls[g.graph['id2path'][cds]]]=v.genotypes[v.calls[g.graph['id2path'][cds]]],v.genotypes[0]
                _calls=dict()
                for sample in v.calls:
                    if v.calls[sample]==v.calls[g.graph['id2path'][cds]]: #same allele as ref, so make 0
                        _calls[sample]=0
                    elif v.calls[sample]==0:
                        _calls[sample]=v.calls[g.graph['id2path'][cds]]
                    else:
                        _calls[sample]=v.calls[sample]
                v.calls=_calls

                sys.stdout.write("%s\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s"% (g.graph['id2path'][cds],
                                                                    startpos,
                                                                    ".",
                                                                    v.genotypes[0],
                                                                    ",".join(v.genotypes[1:]),
                                                                    ".",
                                                                    "PASS",
                                                                    "reveal_diffsize=%s;reveal_source=%s;reveal_sink=%s;reveal_bubbletype=%s;reveal_start=%d;reveal_end=%d"%(maxa-mina, 
                                                                                                            v.source if type(v.source)!=str else '<start>', 
                                                                                                            v.sink if type(v.sink)!=str else '<end>',
                                                                                                            v.vtype,
                                                                                                            startpos,
                                                                                                            endpos),
                                                                    "GT:SZ"
                                                                    ))

                for sample in gori:
                    if sample in v.calls:
                        sys.stdout.write("\t%s:%d"%(v.calls[sample], len(v.genotypes[v.calls[sample]])))
                    else:
                        sys.stdout.write("\t.")

            else:
                sys.stdout.write("%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%s\t%s\t%s\t%s\t%s\t%s"% (g.graph['id2path'][cds],
                                                                    startpos,
                                                                    endpos,
                                                                    sourcelen,
                                                                    sinklen,
                                                                    maxa,
                                                                    mina,
                                                                    maxa-mina,
                                                                    v.source if type(v.source)!=str else '<start>',
                                                                    v.sink if type(v.sink)!=str else '<end>',
                                                                    g.node[v.source]['seq'][-20:] if v.source in g else '-',
                                                                    g.node[v.sink]['seq'][:20] if v.sink in g else '-',
                                                                    v.vtype,
                                                                    genotypestr))
                for sample in gori:
                    if sample in v.calls:
                        sys.stdout.write("\t%s"%v.calls[sample])
                    else:
                        sys.stdout.write("\t-")
            
            sys.stdout.write("\n")

            sys.stdout.flush()
    except IOError:
        pass

class InvalidBubble(Exception):
    pass

class Bubble:
    def __init__(self,G,source,sink,source_idx,sink_idx,nodes):
        self.source=source
        self.sink=sink
        self.source_idx=source_idx
        self.sink_idx=sink_idx
        self.G=G
        self.nodes=nodes
        self.ordD={node:i for i,node in enumerate(nodes)}

        if len(self.nodes)<=2:
            raise InvalidBubble("Not a valid source sink pair as bubble")

        self.simple=None
        
        self.paths=set([k for k in G.node[self.source]['offsets'].keys() if not G.graph['id2path'][k].startswith("*")]) & set([k for k in G.node[self.sink]['offsets'].keys() if not G.graph['id2path'][k].startswith("*")])

        if 'seq' in G.node[self.source]:
            l=len(G.node[self.source]['seq'])
        else:
            l=(self.source[1]-self.source[0])

        self.allelesizes=[G.node[self.sink]['offsets'][p]-(G.node[self.source]['offsets'][p]+l) for p in self.paths]

        self.minsize=min(self.allelesizes)
        # self.minsize=min([len(G.node[node]['seq']) for node in self.nodes[1:-1]])

        assert(self.minsize>=0)

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

class Variant(Bubble):
    def __init__(self,bubble):

        Bubble.__init__(self,bubble.G,bubble.source,bubble.sink,bubble.source_idx,bubble.sink_idx,bubble.nodes)
        
        self.genotypes=[] #list of variant sequence
        self.vtype='undefined' #type definition of the variant
        self.calls=dict() #key is sample, value is index within genotypes
        self.vpos=dict() #key is sample, value is position within sample
        self.spans_gap=False

        gt=list(set(self.G.successors(self.source)) & set(self.nodes))
        gt.sort(key=lambda l: self.ordD[l])
        bsamples=set(self.G.node[self.source]['offsets'].keys())&set(self.G.node[self.sink]['offsets'].keys())

        # bsamplestmp=bsamples.copy()
        # if self.issimple():
        #     for i,v in enumerate(gt):
        #         if v==self.sink:
        #             self.genotypes.append('-')
        #         else:
        #             s=self.G.node[v]['seq']
        #             self.genotypes.append(s)

        #         for sampleid in self.G.node[v]['offsets'].keys():
        #             if sampleid in bsamplestmp:
        #                 self.calls[bubble.G.graph['id2path'][sampleid]]=i
        #                 bsamplestmp.discard(sampleid)
        # else:
        
        self.vtype="complex"

        seqd=dict()
        for sid in bsamples:
            seq=""
            for v in self.nodes[1:-1]: #determine sequence through the complex bubble; use the entire path as genotype
                if sid in self.G.node[v]['offsets']:
                    seq+=self.G.node[v]['seq']

            if seq=="":
                seq="-"

            if seq in seqd:
                seqd[seq].append(sid)
            else:
                seqd[seq]=[sid]

        self.genotypes=list(seqd.keys())
        for i,k in enumerate(self.genotypes):
            for sid in seqd[k]:
                self.calls[bubble.G.graph['id2path'][sid]]=i
        
        if self.issimple():
            if self.G.has_edge(self.source,self.sink):
                self.vtype='indel'
            elif len(self.genotypes)==2:
                if len(self.genotypes[0])==1 and len(self.genotypes[1])==1:
                    self.vtype='snv'
                else:
                    self.vtype='region'
            else:
                self.vtype='multi-allelic'

        for node in self.nodes:
            if 'N' in self.G.node[node]['seq']:
                self.spans_gap=True
                self.vtype="gap"
                break

        v=self.G.node[self.source]
        t=self.G.node[self.sink]
        o=set(v['offsets'].keys())&set(t['offsets'].keys())
        for s in o:
            self.vpos[s]=v['offsets'][s]+len(v['seq'])+1