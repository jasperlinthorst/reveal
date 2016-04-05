# -*- coding: utf-8 -*-
"""
Created on Fri Nov 13 13:24:31 2015

@author: jasperlinthorst
"""

import networkx as nx

def rc(seq):
    d = d = {'A':'T','C':'G','G':'C','T':'A','N':'N','a':'t','c':'g',\
        'g':'c','t':'a','n':'n','Y':'R','R':'Y','K':'M','M':'K',\
        'S':'S','W':'W','B':'V','V':'B','D':'H','H':'D','N':'N',\
        'X':'X','-':'-'}
    return "".join([d[b] for b in reversed(seq)])

class Caller:
    
    def __init__(self, gfafile, logger):
        self.gfafile = gfafile
        self.logging = logger
        self.graph = nx.DiGraph()
        self.parse()
    
    def call(self,minlength=1,maxlength=1000,allvar=False,inversions=False,indels=False,snps=False,multi=False):
        import seqal
        self.logging.info("%d nodes in graph.", self.graph.number_of_nodes())
    
        g=self.graph
        
        try:
            #for g in nx.connected_component_subgraphs(self.graph.to_undirected()):
            #self.logging.info("Iterating through subgraph of %d nodes.", g.number_of_nodes())
            #for node in nx.topological_sort(g):
            for node in self.graph.nodes():
                if len(g[node])>1: #split
                    name='N/A'
                    assert(len(g[node])<=self.graph.graph['samples'])
                    alleles=[]
                    allele_nodes=[]
                    allele_samples=[]
                    allelelengths=[]
                    indel=False
                    for nei in g[node]:
                        #if g.node[nei]['sample'].issubset(g.node[node]['sample']):
                        if not(g.node[node]['sample']==g.node[nei]['sample']):
                            s=g.node[nei]['seq']
                            allele_nodes.append(g.node[nei])
                            allele_samples.append(len(g.node[nei]['sample'].intersection(g.node[node]['sample'])))
                            alleles.append(s)
                            allelelengths.append(len(s))
                        else:
                            indel=True #direct neighbor has the same set of samples, then it has to be an indel
                        
                    #TODO: find join node and only(!) output simple bubbles
                    #if len(alleles)<=0:
                    #    print g.node[node]['sample'], [g.node[nei]['sample'] for nei in g[node]], g.node[node]['seq']
                    #    print g[node]
                    
                    assert(len(alleles)>=1)
                    
                    if allvar:
                        inversions,indels,snps,multi=True,True,True,True
                    
                    output=False
                    
                    if min(allelelengths)>=minlength and max(allelelengths)<=maxlength:
                    
                        if inversions and not(output):
                            if len(alleles)==2 and min(allelelengths)*1.5>max(allelelengths): #only test inversion if two alleles and not too variable in size
                                #print "NW-align:",allelelengths
                                score=seqal.nw_align(alleles[0],alleles[1])[2]
                                al1,al2,rcscore=seqal.nw_align(rc(alleles[0]),alleles[1])
                                #rcscore=1
                                #score=0
                                if rcscore>score:
                                    if (rcscore-score)>(sum(allelelengths)/2):
                                        #alleles[0]=al1
                                        #alleles[1]=al2
                                        output=True
                                        name='inversion'
                        
                        if indels and not(output):
                            if indel:
                                name='insert'
                                for vnode in allele_nodes:
                                    if vnode['coordsample']==g.node[node]['coordsample']:
                                        name='delete'
                                        break
                                    output=True
                        
                        if snps and not(output):
                            if min(allelelengths)==max(allelelengths)==1 and len(alleles)==2:
                                name='SNP'
                                output=True
                        
                        if multi and not(output):
                            if len(alleles)>2:
                                name='multi'
                                output=True
                    
                    if output or allvar:
                        print g.node[node]['start']+len(g.node[node]['seq']),g.node[node]['coordsample'],g.node[node]['coordcontig'],allelelengths,allele_samples,name,alleles
                        
        except IOError:
            pass
    
    def parse(self):
        f=open(self.gfafile,'r')
        self.graph.graph['samples']=[]
        sep=";"
        for line in f:
            if line.startswith('H'):
                h=line.split()
                tag=h[1].split(':')
                if tag[0]=="ORI":
                    for sample in tag[2].rstrip(sep).split(sep):
                        self.graph.graph['samples'].append(sample)
            
            if line.startswith('S'):
                s=line.strip().split('\t')
                nodeid=s[1]
                ann={}
                for v in s[4:]:
                    v=v.split(':')
                    if v[2].find(sep)!=-1: #multi-valued
                        ann[v[0]]=v[2].rstrip(sep).split(sep)
                    else:
                        ann[v[0]]=v[2]
                if "ORI" not in ann: #not a reveal graph, so no metadata on nodes
                    self.graph.add_node(nodeid,sample={gfafile},seq=s[2].upper(),contig={nodeid},aligned=0)
                else: #there are annotations on the nodes, so use them
                    if not(isinstance(ann['ORI'], list)):
                        ann['ORI']=[ann['ORI']]
                    if not(isinstance(ann['CTG'], list)):
                        ann['CTG']=[ann['CTG']]
                    ann['ORI']=set(ann['ORI'])
                    ann['CTG']=set(ann['CTG'])
                    self.graph.add_node(nodeid,sample=ann['ORI'],contig=ann['CTG'],coordsample=ann['CRD'],coordcontig=ann['CRDCTG'],start=int(ann['START']),seq=s[2].upper())
            
            #L    206    +    155    +    0M
            if line.startswith('L'):
                e=line.strip().split()
                assert(not self.graph.has_edge(e[1],e[3]))
                assert(not self.graph.has_edge(e[3],e[1]))
                self.graph.add_edge(e[1],e[3],ofrom=e[2],oto=e[4],cigar=e[5])

