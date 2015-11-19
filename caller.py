# -*- coding: utf-8 -*-
"""
Created on Fri Nov 13 13:24:31 2015

@author: jasperlinthorst
"""

import networkx as nx
import seqal

def rc(seq):
    d = {'A':'T','T':'A','G':'C','C':'G','N':'N'}
    return "".join([d[b] for b in reversed(seq)])

class Caller:
    def __init__(self, gfafile):
        self.gfafile = gfafile
        self.vcffile = gfafile.replace(".gfa",".vcf")
        self.graph = nx.DiGraph()
        self.parse()
    
    def call(self,minlength=1,maxlength=1000,allvar=False,inversions=False,indels=False,snps=False,multi=False):
        g=self.graph
        for node in g.nodes():
            if len(g[node])>1: #split
                assert(len(g[node])<=self.graph.graph['samples'])
                alleles=[]
                allelelengths=[]
                indel=False
                for nei in g[node]:
                    if g.node[nei]['sample'].issubset(g.node[node]['sample']) and not(g.node[node]['sample']==g.node[nei]['sample']):
                        s=g.node[nei]['seq']
                        alleles.append(s)
                        allelelengths.append(len(s))
                #TODO: find join node and only output only simple bubbles
                assert(len(alleles)>=1)
                output=allvar
                if min(allelelengths)>=minlength and max(allelelengths)<=maxlength:
                    if inversions and not(output):
                        if len(alleles)==2: #only test inversion if two alleles
                            score=seqal.nw_align(alleles[0],alleles[1])[2]
                            al1,al2,rcscore=seqal.nw_align(rc(alleles[0]),alleles[1])
                            if rcscore>score:
                                print al1
                                print al2
                                output=True
                    if indels and not(output):
                        output=indel
                    if snps and not(output):
                        if min(allelelengths)==max(allelelengths)==1 and len(alleles)==2:
                            output=True
                    if multi and not(output):
                        if len(alleles)>2:
                            output=True
                    if output:
                        print alleles
    
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
            
            #L	206	+	155	+	0M
            if line.startswith('L'):
                e=line.strip().split()
                assert(not self.graph.has_edge(e[1],e[3]))
                assert(not self.graph.has_edge(e[3],e[1]))
                self.graph.add_edge(e[1],e[3],ofrom=e[2],oto=e[4],cigar=e[5])

