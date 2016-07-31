import reveal
import networkx as nx
import argparse

def main():
    desc="""
    Very simple (and memory naive) script to transform FALCON's p_ctg/a_ctg_all/a_ctg_base structure to a GFA kind of graph.
    """
    
    parser = argparse.ArgumentParser(prog="falcon2gfa", usage="falcon2gfa.py -h for usage", description=desc)
    parser.add_argument("p_ctg", type=str, help="Primary contigs file")
    parser.add_argument("a_ctg_base", type=str, help="Alternative base contigs file")
    parser.add_argument("a_ctg_all", type=str, help="Alternative all contigs file")
    parser.add_argument("--align", dest="align", action="store_true", default=False, help="Whether to align alt to primary (default=False).")
    parser.add_argument("-m", dest="minlength", type=int, default=20, help="Min length of an exact match (default 20).")
    parser.add_argument("-c", dest="minscore", type=int, default=0, help="Min score of an exact match (default 0), exact maches are scored by their length and penalized by the indel they create with respect to previously accepted exact matches.")
    parser.add_argument("-n", dest="minn", type=int, default=2, help="Only align graph on exact matches that occur in at least this many samples.")
    
    args = parser.parse_args()
    
    pctg2seq=dict()
    actg2seq=dict()

    pctg2bases=dict()
    pctg2alts=dict()
    base2alts=dict()

    #read primary contigs into memory
    for name,seq in reveal.fasta_reader(args.p_ctg):
        pctg2seq[name.split()[0]]=seq

    for name,seq in reveal.fasta_reader(args.a_ctg_base):
        name=name.split()[0]
        actg2seq[name]=seq
        pctgname=name.split("-")[0]
        if pctgname in pctg2bases:
            pctg2bases[pctgname].append(name)
        else:
            pctg2bases[pctgname]=[name]

    for name,seq in reveal.fasta_reader(args.a_ctg_all):
        name=name.split()[0]
        actg2seq[name]=seq
        pctgname=name.split("-")[0]
        
        basename=name[:-1]+"0"
        
        if basename in base2alts:
            base2alts[basename].append(name)
        else:
            base2alts[basename]=[name]
        
        if pctgname in pctg2alts:
            pctg2alts[pctgname].append(name)
        else:
            pctg2alts[pctgname]=[name]
    
    if args.align:
        transform_collapse(pctg2seq,actg2seq,pctg2bases,pctg2alts,base2alts,minlength=args.minlength,minn=args.minn,minscore=args.minscore)
    else:
        transform(pctg2seq,actg2seq,pctg2bases,pctg2alts)


#Construct a gfa graph representation, where the associated and primary contigs are aligned using reveal
def transform_collapse(pctg2seq, actg2seq, pctg2bases, pctg2alts, base2alts,minlength=20,minn=2,minscore=0):
    #for every primary contig
    for pctg in pctg2seq:
        a1=pctg2seq[pctg]
        a2=a1
        if pctg in pctg2bases: #are there associated contigs?
            bases=pctg2bases[pctg]
            for base in bases:
                alts=base2alts[base]
                a2=a2.replace(actg2seq[base],actg2seq[alts[0]])
                #print pctg,"replace",base,"by",alts[0]
            assert(a1!=a2)
            G,idx=reveal.align([(pctg,a1),(pctg+"-associated",a2)],minlength=minlength,minn=minn,minscore=minscore)
            reveal.write_gfa(G,idx.T,outputfile=pctg+".gfa")

#Construct a gfa graph representation, where a bubble is formed at every associated contig, without aligning associated contig on top of the primary contig
def transform(pctg2seq,actg2seq,pctg2bases,pctg2alts):

    for pctg in pctg2seq.keys():
        g=nx.DiGraph()
        seq=pctg2seq[pctg]
        base2pos=dict()
        pos2base=dict()
        breaks=[]
        altname=pctg+"_alt"
        primname=pctg
        g.graph['samples']={altname,primname}
        
        if pctg not in pctg2bases:
            #print "No alternative contigs for %s"%pctg
            continue

        for base in pctg2bases[pctg]:
            aseq=actg2seq[base]
            p=seq.find(aseq)
            if p==-1:
                print base,"ERROR: Could not find base sequence",aseq,"in primary contig"
                return
            base2pos[base]=p+len(aseq)
            pos2base[p+len(aseq)]=base
            breaks.append(p)
            breaks.append(p+len(aseq))
        
        breaks.append(len(seq))
        breaks=set(breaks)
        breaks=sorted(list(breaks))
        pp=0
        for p in breaks:
            if p in pos2base: #then node is base
                sample={primname}
            else:
                sample={primname,altname}
            g.add_node(str(p),seq=seq[pp:p],sample=sample,offsets={primname:pp})
            if pp!=0:
                g.add_edge(str(pp),str(p))
            pp=p
        
        for alt in sorted(pctg2alts[pctg]):
            pctg,alti,allele=alt.split('-')

            if int(allele)>1: #expect diploid so skip third alternative allele
                #print "Skipping thrid alt allele"
                continue
            
            base="%s-%s-00"%(pctg,alti)
            
            if alt not in actg2seq:
                print "ERROR", alt, "no link to base"
                continue
            
            g.add_node(alt,seq=actg2seq[alt],sample={altname},offsets={altname:0})
            
            for suc in g.successors(str(base2pos[base])):
                g.add_edge(alt,suc)
            for pre in g.predecessors(str(base2pos[base])):
                g.add_edge(pre,alt)
        
        reveal.write_gml(g,"",outputfile="tmp",partition=False)
        #walk the graph and update offsets
        sg=g.subgraph([n for n,d in g.nodes(data=True) if altname in d['sample']])
        
        o=0
        for node in nx.topological_sort(sg):
            g.node[node]['offsets'][altname]=o
            o+=len(g.node[node]['seq'])
        
        reveal.write_gfa(g,"",outputfile="%s.gfa"%pctg)
        #print "Wrote %s.gfa."%pctg

if __name__ == '__main__':
    main()
