import utils
import bubbles
import networkx as nx
import os
import logging

def unzip(args):
    if not args.graph[0].endswith(".gfa"):
        logging.fatal("Invalid gfa file.")
        return

    # G=nx.MultiDiGraph()
    G=nx.DiGraph()
    utils.read_gfa(args.graph[0], None, None, G, remap=False)

    if args.source==None and args.sink==None:
        unzip_graph(G,args,minunzip=args.minunzip)
    else:
        b=bubbles.Bubble(G,args.source,args.sink)
        unzip_bubble(G,b,minunzip=args.minunzip,idoffset=max([n for n in G.nodes() if type(n)==int])+1)

    if args.output==None:
        of=os.path.splitext(args.graph[0])[0]+".unzipped.gfa"
    else:
        of=args.output+".gfa"

    utils.write_gfa(G,None,outputfile=of)

#determine uncertainty about bubble positions
def unzip_graph(G,args,minunzip=0):
    nid=max([n for n in G.nodes() if type(n)==int])
    nid+=1

    for b in bubbles.bubbles(G):
        
        if b.maxsize-b.minsize<args.mindiff:
            logging.debug("Skipping bubble %s, diff between smallest and largest allele (%dbp) is smaller than mindiff=%d."%(str(b.nodes),b.maxsize-b.minsize,args.mindiff))
            continue

        if args.maxdiff and b.maxsize-b.minsize>args.maxdiff:
            logging.debug("Skipping bubble %s, diff between smallest and largest allele (%dbp) is larger than maxdiff=%d."%(str(b.nodes),b.maxsize-b.minsize,args.maxdiff))
            continue

        if isinstance(b,bubbles.Bubble):
            nid=unzip_bubble(G,b,minunzip=minunzip,idoffset=nid)

def unzip_bubble(G,b,minunzip=0,idoffset=0):
    
    wiggle=b.getwiggle(minwiggle=minunzip)

    if type(b.sink)==str:
        wiggle=(wiggle[0],0)

    if type(b.source)==str:
        wiggle=(0,wiggle[1])

    if wiggle!=(0,0):
        logging.debug("Unzipping bubble between %s and %s"%(b.source,b.sink))
        srcl=len(G.node[b.source]['seq'])
        snkl=len(G.node[b.sink]['seq'])
        maxlw=int(round((srcl-2 if srcl>2 else 0)/float(2)))
        maxrw=int(round((snkl-2 if snkl>2 else 0)/float(2)))

        if wiggle[0]>maxlw:
            wiggle=(maxlw,wiggle[1])

        if wiggle[1]>maxrw:
            wiggle=(wiggle[0],maxrw)

        if wiggle[0]>0:
            ls=G.node[b.source]['seq'][-wiggle[0]:]
            assert(G.node[b.source]['seq'][:-wiggle[0]]!="")
            G.node[b.source]['seq']=G.node[b.source]['seq'][:-wiggle[0]]
        else:
            ls=""
        
        if wiggle[1]>0:
            rs=G.node[b.sink]['seq'][:wiggle[1]]
            assert(G.node[b.sink]['seq'][wiggle[1]:]!="")
            G.node[b.sink]['seq']=G.node[b.sink]['seq'][wiggle[1]:]
            G.node[b.sink]['offsets']={k:G.node[b.sink]['offsets'][k]+len(rs) for k in G.node[b.sink]['offsets']}
        else:
            rs=""

        successors=list(G.successors(b.source))
        predecessors=list(G.predecessors(b.sink))

        if ls!="":
            for n in successors:
                if len(list(G.predecessors(n)))>1:
                    # G.add_node(idoffset,seq=ls if n!=b.sink else ls+rs,offsets={p:(G.node[b.source]['offsets'][p]+srcl)-len(ls) for p in G[b.source][n].values()[0]['paths']})
                    G.add_node(idoffset,seq=ls if n!=b.sink else ls+rs,offsets={p:(G.node[b.source]['offsets'][p]+srcl)-len(ls) for p in G[b.source][n]['paths']})
                    # props=G[b.source][n].values()[0].copy() #TODO: consider possibilty of structural variant paths!
                    props=G[b.source][n]
                    G.remove_edge(b.source,n)
                    G.add_edge(b.source,idoffset,**props)
                    G.add_edge(idoffset,n,**props)
                    idoffset+=1
                else:
                    G.node[n]['seq']=ls+G.node[n]['seq']
                    G.node[n]['offsets']={k:G.node[n]['offsets'][k]-len(ls) for k in G.node[n]['offsets']}

        if rs!="":
            for n in predecessors:
                if n==b.source and ls!="":
                    continue #was already handled by looping over successors
                if len(list(G.successors(n)))>1:
                    # G.add_node(idoffset,seq=rs if n!=b.source else ls+rs,offsets={p:(G.node[b.sink]['offsets'][p])-len(rs) for p in G[n][b.sink].values()[0]['paths']})
                    G.add_node(idoffset,seq=rs if n!=b.source else ls+rs,offsets={p:(G.node[b.sink]['offsets'][p])-len(rs) for p in G[n][b.sink]['paths']})
                    # props=G[n][b.sink].values()[0].copy() #TODO: consider possibilty of structural variant paths!
                    props=G[n][b.sink]
                    G.remove_edge(n,b.sink)
                    G.add_edge(n,idoffset,**props)
                    G.add_edge(idoffset,b.sink,**props)
                    idoffset+=1
                else:
                    G.node[n]['seq']=G.node[n]['seq']+rs

    return idoffset