
import utils
import bubbles
import networkx as nx
import os

def unzip(args):
    if not args.graph[0].endswith(".gfa"):
        logging.fatal("Invalid gfa file.")
        return

    G=nx.MultiDiGraph()
    utils.read_gfa(args.graph[0], None, None, G, remap=False)

    if args.source==None and args.sink==None:
        unzip_graph(G,minunzip=args.minunzip)
    else:
        b=bubbles.Bubble(G,args.source,args.sink)
        unzip_bubble(G,b,minunzip=args.minunzip,idoffset=max([n for n in G.nodes() if type(n)==int])+1)

    if args.output==None:
        of=os.path.splitext(args.graph[0])[0]+".unzipped.gfa"
    else:
        of=args.output+".gfa"

    utils.write_gfa(G,None,outputfile=of)

#determine uncertainty about bubble positions
def unzip_graph(G,minunzip=0):
    nid=G.number_of_nodes()
    while nid in G:
        nid+=1

    # for node,data in G.nodes(data=True):
    #     pred=list(G.predecessors(node))
    #     suc=list(G.predecessors(node))

    #     lw=0
    #     l=len(G.node[node]['seq'])

    #     if len(pred)>1 and len(suc)>1:
    #         maxw=int(round((l-2 if l>2 else 0)/float(2)))
    #         lw=rw=maxw if maxw<minunzip else minunzip
    #     elif len(pred)>1 and len(suc)<=1:
    #         lw=l-1
    #     else:
    #         continue

    #     dele=[]
    #     if lw>0:
    #         negpaths=set()
    #         pospaths=set()
    #         ndata=data.copy()
    #         ndata['seq']=G.node[node]['seq'][:lw]
    #         data['seq']=data['seq'][lw:] #update
    #         pn=nid
    #         G.add_node(nid,**ndata)
    #         nid+=1
    #         for n0,n1,d in G.in_edges(node,data=True):
    #             if d['oto']=='-':
    #                 negpaths=negpaths.union(d['paths'])
    #             else:
    #                 pospaths=pospaths.union(d['paths'])
    #                 G.add_edge(n0,pn,**d)
    #                 dele.append((n0,n1))

    #         # for n0,n1,d in G.out_edges(node,data=True):
    #         #     if d['ofrom']=='-':
    #         #         negpaths=negpaths.union(d['paths'])
    #         #         # G.add_edge(pn,n1,**d)
    #         #         dele.append((n0,n1))
    #         #     else:
    #         #         pospaths=pospaths.union(d['paths'])

    #         if len(pospaths)>0:
    #             G.add_edge(pn,node,ofrom='+',oto='+',paths=pospaths)
            
    #         # if len(negpaths)>0:
    #             # G.add_edge(node,pn,ofrom='-',oto='-',paths=negpaths)
    #     G.remove_edges_from(dele)

    for b in bubbles.bubbles(G):
        # if b.issimple():
        # if type(b.source)==str or type(b.sink)==str:
        #     print "skipping"
        #     continue
        nid=unzip_bubble(G,b,minunzip=minunzip,idoffset=nid)

def unzip_bubble(G,b,minunzip=0,idoffset=0):
    
    wiggle=b.getwiggle(minwiggle=minunzip)

    if type(b.sink)==str:
        wiggle=(wiggle[0],0)

    if type(b.source)==str:
        wiggle=(0,wiggle[1])

    if wiggle!=(0,0):
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
                    G.add_node(idoffset,seq=ls if n!=b.sink else ls+rs,offsets={p:(G.node[b.source]['offsets'][p]+srcl)-len(ls) for p in G[b.source][n].values()[0]['paths']})
                    props=G[b.source][n].values()[0].copy() #TODO: consider possibilty of structural variant paths!
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
                    G.add_node(idoffset,seq=rs if n!=b.source else ls+rs,offsets={p:(G.node[b.sink]['offsets'][p])-len(rs) for p in G[n][b.sink].values()[0]['paths']})
                    props=G[n][b.sink].values()[0].copy() #TODO: consider possibilty of structural variant paths!
                    G.remove_edge(n,b.sink)
                    G.add_edge(n,idoffset,**props)
                    G.add_edge(idoffset,b.sink,**props)
                    idoffset+=1
                else:
                    G.node[n]['seq']=G.node[n]['seq']+rs
    else:
        print "Skipping bubble unzip!"

    return idoffset