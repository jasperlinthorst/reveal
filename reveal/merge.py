
from utils import *

def merge_cmd(args):
    if len(args.graphs)<2:
        logging.fatal("Specify multiple gfa files to merge them.")
        return
    
    G=nx.DiGraph()
    for graph in args.graphs:
        logging.info("Adding %s ..." %graph)
        read_gfa(graph,None,"",G,remap=True)
    
    if args.outprefix!=None:
        write_gfa(G,"",outputfile=args.outprefix+".gfa")
    else:
        write_gfa(G,"",outputfile="_".join([os.path.basename(f)[:os.path.basename(f).rfind('.')] for f in args.graphs])+".gfa")
