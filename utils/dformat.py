#!/usr/bin/env python
import argparse
import sys, errno

def main():
    desc="""
    Script that converts fasta files to fasta files with daligner compatible naming. Eg. >name/well/start_end
    """
    
    parser = argparse.ArgumentParser(prog="dformat", usage="dformat -h for usage", description=desc)
    parser.add_argument("fasta", type=str, help="Fasta file for which names have to be daligner compatible.")
    
    args = parser.parse_args()
    c=0
    
    try:

        name=None
        for line in open(args.fasta,'r'):
            if line[0]==">":
                if name!=None:
                    readname=">%s/%s/%d_%d\n"%(name,c,0,len(seq))
                    sys.stdout.write(readname)
                    for i in range(0,len(seq),100):
                        if i+100<len(seq):
                            sys.stdout.write(seq[i:i+100]+"\n")
                        else:
                            sys.stdout.write(seq[i:len(seq)]+"\n")
                
                name=line[1:].rstrip().replace("/","_").replace(" ","_").replace(".","_")
                c+=1
                seq=""
            else:
                seq+=line.rstrip()

        if name!=None:
            readname=">%s/%s/%d_%d\n"%(name,c,0,len(seq))
            sys.stdout.write(readname)
            for i in range(0,len(seq),100):
                if i+100<len(seq):
                    sys.stdout.write(seq[i:i+100]+"\n")
                else:
                    sys.stdout.write(seq[i:len(seq)]+"\n")

    except IOError as e:
        if e.errno == errno.EPIPE:
            pass

if __name__ == '__main__':
    main()
