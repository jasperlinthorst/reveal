#!/usr/bin/env python
import argparse
import sys, errno
import os

def main():
    desc="""
    Script that converts fasta files to fasta files with daligner compatible naming. Eg. >name/well/start_end
    """
    
    parser = argparse.ArgumentParser(prog="dformat", usage="dformat -h for usage", description=desc)
    parser.add_argument("fasta", type=str, help="Fasta file for which names have to be daligner compatible.")
    
    args = parser.parse_args()
    c=0
    
    template=os.path.basename(args.fasta).replace(".fasta","").replace(".fa","").replace(".fna","").replace(" ","").replace(".","")
    
    try:

        name=None
        for line in open(args.fasta,'r'):
            if line[0]==">":
                if name!=None:
                    l=len(seq)
                    readname=">%s/%s/%d_%d/%d\n"%(template,c,0,l,l)
                    sys.stdout.write(readname)
                    for i in range(0,len(seq),100):
                        if i+100<len(seq):
                            sys.stdout.write(seq[i:i+100]+"\n")
                        else:
                            sys.stdout.write(seq[i:len(seq)]+"\n")
                
                name=line
                c+=1
                seq=""
            else:
                seq+=line.rstrip()

        if name!=None:
            l=len(seq)
            readname=">%s/%s/%d_%d/%d\n"%(template,c,0,l,l)
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
