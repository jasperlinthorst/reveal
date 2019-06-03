
import sys
import os
import logging
import uuid
import gzip

def annotate(args):
    import pysam

    vcf_reader=pysam.VariantFile(args.vcffile,"r")

    aid=0
    aid2variant={}
    vi=0

    vfile=args.vcffile+".fasta"
    with open(vfile,'w') as v:
        for record in vcf_reader:
            if 'reveal_diffsize' in record.info:
                if record.info['reveal_diffsize']>=args.mindiff and record.info['reveal_diffsize']<args.maxdiff:
                    for i,allele in enumerate(record.alleles):
                        allele=str(allele)
                        aid+=1
                        aid2variant[aid]=(record.chrom,record.pos,i,len(allele))
                        # variant2aid[(record.CHROM,record.POS,i)]=aid
                        # v.write(">%s\n"%str(record.CHROM+"_"+str(record.POS)+"_"+str(i)))
                        v.write(">%d\n"%aid)
                        for i in range((len(allele)/50)+1): #write in blocks of 50
                            v.write("%s\n"% allele[i*50:(i+1)*50] )
                    vi+=1

    if vi==0:
        logging.info("No variants in this size range.")
        sys.exit(1)
    
    vcf_reader=pysam.VariantFile(args.vcffile,"r")

    repmd={}
    if args.repm:
        #call repeatmasker
        logging.info("Running repeatmasker...")
        if os.system("RepeatMasker -species %s -pa %d -nolow -nocut %s"%(args.species,args.repmproc,vfile))!=0:
            logging.fatal("RepeatMasker failed, make sure calls to \'RepeatMasker\' are possible.")
            sys.exit(1)
        logging.info("Done.")

        repmd=load_repm_annotations(vfile+".out",aid2variant)

        #clean up
        for ext in [".cat",".cat.all",".out",".tbl",".masked"]:
            try:
                os.remove(vfile+ext)
            except OSError:
                pass

        if 'repm_rfamily' not in vcf_reader.header.info:
            vcf_reader.header.info.add('repm_rfamily', 1, 'String', 'Best match for RepeatMasker - Repeat family.')
        if 'repm_rtype' not in vcf_reader.header.info:
            vcf_reader.header.info.add('repm_rtype', 1, 'String', 'Best match for RepeatMasker - Repeat type.')
        if 'repm_rcov' not in vcf_reader.header.info:
            vcf_reader.header.info.add('repm_rcov', 1, 'Float', 'Best match for RepeatMasker - Fraction of the repeat annotation covered.')
        if 'repm_qcov' not in vcf_reader.header.info:
            vcf_reader.header.info.add('repm_qcov', 1, 'Float', 'Best match for RepeatMasker - Fraction of the indel covered.')
        if 'repm_allele' not in vcf_reader.header.info:
            vcf_reader.header.info.add('repm_allele', 1, 'Integer', 'The allele that contains the best RepeatMasker match.')

    trfd={}
    if args.trf:
        #call trf
        logging.info("Running tandem repeat finder...")
        if os.system("trf %s 2 7 7 80 10 20 500 -ngs -h > %s.trf"%(vfile,vfile))!=0:
            logging.fatal("Tandem Repeat Finder failed, make sure calls to \'trf\' are possible.")
            sys.exit(1)
        logging.info("Done.")

        trfd=load_trf_annotations(vfile+".trf",aid2variant)

        #clean up
        for ext in [".trf"]:
            try:
                os.remove(vfile+ext)
            except OSError:
                pass

        if 'trf_copynumber' not in vcf_reader.header.info:
            vcf_reader.header.info.add('trf_copynumber', 1, 'Float', 'Best match for TRF - copynumber.')
        if 'trf_conssize' not in vcf_reader.header.info:
            vcf_reader.header.info.add('trf_conssize', 1, 'Integer', 'Best match for TRF - concensus size.')
        if 'trf_entropy' not in vcf_reader.header.info:
            vcf_reader.header.info.add('trf_entropy', 1, 'Float', 'Best match for TRF - entropy.')
        if 'trf_pattern' not in vcf_reader.header.info:
            vcf_reader.header.info.add('trf_pattern', 1, 'String', 'Best match for TRF - pattern.')
        if 'trf_start' not in vcf_reader.header.info:
            vcf_reader.header.info.add('trf_start', 1, 'String', 'Best match for TRF - start position of tr.')
        if 'trf_end' not in vcf_reader.header.info:
            vcf_reader.header.info.add('trf_end', 1, 'String', 'Best match for TRF - end position of tr.')
        if 'trf_gccontent' not in vcf_reader.header.info:
            vcf_reader.header.info.add('trf_gccontent', 1, 'Float', 'Best match for TRF - Fraction of GC bases in the repeat pattern.')
        if 'trf_percent_indel' not in vcf_reader.header.info:
            vcf_reader.header.info.add('trf_percent_indel', 1, 'Integer', 'Best match for TRF - percentage of indels within the aligned repeat pattern.')
        if 'trf_percent_match' not in vcf_reader.header.info:
            vcf_reader.header.info.add('trf_percent_match', 1, 'Integer', 'Best match for TRF - percentage of matches within the aligned repeat pattern.')
        if 'trf_cov' not in vcf_reader.header.info:
            vcf_reader.header.info.add('trf_cov', 1, 'Float', 'Best match for TRF - fraction of the allele that is masked by this tandem repeat pattern')
        if 'trf_allele' not in vcf_reader.header.info:
            vcf_reader.header.info.add('trf_allele', 1, 'Integer', 'The allele that contains the best TRF match.')

    #clean up
    os.remove(vfile)
    
    try: #if its already there, leave, just try to relabel
        vcf_reader.header.info.add('reveal_type', 1, 'String', 'REVEAL\'s best guess at the type of variant.')
    except:
        pass

    if args.vcffile.endswith('.gz'):
        outputfile=gzip.open(args.vcffile[:-7]+".annotated"+args.vcffile[-7:],'w')
    else:
        outputfile=open(args.vcffile[:-4]+".annotated"+args.vcffile[-4:],'w')

    vcf_writer = pysam.VariantFile(outputfile, 'w', header=vcf_reader.header)

    try:
        for record in vcf_reader:

            if 'reveal_diffsize' in record.info:
                if record.info['reveal_diffsize']>=args.mindiff and record.info['reveal_diffsize']<args.maxdiff:
                    key=(record.chrom,record.pos)

                    if key in repmd: #we have a repeat masker annotation for this allele
                        record.info['repm_rfamily']=repmd[key]['rfamily']
                        record.info['repm_rtype']=repmd[key]['rtype']
                        record.info['repm_rcov']=repmd[key]['rcov']
                        record.info['repm_qcov']=repmd[key]['qcov']
                        record.info['repm_allele']=repmd[key]['allele'] #numeric representation of the allele that this annotation was based on

                    if key in trfd: #we have a trf annotation for this allele
                        record.info['trf_copynumber']=trfd[key]['copynumber']
                        record.info['trf_conssize']=trfd[key]['cons_size']
                        record.info['trf_entropy']=trfd[key]['entropy']
                        record.info['trf_start']=trfd[key]['start']
                        record.info['trf_end']=trfd[key]['end']
                        record.info['trf_pattern']=trfd[key]['pattern']
                        record.info['trf_percent_indel']=trfd[key]['percent_indel']
                        record.info['trf_percent_match']=trfd[key]['percent_match']
                        record.info['trf_gccontent']=(trfd[key]['G']+trfd[key]['C'])/float(trfd[key]['A']+trfd[key]['C']+trfd[key]['G']+trfd[key]['T'])
                        record.info['trf_cov']=len(trfd[key]['masked'])/float(trfd[key]['allelesize'])
                        record.info['trf_allele']=trfd[key]['allele'] #numeric representation of the allele that this annotation was based on

                    #add custom reveal annotation derived from repeatmasker and trf annotations
                    if key in repmd and repmd[key]['rcov']>0.8 and repmd[key]['qcov']>0.8 and not repmd[key]['rfamily'].startswith('Satellite'): #repeat annotation and allele have reciprocal overlap >0.8
                        record.info['reveal_type']='mei'
                    elif key in trfd and record.info['trf_cov']>0.5: #no mei, but more than 50% of the indel size can be attributed to tandemly repeated sequence
                        if trfd[key]['cons_size']==1:
                            record.info['reveal_type']='homopolymer'
                        elif trfd[key]['cons_size']<=6:
                            record.info['reveal_type']='micro-satellite'
                        elif trfd[key]['cons_size']<100:
                            record.info['reveal_type']='mini-satellite'
                        elif trfd[key]['cons_size']<1000:
                            record.info['reveal_type']='macro-satellite'
                        elif trfd[key]['cons_size']>1000:
                            record.info['reveal_type']='mega-satellite'
                    else:
                        record.info['reveal_type']='other'

            vcf_writer.write(record)

        vcf_writer.close()
        
    except IOError:
        vcf_writer.close()


def load_trf_annotations(trffile,aid2variant):
    trfd={}
    pvariant=None
    trfcolnames=['start','end','period_size','copynumber','cons_size','percent_match','percent_indel','score','A','C','G','T','entropy','pattern','masked']
    trfcoltypes=[str,str,int,float,int,int,int,int,int,int,int,int,float,str,str]
    #13 47 1 35.0 1 100 0 70 0 0 0 100 0.00 T TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT

    with open(trffile) as trfout:
        for line in trfout:
            line=line.strip()
            if line[0]=='@':
                aid=int(line[1:])
                allele=aid2variant[aid][2]
                allelesize=aid2variant[aid][3]
                variant=aid2variant[aid][:2]
                continue
            else:
                cols=line.rstrip().split()
                repeat={k:t(v) for k,t,v in zip(trfcolnames[:15],trfcoltypes[:15],cols[:15])}
                repeat['allelesize']=allelesize
                repeat['allele']=allele
            
            if pvariant==variant:
                if repeat['score']>trfd[variant]['score']:
                    trfd[variant]=repeat
            else:
                trfd[variant]=repeat
            
            pvariant=variant

    return trfd


def load_repm_annotations(repmfile,aid2variant):
    repmd={}
    pvariant=None

    repmcolnames=['score','div','del','ins','qsequence','qbegin','qend','qleft','C_','rtype','rfamily','rbegin','rend','rleft','vid','rcov','qcov']
    with open(repmfile) as repmfile:
        for i in range(3):
            h=repmfile.readline()
        
        for line in repmfile:
            cols=line.split()
            
            repeat={k:v for k,v in zip(repmcolnames,cols)}
            
            repeat['score']=float(repeat['score'])
            
            aid=int(cols[4])
            allele=aid2variant[aid][2]
            variant=aid2variant[aid][:2]
            
            qcov=(int(repeat['qend'])-int(repeat['qbegin']))/float(int(repeat['qleft'][1:-1])+int(repeat['qend']))
            
            if repeat['rbegin'][0]=='(':
                rcov=(int(repeat['rend'])-int(repeat['rleft']))/float(int(repeat['rbegin'][1:-1])+int(repeat['rend']))
            else:
                rcov=(int(repeat['rend'])-int(repeat['rbegin']))/float(int(repeat['rleft'][1:-1])+int(repeat['rend']))
                    
            assert(rcov>0 and rcov<=1)
            assert(qcov>0 and qcov<=1)
            
            repeat['rcov']=rcov
            repeat['qcov']=qcov
            repeat['allele']=allele

            if pvariant==variant:
                if repeat['score']>repmd[variant]['score']:
                    repmd[variant]=repeat
            else:
                repmd[variant]=repeat
            
            pvariant=variant

    return repmd













