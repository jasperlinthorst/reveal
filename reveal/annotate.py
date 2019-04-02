
import sys
import os
import logging
import uuid
import gzip


def annotate(args):
    import vcf
    from vcf.parser import _Info as VcfInfo

    vcf_reader = vcf.Reader(filename=args.vcffile)

    aid=0
    # variant2aid={}
    aid2variant={}
    vi=0

    vfile=args.vcffile+".fasta"
    with open(vfile,'w') as v:
        for record in vcf_reader:
            if 'diffsize' in record.INFO:
                if record.INFO['diffsize']>=args.mindiff and record.INFO['diffsize']<args.maxdiff:
                    for i,allele in enumerate(record.alleles):
                        allele=str(allele)
                        aid+=1
                        aid2variant[aid]=(record.CHROM,record.POS,i)
                        # variant2aid[(record.CHROM,record.POS,i)]=aid
                        # v.write(">%s\n"%str(record.CHROM+"_"+str(record.POS)+"_"+str(i)))
                        v.write(">%d\n"%aid)
                        for i in range((len(allele)/50)+1): #write in blocks of 50
                            v.write("%s\n"% allele[i*50:(i+1)*50] )
                    vi+=1

    if vi==0:
        logging.info("No variants in this size range.")
        sys.exit(1)

    vcf_reader = vcf.Reader(filename=args.vcffile)

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
        try:
            os.remove(vfile+".out")
            os.remove(vfile+".tbl")
            os.remove(vfile+".masked")
            os.remove(vfile+".cat")
        except OSError:
            pass

        vcf_reader.infos['repm_rfamily']=VcfInfo('repm_rfamily', 1, 'String', 'Best match for RepeatMasker - Repeat family.',None,None)
        vcf_reader.infos['repm_rtype']=VcfInfo('repm_rtype', 1, 'String', 'Best match for RepeatMasker - Repeat type.',None,None)
        vcf_reader.infos['repm_rcov']=VcfInfo('repm_rcov', 1, 'Float', 'Best match for RepeatMasker - Fraction of the repeat annotation covered.',None,None)
        vcf_reader.infos['repm_qcov']=VcfInfo('repm_qcov', 1, 'Float', 'Best match for RepeatMasker - Fraction of the indel covered.',None,None)
        vcf_reader.infos['repm_allele']=VcfInfo('repm_allele', 1, 'Integer', 'The allele that contains the best RepeatMasker match.',None,None)

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
        try:
            os.remove(vfile+".trf")
        except OSError:
            pass

        vcf_reader.infos['trf_copynumber']=VcfInfo('trf_copynumber', 1, 'Float', 'Best match for TRF - copynumber.',None,None)
        vcf_reader.infos['trf_conssize']=VcfInfo('trf_conssize', 1, 'Integer', 'Best match for TRF - concensus size.',None,None)
        vcf_reader.infos['trf_entropy']=VcfInfo('trf_entropy', 1, 'Float', 'Best match for TRF - entropy.',None,None)
        vcf_reader.infos['trf_pattern']=VcfInfo('trf_pattern', 1, 'String', 'Best match for TRF - pattern.',None,None)
        vcf_reader.infos['trf_start']=VcfInfo('trf_start', 1, 'String', 'Best match for TRF - start position of tr.',None,None)
        vcf_reader.infos['trf_end']=VcfInfo('trf_end', 1, 'String', 'Best match for TRF - end position of tr.',None,None)
        vcf_reader.infos['trf_gccontent']=VcfInfo('trf_gccontent', 1, 'Float', 'Best match for TRF - Fraction of GC bases in the repeat pattern.',None,None)
        vcf_reader.infos['trf_percent_indel']=VcfInfo('trf_percent_indel', 1, 'Integer', 'Best match for TRF - percentage of indels within the aligned repeat pattern.',None,None)
        vcf_reader.infos['trf_percent_match']=VcfInfo('trf_percent_match', 1, 'Integer', 'Best match for TRF - percentage of matches within the aligned repeat pattern.',None,None)
        vcf_reader.infos['trf_allele']=VcfInfo('trf_allele', 1, 'Integer', 'The allele that contains the best TRF match.',None,None)

    #clean up
    os.remove(vfile)
    
    vcf_reader.infos['reveal_type']=VcfInfo('reveal_type', 1, 'String', 'REVEAL\'s best guess at the type of variant.',None,None)

    if args.vcffile.endswith('.vcf.gz'):
        outputfile=gzip.open(args.vcffile[:-7]+".annotated.vcf.gz",'w')
    else:
        outputfile=open(args.vcffile[:-4]+".annotated.vcf",'w')

    vcf_writer = vcf.Writer(outputfile, vcf_reader)

    try:
        for record in vcf_reader:

            if record.INFO['diffsize']>=args.mindiff and record.INFO['diffsize']<args.maxdiff:
                # key=str(record.CHROM+"_"+str(record.POS))
                key=(record.CHROM,record.POS)

                if key in repmd: #we have a repeat masker annotation for this allele
                    record.INFO['repm_rfamily']=repmd[key]['rfamily']
                    record.INFO['repm_rtype']=repmd[key]['rtype']
                    record.INFO['repm_rcov']=repmd[key]['rcov']
                    record.INFO['repm_qcov']=repmd[key]['qcov']
                    record.INFO['repm_allele']=repmd[key]['allele'] #numeric representation of the allele that this annotation was based on

                if key in trfd: #we have a trf annotation for this allele
                    record.INFO['trf_copynumber']=trfd[key]['copynumber']
                    record.INFO['trf_conssize']=trfd[key]['cons_size']
                    record.INFO['trf_entropy']=trfd[key]['entropy']
                    record.INFO['trf_start']=trfd[key]['start']
                    record.INFO['trf_end']=trfd[key]['end']
                    record.INFO['trf_pattern']=trfd[key]['pattern']
                    record.INFO['trf_percent_indel']=trfd[key]['percent_indel']
                    record.INFO['trf_percent_match']=trfd[key]['percent_match']
                    record.INFO['trf_gccontent']=(trfd[key]['G']+trfd[key]['C'])/float(trfd[key]['A']+trfd[key]['C']+trfd[key]['G']+trfd[key]['T'])
                    record.INFO['trf_allele']=trfd[key]['allele'] #numeric representation of the allele that this annotation was based on

                #add custom reveal annotation derived from repeatmasker and trf annotations
                if key in repmd and repmd[key]['rcov']>0.8 and repmd[key]['qcov']>0.8 and repmd[key]['rfamily'].startswith('Satellite'): #repeat annotation and allele have reciprocal overlap >0.8
                    record.INFO['reveal_type']='mei'
                elif key in trfd and (len(trfd[key]['masked'])/record.INFO['diffsize'])>0.5: #no mei, but more than 50% of the indel size can be attributed to tandemly repeated sequence
                    if trfd[key]['cons_size']==1:
                        record.INFO['reveal_type']='homopolymer'
                    elif trfd[key]['cons_size']<=6:
                        record.INFO['reveal_type']='micro-satellite'
                    elif trfd[key]['cons_size']<100:
                        record.INFO['reveal_type']='mini-satellite'
                    elif trfd[key]['cons_size']<1000:
                        record.INFO['reveal_type']='macro-satellite'
                    elif trfd[key]['cons_size']>1000:
                        record.INFO['reveal_type']='mega-satellite'
                else:
                    record.INFO['reveal_type']='other'

            vcf_writer.write_record(record)

    except IOError:
        print "IOError"
        pass


def load_trf_annotations(trffile,aid2variant):
    trfd={}
    pvariant=None
    trfcolnames=['start','end','period_size','copynumber','cons_size','percent_match','percent_indel','score','A','C','G','T','entropy','pattern','masked']
    trfcoltypes=[int,int,int,float,int,int,int,int,int,int,int,int,float,str,str]
    #13 47 1 35.0 1 100 0 70 0 0 0 100 0.00 T TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT

    with open(trffile) as trfout:
        for line in trfout:
            line=line.strip()
            if line[0]=='@':
                aid=int(line[1:])
                allele=aid2variant[aid][2]
                variant=aid2variant[aid][:2]
                continue
            else:
                cols=line.rstrip().split()
                repeat={k:t(v) for k,t,v in zip(trfcolnames[:15],trfcoltypes[:15],cols[:15])}
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













