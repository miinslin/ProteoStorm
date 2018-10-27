# -*- coding: utf-8 -*-
import os, time, itertools
from ComputeFDR import Compute_FDR
from operator import itemgetter

def CreateRefinedDBGenus(ProteoStorm_dir, subfoldername,
                         REFINED_DB_PEPLEVEL_FDR,
                         PSlogfile, fastaIDmap, FASTA_dir, save_space):

    print '\tCreating refined database (genera restriction approach)...'
    start_time = time.time()

    max_fidx = max(int(x) for x in fastaIDmap.values())
#    print 'max fasta idx:', max_fidx 
    
    from itertools import izip
    inv_fastaIDmap = dict(izip(fastaIDmap.itervalues( ), fastaIDmap.iterkeys( )))
    
    ## PARAMS
    INITIAL_PEPLEVEL_FDR = 0.01
    
    ABUNDANCE_DECOY_BASED = 1 # only accept taxon with abundance better than the first decoy taxon
    TAXON_ABUNDANCE_THRESHOLD = '' # 0.5%
    TAXON_ABUNDANCE_TYPE = 'PEPTIDE' # either 1% peptides, or psms mapping to peptides. PEPTIDE or PSM
    
    ## DIRECTORIES/FILES
    s1preprocessing_output = os.path.join(ProteoStorm_dir, 'S1_PreprocessingOutput')
    S1_PS_dir = os.path.join(s1preprocessing_output, 'ProteoStorm_input') # contains mappings

    
    S1_OutputFiles = os.path.normpath(os.path.join(ProteoStorm_dir, subfoldername, 'S1_OutputFiles'))
    tsv_out = os.path.normpath(os.path.join(S1_OutputFiles, 'PVAL_computations'))
    S2_InputFiles = os.path.join(ProteoStorm_dir, subfoldername,'S2_InputFiles')
    os.mkdir(S2_InputFiles)
    
    ## INTERMEDIATE FILES
    BESTPEPTIDEforSPECTRUMfile = os.path.join(S1_OutputFiles,'S1_No_cutoff_AllPSMs.txt')
    BESTPSMforPEPTIDEfile = os.path.join(S1_OutputFiles, 'S1_No_cutoff_bestPSMforPeptide.txt')
    psm_taxon_file = os.path.join(S1_OutputFiles,'S1_No_cutoff_taxonIDs.txt')
    
    ## COLUMN Variables/Parameters
    title_column = 1
    unrolled_prot_col = 4
    peptideseq_column = 3 # note should be peptide sequences only (can have mods, will remove)
    decoy_prefix = 'XXX_'
    betweenscans = '5,-1' # ex: 12,1  [col, 1 if greater is betterm -1 if smaller is better]
    FDRcalc = '5,0,-1'  # ex: 13,1,1 [col, uselog, direction]         
    
    decidePSMcolumn = int(betweenscans.split(',')[0]) # MSGFScore used to choose between PSMs to the same scan but different pep seq
    decidePSMdirection = int(betweenscans.split(',')[1]) # 1 if greater is better, -1 if smaller is better hit
    FDRcalc_score_column = int(FDRcalc.split(',')[0])  # score to be used for fdr calculation
    FDR_calc_score_direction = int(FDRcalc.split(',')[2])
    
    specfilename_col = 0
    scan_col = 1
    prepost_pepseq_col = 2   
    pepseq_col = 3 
    protein_col = 4
    rawscore_col = 5
    pvalue_col = 6
    mapping_col = 7  
    
    #################################################################################################################
    
    ## ============= CHOOSE BEST PEPTIDE for spectrum =====================##
    combined_scans = {}
    tsv_out_list = sorted(os.listdir(tsv_out))
    for tsv in tsv_out_list:
        partition_idx = tsv.split('_')[1]
        proteinmappingfile = os.path.join(S1_PS_dir, 'DBPart_'+partition_idx+'.txt') # split by '\t', i==3
        
        protmap = {}
        with open(proteinmappingfile,'r') as mfile:
            for line in mfile:
                splitline = line.strip().split('\t')
                protmap[splitline[1][2:-2]] = splitline[3]
        
        with open(os.path.join(tsv_out,tsv),'r') as infile:
            for line in infile:
                if line[0]!='#':
                    splitlines = line.strip().split('\t') 
                    pepsequence = ''.join([aa for aa in splitlines[peptideseq_column][2:-2] if aa.isalpha()])
                    
                    if pepsequence not in protmap:
                        continue

                    mapping = protmap[pepsequence]
                    protein = splitlines[unrolled_prot_col]
                    title = splitlines[title_column]
                    scannum = title.split('scan=')[1].split('"')[0]
                    specfilename = os.path.splitext(title.split('File:"')[1].split('"')[0])[0] 
                    if specfilename not in combined_scans:
                        combined_scans[specfilename] = {}
                    confidencescore = splitlines[FDRcalc_score_column]
                    decidePSM = splitlines[decidePSMcolumn]
    
                    prepend = [specfilename, scannum, splitlines[peptideseq_column], \
                               pepsequence, protein, decidePSM, confidencescore, mapping]  
    
                    if scannum in combined_scans[specfilename]:
                        PSMovw = float(decidePSM)-float(combined_scans[specfilename][scannum][rawscore_col])
                        # if raw score same between two PSMS from same specfile and scannum
                        if PSMovw==0:
                            if decoy_prefix in combined_scans[specfilename][scannum][protein_col] and decoy_prefix not in protein:
                                combined_scans[specfilename][scannum] = prepend[0:]
                            continue
                        # if PSM has better raw score
                        if PSMovw/abs(PSMovw) == decidePSMdirection:
                            combined_scans[specfilename][scannum] = prepend[0:]
                    elif scannum not in combined_scans[specfilename]:
                        combined_scans[specfilename][scannum] = prepend[0:]               

    outfile = open(BESTPEPTIDEforSPECTRUMfile,'w')
    outfile.write('#'+'\t'.join(['Specfilename','Scannum','Peptide','Pepgroup','Protein','RawScore','pValue'])+'\n')
    # pepgroup is amino acid sequence                     
    peplevel = {}
    unique_peplist = set()
    #'Specfilename','Scannum','Peptide','Pepgroup','Protein','RawScore','pValue', add 'mapping'
    for specf in combined_scans:
        for entry in combined_scans[specf].values():
            outfile.write('\t'.join([z for zi, z in enumerate(entry) if zi != mapping_col])+'\n')
            scannum = entry[scan_col]
            pepgroup = entry[pepseq_col]
            unique_peplist.add((pepgroup, entry[mapping_col]))
            PSMscore = float(entry[pvalue_col])
            protein = entry[protein_col]            
            # choose best psm for pepgroup
            if pepgroup in peplevel:
                PSMovw = PSMscore-float(peplevel[pepgroup][pvalue_col])
                if PSMovw==0:
                    if decoy_prefix in peplevel[pepgroup][protein_col] and decoy_prefix not in protein:
                        peplevel[pepgroup] = entry
                    continue
                if PSMovw/abs(PSMovw) == FDR_calc_score_direction:
                    peplevel[pepgroup] = entry
            else:
                peplevel[pepgroup] = entry
    
    outfile.close()
    unique_peplist = sorted(unique_peplist)                       
    del combined_scans
    
    with open(BESTPSMforPEPTIDEfile,'w') as outfile:
        outfile.write('#'+'\t'.join(['Specfilename','Scannum','Peptide', \
                                     'Pepgroup', 'Protein','RawScore',\
                                     'pValue', 'mappings'])+'\n')
        for pep in peplevel:
            outfile.write('\t'.join(peplevel[pep])+'\n')

    # ============= detemine genera for peptides =====================##
    # includes both decoys and targets
    m_peptides_g = {z:{} for z in [x[0] for x in unique_peplist]}
    for f in xrange(max_fidx+1):
        f_idx = str(f)
        # header: >protein_name \t genus \t taxonID
        # user can parse ncbi taxonomy xml for other taxon level
        with open(os.path.join(FASTA_dir, inv_fastaIDmap[f_idx]+'.fasta'),'r') as infile:
            generalist = [line.strip().split('\t')[1] for line in infile.readlines() if line[0]=='>'] 
        
        # peptide = (pepgroup, mappings)
        for peptide in unique_peplist:
            matches = [z.replace('d_','').split('|')[1] for z in peptide[1].split(';') if z.split('|')[0]==f_idx]
            if matches:
                pepseq = peptide[0]
                # given protein index in chunk, append taxon id to peptide dictionary
                for prot_idx in [int(z) for z in matches]:
                    taxon = generalist[prot_idx]
                    # fasta chunk idx, protein idx
                    cidx_protidx_td = str(f_idx)+'|'+str(prot_idx)
                    if taxon in m_peptides_g[pepseq]:
                        m_peptides_g[pepseq][taxon].add(cidx_protidx_td)
                    else:
                        m_peptides_g[pepseq][taxon] = set([cidx_protidx_td])
                    
    del matches

    # ============= write file with psms and taxon ids =====================##
    # pepgroup does not include charge
    # no cutoff applied yet
    with open(psm_taxon_file,'w') as outfile:
        outfile.write('#'+'\t'.join(['Specfilename','Scannum','Peptide', 'Pepgroup', 'Protein','RawScore','pValue', 'taxonMappings'])+'\n')
        for pepgroup in peplevel:
            taxons_mapping = [t+':'+','.join(m_peptides_g[pepgroup][t]) for t in m_peptides_g[pepgroup]]
            outfile.write('\t'.join([str(x) for xi, x in enumerate(peplevel[pepgroup]) if xi != mapping_col])\
                          +'\t'+';'.join(taxons_mapping)+'\n')

    del m_peptides_g
    del unique_peplist
    del peplevel
    
    #################################################################################################################          
    peptide_psm_nocutoff = []
    with open(psm_taxon_file,'r') as infile:
        for line in infile:
            if line[0]!='#':
                sp = line.strip().split('\t')
                sp = [float(x) if xi == pvalue_col else x for xi, x in enumerate(sp)]
                peptide_psm_nocutoff.append(sp)

    # ============ COMPUTE initial PEP LEVEL FDR ===========================================#            
    # entry: [Specfilename, Scannum, Peptide, Protein, RawScore, pValue, taxonlist]
    #compute peptidelevel FDR
    FDRcompute_targetpass_col = 0
    FDRcompute_decoypass_col = 1
    FDRcompute_scorethreshold = 2
    FDRcompute_FDRscore = 4 

    peptidelevelFDR = Compute_FDR(peptide_psm_nocutoff, INITIAL_PEPLEVEL_FDR, \
                                  pvalue_col, FDR_calc_score_direction, \
                                  protein_col, decoy_prefix)    

    if peptidelevelFDR == 'cannot reach specified FDR':
        PSlogfile.write('Could not reach specified FDR. Choose a different FDR.'+'\n')
        raise ValueError('Could not reach specified FDR. Choose a different FDR.')    
    
    target_peptides = peptidelevelFDR[FDRcompute_targetpass_col]
    decoy_peptides = peptidelevelFDR[FDRcompute_decoypass_col]

    PSlogfile.write(str(len(peptide_psm_nocutoff))+' peptides'+'\n')
    PSlogfile.write(str(len(target_peptides))+' target peptides passing FDR threshold'+'\n')
    PSlogfile.write(str(len(decoy_peptides))+' decoy peptides passing FDR threshold'+'\n')
    PSlogfile.write(str(peptidelevelFDR[FDRcompute_scorethreshold])+' cut off score at '+\
            str(round(100*(float(peptidelevelFDR[FDRcompute_FDRscore])),3))+'% FDR'+'\n')

    all_pep_passing = target_peptides
    if ABUNDANCE_DECOY_BASED==1:
        all_pep_passing.extend(decoy_peptides)
    del target_peptides
    del decoy_peptides
    del peptidelevelFDR

    total_unique_genus_pep = {'uniquegenuspsms':0,'totalpsm':0,'uniquegenuspeptides':0,'totalpeptides':0}
    # counting all psms where peptide passes 1% pep level FDR
    # whether a peptide is target or decoy in our database is already predetermined 
    psms_fromallpassingpep = {x[pepseq_col]:0 for x in all_pep_passing}
    allpassingpeptidelist = set([x[pepseq_col] for x in all_pep_passing])
    #'Specfilename','Scannum','Peptide','Pepgroup','Protein','RawScore','pValue', add 'mapping'
    with open(BESTPEPTIDEforSPECTRUMfile,'r') as infile:
        for line in infile:
            if line[0]!='#':
                entry = line.strip().split('\t')
                pepgroup = entry[pepseq_col]
                if pepgroup in allpassingpeptidelist:
                    psms_fromallpassingpep[pepgroup]+=1

    total_unique_genus_pep['totalpsm'] = sum(psms_fromallpassingpep.values())    
    PSlogfile.write(str(total_unique_genus_pep['totalpsm'])+' number of 1% pooled psm level psms with peptides passing peptide-level FDR'+'\n')
    PSlogfile.write(str(len(psms_fromallpassingpep))+' number of peptides passing 1% pep level FDR'+'\n')

    taxon_pepcounts = {}
    taxon_psmcounts = {}
    with open(os.path.join(S1_OutputFiles, 'S1_passing_peptides.txt'),'w') as all_passing_peptides_file:
        for peptide in all_pep_passing:
            total_unique_genus_pep['totalpeptides']+=1
            all_passing_peptides_file.write(peptide[3]+'\t'+peptide[4]+'\t'+str(peptide[6])+'\n')
            # taxon is the genusname
            # fasta chunk idx| protein idx
            # taxon:0|2342,1|32432;taxon2:3|2343,3|2222
            taxonlist = {z.split(':')[0]:z.split(':')[1].split(',') for z in peptide[mapping_col].split(';')}
            taxon_num = set(taxonlist.keys())
            if len(taxon_num) ==1:
                total_unique_genus_pep['uniquegenuspeptides']+=1
                pepseq = peptide[pepseq_col]
                TAXON = list(taxon_num)[0]
                istarget = bool(decoy_prefix not in peptide[protein_col])
                if istarget ==False:
                    TAXON = decoy_prefix+list(taxon_num)[0]
                
                psm_counts_for_peptide = 1.0*psms_fromallpassingpep[pepseq]
                del psms_fromallpassingpep[pepseq]
                
                if TAXON in taxon_psmcounts:
                    taxon_psmcounts[TAXON] += psm_counts_for_peptide
                    taxon_pepcounts[TAXON]+=1
                    total_unique_genus_pep['uniquegenuspsms']+=psm_counts_for_peptide
                if TAXON not in taxon_psmcounts:
                    taxon_psmcounts[TAXON] = psm_counts_for_peptide
                    taxon_pepcounts[TAXON] = 1
                    total_unique_genus_pep['uniquegenuspsms']+=psm_counts_for_peptide
                    
    PSlogfile.write(str(total_unique_genus_pep['uniquegenuspsms'])+' psms out of '+str(1.0*total_unique_genus_pep['totalpsm'])+\
    ' or'+str(total_unique_genus_pep['uniquegenuspsms']/(1.0*total_unique_genus_pep['totalpsm'])*100)+'% from a single genus...'+'\n')

    PSlogfile.write(str(total_unique_genus_pep['uniquegenuspeptides'])+' peptides out of '+str(1.0*total_unique_genus_pep['totalpeptides'])+\
    ' or'+str(total_unique_genus_pep['uniquegenuspeptides']/(1.0*total_unique_genus_pep['totalpeptides'])*100)+'% from a single genus...'+'\n') 
    
    # =========== USE TAXON ABUNDANCE (Unique PEPTIDES) THRESHOLD ============================#
    TOTAL_uniquetaxonpeps = 1.0*sum(taxon_pepcounts.values())
    TOTAL_uniquetaxonPSMs = 1.0*sum(taxon_psmcounts.values())
    
    sort_taxon_bypvalue = []
    for genus in taxon_pepcounts:
        PEPcounts = 1.0*taxon_pepcounts[genus]
        PSMcounts = 1.0*taxon_psmcounts[genus]
        if TAXON_ABUNDANCE_TYPE=='PEPTIDE':
            abundance = PEPcounts/TOTAL_uniquetaxonpeps
            counts = PEPcounts
        if TAXON_ABUNDANCE_TYPE == 'PSM':
            abundance = PSMcounts/TOTAL_uniquetaxonPSMs
            counts = PSMcounts
        if ABUNDANCE_DECOY_BASED==0:
            if abundance>=TAXON_ABUNDANCE_THRESHOLD:
                sort_taxon_bypvalue.append((genus, counts, abundance))
        if ABUNDANCE_DECOY_BASED==1:
            sort_taxon_bypvalue.append((genus, counts, abundance))
        
    sort_taxon_bypvalue = sorted(sort_taxon_bypvalue, key=itemgetter(2), reverse = True)
    with open(os.path.join(S1_OutputFiles, 'S1_taxon_peplevel_'+str(INITIAL_PEPLEVEL_FDR)+'_all.txt'),'w') as outfile:
        outfile.write('\t'.join(['#Genus','PeptideCount','Abundance'])+'\n')
        for item in sort_taxon_bypvalue:
            outfile.write('\t'.join([str(z) for z in item])+'\n')

    if ABUNDANCE_DECOY_BASED ==1:
        # find abundance for first decoy
        first_decoy_idx = [xi for xi,x in enumerate(sort_taxon_bypvalue) if decoy_prefix in x[0]]
        if first_decoy_idx != []:
            first_decoy_idx = first_decoy_idx[0]
            TAXON_ABUNDANCE_THRESHOLD = sort_taxon_bypvalue[first_decoy_idx][2]
            PSlogfile.write('First decoy '+str(sort_taxon_bypvalue[first_decoy_idx][0])+\
                            ' using abundance cutoff '+str(TAXON_ABUNDANCE_THRESHOLD)+'\n')
            sort_taxon_bypvalue = [x for xi,x in enumerate(sort_taxon_bypvalue) if xi<first_decoy_idx]
            
    if TAXON_ABUNDANCE_TYPE=='PEPTIDE':
        with open(os.path.join(S1_OutputFiles, 'S1_taxon_peplevel_'+str(INITIAL_PEPLEVEL_FDR)+'_PEPabundance'+str(TAXON_ABUNDANCE_THRESHOLD)+'.txt'),'w') as outfile:
            outfile.write('\t'.join(['#Genus','PeptideCount','Abundance'])+'\n')
            for item in sort_taxon_bypvalue:
                outfile.write('\t'.join([str(z) for z in item])+'\n')
                   
    ################################################################################################################## 
    # ============ COMPUTE refined prot database pep level FDR ================================#
    PEPLEVEL_FDR_DB = Compute_FDR(peptide_psm_nocutoff, REFINED_DB_PEPLEVEL_FDR, pvalue_col, FDR_calc_score_direction, protein_col, decoy_prefix)
    if PEPLEVEL_FDR_DB == 'cannot reach specified FDR':
        PSlogfile.write('Could not reach specified FDR for refined protein database.'+'\n')
        raise ValueError('Could not reach specified FDR for refined protein database.')

    all_passing_peptides = PEPLEVEL_FDR_DB[FDRcompute_targetpass_col]
    all_passing_peptides.extend(PEPLEVEL_FDR_DB[FDRcompute_decoypass_col])

    PSlogfile.write(str(len(PEPLEVEL_FDR_DB[FDRcompute_targetpass_col]))+' target peptides passing FDR threshold'+'\n')
    PSlogfile.write(str(len(PEPLEVEL_FDR_DB[FDRcompute_decoypass_col]))+' decoy peptides passing FDR threshold'+'\n')
    PSlogfile.write(str(PEPLEVEL_FDR_DB[FDRcompute_scorethreshold])+' cut off score at '+str(round(100*(float(PEPLEVEL_FDR_DB[FDRcompute_FDRscore])),3))+'% FDR'+'\n')
    del PEPLEVEL_FDR_DB
    
    target_proteins = {str(chidx):set() for chidx in xrange(max_fidx+1)}
    decoy_proteins = {str(chidx):set() for chidx in xrange(max_fidx+1)}
    
    chosen_taxons = set([z[0] for z in sort_taxon_bypvalue])
    
    # ============ INCLUDE all protein mappings from x% pep level psms if protein from chosen taxon ===============#
    for pi, peptide in enumerate(all_passing_peptides):
        istarget = bool(decoy_prefix not in peptide[protein_col])
        # taxon is genusname
        # fasta chunk idx| protein idx
        # taxon:0|2342,1|32432;taxon2:3|2343,3|2222
        taxonlist = {z.split(':')[0]:z.split(':')[1].split(',') for z in peptide[-1].split(';')}
        # take mapping if taxon has genus part of chosen_taxons (x% taxon level)
        mappings = [taxonlist[z] for z in taxonlist if z in chosen_taxons]            
        mappings = list(itertools.chain.from_iterable(mappings))
        for m in mappings:
            cidx = m.split('|')[0]
            protidx = int(m.split('|')[1])
            if istarget:
                target_proteins[cidx].add((pi, protidx))
            else:
                decoy_proteins[cidx].add((pi, protidx))
    
    del all_passing_peptides
    
    # ============ CREATE refined protein DB ========================================================================#
#    secondlevelproteins = set()
    secondlevelproteins = {}
    RefinedDBfile_t_dir = os.path.join(S1_OutputFiles,'RefinedProteinDB_t')
    os.mkdir(RefinedDBfile_t_dir)
    RefinedDBfile_d_dir = os.path.join(S1_OutputFiles,'RefinedProteinDB_d')
    os.mkdir(RefinedDBfile_d_dir)
    RefinedDBfile_target = os.path.join(RefinedDBfile_t_dir,'RefinedProteinDB_target.fasta')
    RefinedDBfile_decoy = os.path.join(RefinedDBfile_d_dir,'RefinedProteinDB_decoy.fasta')
    RefinedDBfile_combined = os.path.join(S1_OutputFiles, 'RefinedProteinDB.fasta')
    write_db_t = open(RefinedDBfile_target, 'w')
    write_db_d = open(RefinedDBfile_decoy, 'w')
    combined_db = open(RefinedDBfile_combined,'w')
    protnum = -1
    
    for f in xrange(max_fidx+1):
        fidx = str(f)
        #read chunk fasta into memory
        with open(os.path.join(FASTA_dir, inv_fastaIDmap[fidx]+'.fasta'),'r') as infile:
            records = ''.join(['$' if line[0]=='>' else line.strip() for line in infile.readlines()])
        records = records.split('$')[1:]

        with open(os.path.join(FASTA_dir, inv_fastaIDmap[fidx]+'.fasta'),'r') as infile:
            records_name = [line.strip().split()[0][1:] for line in infile.readlines() if line[0]=='>']
            
        # targets
        protidexes = [s[1] for s in target_proteins[fidx]]
        
        for entry in protidexes: 
            protseq = records[entry]
            name = records_name[entry]
            decoyseq = protseq[::-1]
            if protseq in secondlevelproteins:
                secondlevelproteins[protseq]['target'].add(name)
                secondlevelproteins[protseq]['decoy'].add(name)
                continue
            #secondlevelproteins.add(protseq)
            secondlevelproteins[protseq] = {'target':set([name]),'decoy':set([name])}

            sequence_split = [protseq[i:i+60] for i in range(0, len(protseq), 60)]
            dec_sequence_split = [decoyseq[i:i+60] for i in range(0, len(decoyseq), 60)]
            protnum+=1
            header = 'Protein_'+str(protnum)
        
            write_db_t.write('>'+header+'\n'+'\n'.join(sequence_split)+'\n')
            write_db_d.write('>XXX_'+header+'\n'+'\n'.join(dec_sequence_split)+'\n')
                   
        del protidexes

        # decoys
        protidexes = [s[1] for s in decoy_proteins[fidx]]
        
        for entry in protidexes:
            protseq = records[entry]
            name = records_name[entry]
            # if protseq in secondlevelproteins, 1) target prot seq already included, decoy already included, 2) no tarrget, but decoy already included
            if protseq in secondlevelproteins:
                secondlevelproteins[protseq]['decoy'].add(name)
                continue
            #secondlevelproteins.add(protseq)     
            secondlevelproteins[protseq]= {'target':set(),'decoy':set([name])}
            
            xxx_protseq = protseq[::-1] #decoy sequence                
            sequence_split = [xxx_protseq[i:i+60] for i in range(0, len(xxx_protseq), 60)]
            protnum+=1
            header = 'Protein_'+str(protnum)
            
            write_db_d.write('>XXX_'+header+'\n'+'\n'.join(sequence_split)+'\n')

        del protidexes
        del records

    write_db_t.close()
    write_db_d.close()
    # write to combined db and include all prot names
    for protseq in secondlevelproteins:
        if secondlevelproteins[protseq]['target']:
            header = '>'+';'.join(secondlevelproteins[protseq]['target'])
            sequence_split = [protseq[i:i+60] for i in range(0, len(protseq), 60)]
            combined_db.write(header+'\n'+'\n'.join(sequence_split)+'\n')
        if secondlevelproteins[protseq]['decoy']:
            xxx_header = '>XXX_'+';XXX_'.join(secondlevelproteins[protseq]['decoy'])
            xxx_protseq = protseq[::-1]
            xxx_sequence_split = [xxx_protseq[i:i+60] for i in range(0, len(xxx_protseq), 60)]
            combined_db.write(xxx_header+'\n'+'\n'.join(xxx_sequence_split)+'\n')
    
    del secondlevelproteins
    combined_db.close()
    
    if save_space ==1:
        os.remove(BESTPEPTIDEforSPECTRUMfile)
        os.remove(BESTPSMforPEPTIDEfile)
        os.remove(os.path.join(S1_OutputFiles, 'S1_passing_peptides.txt'))
    
    return time.time()-start_time