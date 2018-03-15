# -*- coding: utf-8 -*-
import os, re, time, itertools
import numpy as np
from Bio import SeqIO
from ComputeFDR import Compute_FDR
from operator import itemgetter
from InSilicoDigest import insilicodigest
from CoreModule1_CreateDBPartitions import CreateDBPartitions
import shutil

def CreateRefinedDBGenus(ProteoStorm_dir, subfoldername,
                         max_massdiff,
                         min_pep_len, max_pep_len,
                         miscleavages,
                         REFINED_DB_PEPLEVEL_FDR,
                         enzyme, PSlogfile,
                         RAMgb, parallel_n, cygwinpath):

    print '\tCreating refined database (genera restriction approach)...'
    start_time = time.time()
    
    ## PARAMS
    INITIAL_PEPLEVEL_FDR = 0.01
    CREATE_DBPARTS_for_SEMITRYP_S2 = 1
#    CREATE_MODIFICATION_DBPARTS = 0
    
    REMOVE_KRP = 1
    ABUNDANCE_DECOY_BASED = 1 # only accept taxon with abundance better than the first decoy taxon
    TAXON_ABUNDANCE_THRESHOLD = '' # 0.5%
    TAXON_ABUNDANCE_TYPE = 'PEPTIDE' # either 1% peptides, or psms mapping to peptides. PEPTIDE or PSM
    
    ## DIRECTORIES/FILES
    s1preprocessing_output = os.path.join(ProteoStorm_dir, 'S1_PreprocessingOutput')
    S1_fastachunk_dir = os.path.join(s1preprocessing_output, 'FastaChunks')
    S1_protein_mappingsdir = os.path.join(s1preprocessing_output, 'ProteinMappings')
    S1_OutputFiles = os.path.normpath(os.path.join(ProteoStorm_dir, subfoldername, 'S1_OutputFiles'))
    tsv_out = os.path.normpath(os.path.join(S1_OutputFiles, 'PVAL_computations'))
    S2_InputFiles = os.path.join(ProteoStorm_dir, subfoldername,'S2_InputFiles')
    matchfilesdir = os.path.normpath(os.path.join(S2_InputFiles ,'MatchFilesDir_genus'))

    if os.path.exists(matchfilesdir):
        if len(os.listdir(matchfilesdir))>0:
            print 'Match files already exist...deleting old files...'
            shutil.rmtree(matchfilesdir)
    if not os.path.exists(matchfilesdir):
        os.makedirs(matchfilesdir)

    if CREATE_DBPARTS_for_SEMITRYP_S2 == 1:
        tempfiledir = os.path.normpath(os.path.join(S2_InputFiles, 'temp_DBfiles'))
        if not os.path.exists(tempfiledir):
            os.makedirs(tempfiledir)
    
    ## CHECK DIRECTORIES
    for dirc in [S1_fastachunk_dir, 
                 S1_protein_mappingsdir,
                 tsv_out]:
        if not dirc:
            PSlogfile.write('Check directory: '+dirc+'\n')
            raise ValueError('Check directory: ', dirc)
    
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
        partition_idx = tsv[:-4].split('_')[-1]
        proteinmappingfile = [z for z in os.listdir(S1_protein_mappingsdir) if z[:-8].split('_')[-1]==partition_idx][0]
        
        protmap = {}
        with open(os.path.join(S1_protein_mappingsdir,proteinmappingfile),'r') as mfile:
            for line in mfile:
                splitline = line.strip().split()
                protmap[splitline[0]] = splitline[1]
        
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

    # ============= local write matching files =====================##
    totalnumchunks = len([x for x in os.listdir(S1_fastachunk_dir) if x.endswith('.fasta') and 'revCat.fasta' not in x])

    for CIDX in range(0,totalnumchunks):
        Fastafile = os.path.join(S1_fastachunk_dir,\
                                 [x for x in os.listdir(S1_fastachunk_dir) if x[:-6].split('_')[-1]==str(CIDX)\
                                  and x.endswith('.fasta') and 'revCat.fasta' not in x][0])     
    
        concatenated_prots = ''.join(['$'+str(record.seq)+'*' for record in SeqIO.parse(Fastafile, "fasta")])
        all_end_pos = np.array([m.start() for m in re.finditer(r'[*]', concatenated_prots)])
    
        with open(os.path.join(matchfilesdir,'combinedDB_part_'+str(CIDX)+'.match'),'w') as outfile:
            for peptide in unique_peplist:
                # already having correct mapping. 
                matches = [z.replace('d_','').split('|')[1] for z in peptide[1].split(';') if z.split('|')[0]==str(CIDX)]
                if not matches:
                    outfile.write('-1'+'\n')
                elif matches:
                    outfile.write('\t'.join([str(np.searchsorted(all_end_pos,int(m),side='right')) for m in matches])+'\n') 
    
    del concatenated_prots
    del all_end_pos
    del matches
    
    # ============= write file with psms and taxon ids =====================##
    # includes both decoys and targets
    m_peptides_g = {z:{} for z in [x[0] for x in unique_peplist]}
    
    for chidx in range(0,totalnumchunks):
        mfile = os.path.join(matchfilesdir, 'combinedDB_part_'+str(chidx)+'.match')
        fastafile = os.path.join(S1_fastachunk_dir,'combinedDB_part_'+str(chidx)+'.fasta')
        
        generalist = [str(record.description).split()[1] for record in SeqIO.parse(fastafile, "fasta")]
        
        with open(mfile,'r') as infile:
            for pi, p in enumerate(infile.readlines()):         
                if p.strip()!= '-1':
                    #unmodified peptide
                    pepseq = unique_peplist[pi][0]
                    # given protein index in chunk, append taxon id to peptide dictionary
                    for prot_idx in [int(z) for z in p.strip().split('\t')]:
                        taxon = generalist[prot_idx]
                        # fasta chunk idx, protein idx
                        cidx_protidx_td = str(chidx)+'|'+str(prot_idx)
                        if taxon in m_peptides_g[pepseq]:
                            m_peptides_g[pepseq][taxon].add(cidx_protidx_td)
                        else:
                            m_peptides_g[pepseq][taxon] = set([cidx_protidx_td])
    
    # pepgroup does not include charge
    with open(psm_taxon_file,'w') as outfile:
        outfile.write('#'+'\t'.join(['Specfilename','Scannum','Peptide', 'Pepgroup', 'Protein','RawScore','pValue', 'taxonMappings'])+'\n')
        for pepgroup in peplevel:
            taxons_mapping = [t+':'+','.join(m_peptides_g[pepgroup][t]) for t in m_peptides_g[pepgroup]]
            outfile.write('\t'.join([str(x) for xi, x in enumerate(peplevel[pepgroup]) if xi != mapping_col])\
                          +'\t'+';'.join(taxons_mapping)+'\n')

    del m_peptides_g
    del unique_peplist
    del peplevel
    shutil.rmtree(matchfilesdir)
    
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
            # taxon:0|2342|len,1|32432|len;taxon2:3|2343|len,3|2222|len
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
    peplist = [z[pepseq_col] for z in all_passing_peptides]

    PSlogfile.write(str(len(PEPLEVEL_FDR_DB[FDRcompute_targetpass_col]))+' target peptides passing FDR threshold'+'\n')
    PSlogfile.write(str(len(PEPLEVEL_FDR_DB[FDRcompute_decoypass_col]))+' decoy peptides passing FDR threshold'+'\n')
    PSlogfile.write(str(PEPLEVEL_FDR_DB[FDRcompute_scorethreshold])+' cut off score at '+str(round(100*(float(PEPLEVEL_FDR_DB[FDRcompute_FDRscore])),3))+'% FDR'+'\n')
    del PEPLEVEL_FDR_DB
    
    target_proteins = {str(chidx):set() for chidx in range(0,totalnumchunks)}
    decoy_proteins = {str(chidx):set() for chidx in range(0,totalnumchunks)}
    
    chosen_taxons = set([z[0] for z in sort_taxon_bypvalue])
    
    # ============ INCLUDE all protein mappings from x% pep level psms if protein from chosen taxon ===============#
    for pi, peptide in enumerate(all_passing_peptides):
        istarget = bool(decoy_prefix not in peptide[protein_col])
        # taxon is genusname
        # taxon:0|2342|t,1|32432|d;taxon2:3|2343|t,3|2222|t
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
    num_proteins_per_peptide = {z:{'d':0,'t':0} for z in range(len(peplist))}  
    secondlevelproteins = set()
    RefinedDBfile = os.path.join(S1_OutputFiles,'SecondLevelDB_INIpeplevel_'+str(INITIAL_PEPLEVEL_FDR)+'_genusabundance_'+str(TAXON_ABUNDANCE_THRESHOLD)+'_peplevel_'+str(REFINED_DB_PEPLEVEL_FDR)+'.fasta')
    write_db = open(RefinedDBfile, 'w')

    enzymerules = {'trypsin':r'([KR*](?=[^P]))'}
    regex_enz = re.compile(enzymerules[enzyme])

    if CREATE_DBPARTS_for_SEMITRYP_S2 == 1:
        with open(os.path.join(tempfiledir,'Unsorted_combined.txt'),'w') as outfile:
            outfile.write('ZZZZEND'+'\n')
           
    for f in [chunk for chunk in os.listdir(S1_fastachunk_dir) if chunk.endswith('.fasta') and 'revCat.fasta' not in chunk]: 

        if CREATE_DBPARTS_for_SEMITRYP_S2 == 1:
            t_sequence = ''
            d_sequence = '' 

        records = list(SeqIO.parse(os.path.join(S1_fastachunk_dir,f), "fasta"))
        chunkID = f[:-6].split('_')[-1]

        # targets
        protidexes = [s[1] for s in target_proteins[chunkID]]
        peptideindexes = [s[0] for s in target_proteins[chunkID]]
        
        for index, entry in enumerate(protidexes): 
            protseq = str(records[entry].seq)
            decoyseq = protseq[::-1]
            if protseq in secondlevelproteins:
                continue
            secondlevelproteins.add(protseq)
            num_proteins_per_peptide[peptideindexes[index]]['t']+=1
            num_proteins_per_peptide[peptideindexes[index]]['d']+=1

            sequence_split = [protseq[i:i+60] for i in range(0, len(protseq), 60)]
            dec_sequence_split = [decoyseq[i:i+60] for i in range(0, len(decoyseq), 60)]
            header = str(records[entry].name)
            
            write_db.write('>'+header+'\n'+'\n'.join(sequence_split)+'\n'+'>XXX_'+header+'\n'+'\n'.join(dec_sequence_split)+'\n')

            if CREATE_DBPARTS_for_SEMITRYP_S2 == 1:
                t_sequence += '$'+protseq+'*'   
                d_sequence += '$'+decoyseq+'*'   
                   
        del protidexes
        
        # decoys
        protidexes = [s[1] for s in decoy_proteins[chunkID]]
        peptideindexes = [s[0] for s in decoy_proteins[chunkID]]
        
        for index, entry in enumerate(protidexes):
            protseq = str(records[entry].seq)
            # if target prot seq already included, decoy already included, so skip
            if protseq in secondlevelproteins:
                continue
            protseq = protseq[::-1]
            if protseq in secondlevelproteins:
                continue
            secondlevelproteins.add(protseq)     
            num_proteins_per_peptide[peptideindexes[index]]['d']+=1                   
            
            sequence_split = [protseq[i:i+60] for i in range(0, len(protseq), 60)]
            header = str(records[entry].name)
            
            write_db.write('>XXX_'+header+'\n'+'\n'.join(sequence_split)+'\n')

            if CREATE_DBPARTS_for_SEMITRYP_S2 == 1:
                d_sequence += '$'+protseq+'*'

        del protidexes
        del records

        if CREATE_DBPARTS_for_SEMITRYP_S2 == 1:
            # $ indicates start, * indicates end of protein
            result = insilicodigest(t_sequence, regex_enz, min_pep_len, max_pep_len, miscleavages,'semi','', REMOVE_KRP)
            del t_sequence
            result_d = insilicodigest(d_sequence, regex_enz, min_pep_len, max_pep_len, miscleavages,'semi','d_', REMOVE_KRP)
            del d_sequence
        
            #write to file
            if os.path.exists(os.path.join(tempfiledir,'Unsorted_combined.txt')):
                outfile = open(os.path.join(tempfiledir,'Unsorted_combined.txt'),'a')
            elif not os.path.exists(os.path.join(tempfiledir,'Unsorted_combined.txt')):
                outfile = open(os.path.join(tempfiledir,'Unsorted_combined.txt'),'w')
    
            for pep in result:
                outfile.write(pep+'\t'+';'.join(['|'+k[0]+'|'+k[1] for k in result[pep]]))
                if pep in result_d:
                    outfile.write(';'+';'.join(['|'+k[0]+'|'+k[1] for k in result_d[pep]]))
                outfile.write('\n')
            for pep in result_d:
                if pep not in result:
                    outfile.write(pep+'\t'+';'.join(['|'+k[0]+'|'+k[1] for k in result_d[pep]])+'\n')
        
            outfile.close()
            del result
            del result_d

    del secondlevelproteins
    write_db.close()

    if CREATE_DBPARTS_for_SEMITRYP_S2 == 1:
        print 'Finished refined database creation...'
        print 'Creating semi-tryptic peptide partitions...'
        CreateDBPartitions(ProteoStorm_dir, subfoldername, miscleavages, 
                                 min_pep_len, max_pep_len, REMOVE_KRP, 
                                 max_massdiff, 'S2', enzyme,
                                 RAMgb, parallel_n, cygwinpath)
    
#    if CREATE_MODIFICATION_DBPARTS ==1:
#        print 'creating modification s2 input files...'
#        from InputFileCreation_ModPeptides import S2_Partitions_Mods
#        print S2_Partitions_Mods(ProteoStorm_dir, subfoldername)
    
    return time.time()-start_time