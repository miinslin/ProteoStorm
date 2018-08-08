# -*- coding: utf-8 -*-
import os, time
from ComputeFDR import Compute_FDR

def CreateRefinedDBPeptide(ProteoStorm_dir, subfoldername,
                         REFINED_DB_PEPLEVEL_FDR,
                         PSlogfile, fastaIDmap, FASTA_dir, save_space):

    start_time = time.time()

    max_fidx = max(int(x) for x in fastaIDmap.values())
#    print 'max fasta idx:', max_fidx 
    
    from itertools import izip
    inv_fastaIDmap = dict(izip(fastaIDmap.itervalues( ), fastaIDmap.iterkeys( )))
    
    s1preprocessing_output = os.path.join(ProteoStorm_dir, 'S1_PreprocessingOutput')
    S1_PS_dir = os.path.join(s1preprocessing_output, 'ProteoStorm_input') # contains mappings

    S1output = os.path.join(ProteoStorm_dir, subfoldername,'S1_OutputFiles')
    tsv_out = os.path.normpath(os.path.join(S1output, 'PVAL_computations'))
    BESTPEPTIDEforSPECTRUMfile = os.path.join(S1output,'S1_No_cutoff_AllPSMs.txt')
    BESTPSMforPEPTIDEfile = os.path.join(S1output, 'S1_No_cutoff_bestPSMforPeptide.txt')

    title_column = 1
    unrolled_prot_col = 4
    peptideseq_column = 3 # note should be peptide sequences only (can have mods, will remove)
    decoy_prefix = 'XXX_'
    betweenscans = '5,-1' # ex: 12,1  [col, 1 if greater is betterm -1 if smaller is better]
    FDRcalc = '5,0,-1'  # ex: 13,1,1 [col, uselog, direction]  
    
    S2_InputFiles = os.path.join(ProteoStorm_dir, subfoldername,'S2_InputFiles')
    os.mkdir(S2_InputFiles)
        
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

    ###############################################################
    if not os.path.exists(BESTPSMforPEPTIDEfile):
        ## ============= choose best peptide for spectra =====================##   
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
                        
                        # msgf+ can assign wrong ID. Issue with reading input fasta file sequence and creating additional non-existing peptides...
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
                        
                        prepend = [specfilename, scannum, splitlines[peptideseq_column], pepsequence, protein, decidePSM, confidencescore, mapping]  
        
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
        
        with open(BESTPEPTIDEforSPECTRUMfile,'w') as outfile:
            outfile.write('#'+'\t'.join(['Specfilename','Scannum','Peptide','Pepgroup','Protein','RawScore','pValue'])+'\n')
            for specfilename in combined_scans:
                for entry in combined_scans[specfilename].values():
                    outfile.write('\t'.join([x for xi, x in enumerate(entry) if xi!=mapping_col])+'\n')
    
        ## ============= choose best psm for peptide =====================##
        peplevel = {}
        # group by peptide and report only best PSM (by specEvalue)
        for specf in combined_scans:
            for entry in combined_scans[specf].values():
                pepgroup= entry[pepseq_col]
                PSMscore = float(entry[pvalue_col])
                protein = entry[protein_col]
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
        del combined_scans
    
        with open(BESTPSMforPEPTIDEfile,'w') as outfile:
            outfile.write('#'+'\t'.join(['Specfilename','Scannum','Peptide', 'Pepgroup', 'Protein','RawScore','pValue', 'mappings'])+'\n')
            for pep in peplevel:
                outfile.write('\t'.join(peplevel[pep])+'\n')
    
    ##############################################################################################################################################   
    if os.path.exists(BESTPSMforPEPTIDEfile):
        # ================= Compute Peptide level FDR ================= #
        FDRcompute_targetpass_col = 0
        FDRcompute_decoypass_col = 1
        FDRcompute_scorethreshold = 2
        FDRcompute_FDRscore = 4 

        peplevel = []
        with open(BESTPSMforPEPTIDEfile,'r') as infile:
            for line in infile:
                if line[0]!='#':
                    sp = line.strip().split('\t')
                    # make pvalue a float
                    sp = [float(x) if xi ==pvalue_col else x for xi, x in enumerate(sp)]
                    peplevel.append(sp)       
        
        #psmlist, user_FDR ,PSMscore_col ,FDR_calc_score_direction, prot_col,decoy_prefix
        peplevelFDR = Compute_FDR(peplevel, REFINED_DB_PEPLEVEL_FDR, pvalue_col, FDR_calc_score_direction, protein_col, decoy_prefix)
        if peplevelFDR == 'cannot reach specified FDR':
            PSlogfile.write('Could not reach specified FDR for refined protein database.'+'\n')
            raise ValueError('Could not reach specified FDR for refined protein database.')
        del peplevel
    
        unique_peplist = set([(m[pepseq_col],m[mapping_col],'t') for m in peplevelFDR[FDRcompute_targetpass_col]])
        unique_peplist.update([(m[pepseq_col],m[mapping_col],'d') for m in peplevelFDR[FDRcompute_decoypass_col]])
        unique_peplist = sorted(list(unique_peplist),key=lambda x: x[0])
    
        # ================= Write peptides passing fdr cutoff score to file ==============================#
        search_pepindx = set()
        search_pepindx_decoy = set()
        peplistoutfile = os.path.join(S1output,'S1_uniquepeplist_'+str(REFINED_DB_PEPLEVEL_FDR)+'.txt')
        pepmappingsfile = os.path.join(S1output,'S1_uniquepeplist_Mappings_'+str(REFINED_DB_PEPLEVEL_FDR)+'.txt')
        with open(peplistoutfile,'w') as peplist_out, open(pepmappingsfile,'w') as mappings_out:
            for upepi, upep in enumerate(unique_peplist):
                peplist_out.write(upep[0]+'\t'+upep[2]+'\n')
                mappings_out.write(upep[1]+'\n')
                #if target
                if upep[2] =='t':
                    search_pepindx.add(upepi)
                #if decoy
                elif upep[2]=='d':
                    search_pepindx_decoy.add(upepi)
        
        PSlogfile.write('FDR cutoff score '+str(peplevelFDR[FDRcompute_scorethreshold])+'\n')
        del peplevelFDR
        unique_peplist = [x[0] for x in unique_peplist]
        
        # ================= CREATE MATCH FILES ================= ##
        with open(pepmappingsfile,'r') as infile:
            mappingslist = [z.strip() for z in infile.readlines()]
        
        protmapfile = open(os.path.join(S2_InputFiles,'PepProt.match'),'w')
        # write matching file - 
        for f in xrange(max_fidx+1):
            f_idx = str(f)
            protmapfile.write(f_idx+'\t')
            for idx,pepmappings in enumerate(mappingslist):
                #6476|d_785|RL
                matches = [z.replace('d_','').split('|')[1] for z in pepmappings.split(';') if z.split('|')[0]==f_idx]
                if not matches:
                    protmapfile.write('-1')
                    if idx != len(mappingslist)-1:
                        protmapfile.write('\t')
                elif matches:
                    protmapfile.write(','.join(matches))
                    if idx != len(mappingslist)-1:
                        protmapfile.write('\t')
            protmapfile.write('\n')
        
        protmapfile.close()
        del matches
        
        # ================= CREATE REFINED prot DB ================= #
        secondlevelproteins = set()
        RefinedDBfile_t_dir = os.path.join(S1output,'RefinedProteinDB_t')
        os.mkdir(RefinedDBfile_t_dir)
        RefinedDBfile_d_dir = os.path.join(S1output,'RefinedProteinDB_d')
        os.mkdir(RefinedDBfile_d_dir)
        RefinedDBfile_target = os.path.join(RefinedDBfile_t_dir,'RefinedProteinDB_target.fasta')
        RefinedDBfile_decoy = os.path.join(RefinedDBfile_d_dir,'RefinedProteinDB_decoy.fasta')
        RefinedDBfile_combined = os.path.join(S1output, 'RefinedProteinDB.fasta')
        write_db_t = open(RefinedDBfile_target, 'w')
        write_db_d = open(RefinedDBfile_decoy, 'w')
        combined_db = open(RefinedDBfile_combined,'w')
        protnum = -1
        
        num_proteins_per_peptide = {z:{'d':0,'t':0} for z in range(len(unique_peplist))}
        
        with open(os.path.join(S2_InputFiles,'PepProt.match'),'r') as infile:
            for line in infile:
                sp = line.strip().split('\t')
                f = inv_fastaIDmap[sp[0]]
                #read chunk fasta into memory
                with open(os.path.join(FASTA_dir,f+'.fasta'),'r') as infile:
                    records = ''.join(['$' if line[0]=='>' else line.strip() for line in infile.readlines()])
                records = records.split('$')[1:]
                alllines = sp[1:]
            
                #targets
                pep_prot_match = [(z.split(','), zi) for zi, z in enumerate(alllines) if zi in search_pepindx and z.strip()!='-1']
                
                for pep in pep_prot_match:
                    pep_idx = pep[1]
                    for entry in [int(i) for i in pep[0]]: 
                        protseq = records[entry]
                        decoyseq = protseq[::-1]
                        if protseq in secondlevelproteins:
                            continue
                        secondlevelproteins.add(protseq)
                        num_proteins_per_peptide[pep_idx]['t']+=1
                        num_proteins_per_peptide[pep_idx]['d']+=1
            
                        sequence_split = [protseq[i:i+60] for i in range(0, len(protseq), 60)]
                        dec_sequence_split = [decoyseq[i:i+60] for i in range(0, len(decoyseq), 60)]
                        protnum+=1
                        header = 'Protein_'+str(protnum)
                        
                        write_db_t.write('>'+header+'\n'+'\n'.join(sequence_split)+'\n')
                        write_db_d.write('>XXX_'+header+'\n'+'\n'.join(dec_sequence_split)+'\n')
                        combined_db.write('>'+header+'\n'+'\n'.join(sequence_split)+'\n'+'>XXX_'+header+'\n'+'\n'.join(dec_sequence_split)+'\n')
                
                del pep_prot_match
                
                #decoys
                pep_prot_match = [(z.split(','), zi) for zi, z in enumerate(alllines) if zi in search_pepindx_decoy and z.strip()!='-1']
                
                for pep in pep_prot_match:
                    pep_idx = pep[1]
                    for entry in [int(i) for i in pep[0]]: 
                        protseq = records[entry]
                        # if target prot seq already included, decoy already included, so skip
                        if protseq in secondlevelproteins:
                            continue
                        protseq = protseq[::-1]
                        if protseq in secondlevelproteins:
                            continue
                        secondlevelproteins.add(protseq)
                        num_proteins_per_peptide[pep_idx]['d']+=1
                               
                        sequence_split = [protseq[i:i+60] for i in range(0, len(protseq), 60)]
                        protnum+=1
                        header = 'Protein_'+str(protnum)
                        
                        write_db_d.write('>XXX_'+header+'\n'+'\n'.join(sequence_split)+'\n')
                        combined_db.write('>XXX_'+header+'\n'+'\n'.join(sequence_split)+'\n')
                
                del pep_prot_match

        del secondlevelproteins
        write_db_t.close()
        write_db_d.close()
        combined_db.close()
        os.remove(os.path.join(S2_InputFiles,'PepProt.match'))
        
        if save_space ==1:
            os.remove(BESTPEPTIDEforSPECTRUMfile)
            os.remove(BESTPSMforPEPTIDEfile)
            os.remove(peplistoutfile)
            os.remove(pepmappingsfile)
        
    return time.time()-start_time
