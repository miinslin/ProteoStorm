# -*- coding: utf-8 -*-
import os, time, re
import shutil
import numpy as np
from Bio import SeqIO
from ComputeFDR import Compute_FDR
from InSilicoDigest import insilicodigest
from CoreModule1_CreateDBPartitions import CreateDBPartitions

def CreateRefinedDBPeptide(ProteoStorm_dir, subfoldername,
                         max_massdiff,
                         min_pep_len, max_pep_len,
                         miscleavages,
                         REFINED_DB_PEPLEVEL_FDR,
                         enzyme, PSlogfile,
                         RAMgb, parallel_n, cygwinpath):

    print '\tCreating refined database...'
    start_time = time.time()
    
    s1preprocessing_output = os.path.join(ProteoStorm_dir, 'S1_PreprocessingOutput')
    S1_fastachunk_dir = os.path.join(s1preprocessing_output,'FastaChunks')
    S1_protein_mappingsdir = os.path.join(s1preprocessing_output,'ProteinMappings')

    CREATE_DBPARTS_for_SEMITRYP_S2 = 1
    REMOVE_KRP = 1

    S1output = os.path.join(ProteoStorm_dir, subfoldername,'S1_OutputFiles')
    BESTPEPTIDEforSPECTRUMfile = os.path.join(S1output,'S1_No_cutoff_AllPSMs.txt')
    BESTPSMforPEPTIDEfile = os.path.join(S1output, 'S1_No_cutoff_bestPSMforPeptide.txt')
    
    tsv_out = os.path.normpath(os.path.join(S1output, 'PVAL_computations'))
    title_column = 1
    unrolled_prot_col = 4
    peptideseq_column = 3 # note should be peptide sequences only (can have mods, will remove)
    decoy_prefix = 'XXX_'
    betweenscans = '5,-1' # ex: 12,1  [col, 1 if greater is betterm -1 if smaller is better]
    FDRcalc = '5,0,-1'  # ex: 13,1,1 [col, uselog, direction]  
    
    S2_InputFiles = os.path.join(ProteoStorm_dir, subfoldername,'S2_InputFiles')
    matchfilesdir = os.path.normpath(os.path.join(S2_InputFiles ,'MatchFilesDir_peplevel'))
    if not os.path.exists(matchfilesdir):
        os.makedirs(matchfilesdir)
    
    if CREATE_DBPARTS_for_SEMITRYP_S2 == 1:
        tempfiledir = os.path.normpath(os.path.join(S2_InputFiles, 'temp_DBfiles'))
        if not os.path.exists(tempfiledir):
            os.makedirs(tempfiledir)

    decidePSMcolumn = int(betweenscans.split(',')[0]) # MSGFScore used to choose between PSMs to the same scan but different pep seq
    decidePSMdirection = int(betweenscans.split(',')[1]) # 1 if greater is better, -1 if smaller is better hit
    FDRcalc_score_column = int(FDRcalc.split(',')[0])  # score to be used for fdr calculation
    FDR_calc_score_direction = int(FDRcalc.split(',')[2])
    totalnumchunks = len([x for x in os.listdir(S1_fastachunk_dir) if x.endswith('.fasta') and 'revCat.fasta' not in x])
    
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
        
        # write matching files 
        for CIDX in range(0,totalnumchunks):
            Fastafile = os.path.join(S1_fastachunk_dir,[x for x in os.listdir(S1_fastachunk_dir) if x[:-6].split('_')[-1]==str(CIDX)\
            and x.endswith('.fasta') and 'revCat.fasta' not in x][0])        
        
            concatenated_prots = ''.join(['$'+str(record.seq)+'*' for record in SeqIO.parse(Fastafile, "fasta")])
            all_end_pos = np.array([m.start() for m in re.finditer(r'[*]', concatenated_prots)])
        
            with open(os.path.join(matchfilesdir,'combinedDB_part_'+str(CIDX)+'.match'),'w') as outfile:
                for pepmappings in mappingslist:
                    # already having correct mapping. 
                    matches = [z.replace('d_','').split('|')[1] for z in pepmappings.split(';') if z.split('|')[0]==str(CIDX)]
                    if not matches:
                        outfile.write('-1'+'\n')
                    elif matches:
                        outfile.write('\t'.join([str(np.searchsorted(all_end_pos,int(m),side='right')) for m in matches])+'\n') 
        
        num_proteins_per_peptide = {z:{'d':0,'t':0} for z in range(len(unique_peplist))}    
        
        del concatenated_prots
        del all_end_pos
        del matches
        
        # ================= CREATE REFINED prot DB ================= #
        secondlevelproteins = set()
        RefinedDBfile = os.path.join(S1output,'SecondLevelDB_peplevel_FDR'+str(REFINED_DB_PEPLEVEL_FDR)+'.fasta')
        write_db = open(RefinedDBfile, 'w')
        
        enzymerules = {'trypsin':r'([KR*](?=[^P]))'}
        regex_enz = re.compile(enzymerules[enzyme])
          
        if CREATE_DBPARTS_for_SEMITRYP_S2 == 1:
            with open(os.path.join(tempfiledir,'Unsorted_combined.txt'),'w') as outfile:
                outfile.write('ZZZZEND'+'\n')
          
        for f in [chunk for chunk in os.listdir(S1_fastachunk_dir) if chunk.endswith('.fasta') and 'revCat.fasta' not in chunk]: 
            #read chunk fasta into memory
            records = list(SeqIO.parse(os.path.join(S1_fastachunk_dir,f), "fasta"))
            chunkID = f[:-6].split('_')[-1]
            #combinedDB_part_0.match
            with open(os.path.join(matchfilesdir,'combinedDB_part_'+chunkID+'.match'),'r') as infile:
                alllines = infile.readlines()
        
            if CREATE_DBPARTS_for_SEMITRYP_S2 == 1:
                t_sequence = ''
                d_sequence = '' 
        
            #targets
            pep_prot_match = [(z.strip().split('\t'), zi) for zi, z in enumerate(alllines) if zi in search_pepindx and z.strip()!='-1']
            
            for pep in pep_prot_match:
                pep_idx = pep[1]
                for entry in [int(i) for i in pep[0]]: 
                    protseq = str(records[entry].seq)
                    decoyseq = protseq[::-1]
                    if protseq in secondlevelproteins:
                        continue
                    secondlevelproteins.add(protseq)
                    num_proteins_per_peptide[pep_idx]['t']+=1
                    num_proteins_per_peptide[pep_idx]['d']+=1
        
                    sequence_split = [protseq[i:i+60] for i in range(0, len(protseq), 60)]
                    dec_sequence_split = [decoyseq[i:i+60] for i in range(0, len(decoyseq), 60)]
                    header = str(records[entry].name)
                    
                    write_db.write('>'+header+'\n'+'\n'.join(sequence_split)+'\n'+'>XXX_'+header+'\n'+'\n'.join(dec_sequence_split)+'\n')
                    
                    if CREATE_DBPARTS_for_SEMITRYP_S2 == 1:
                        t_sequence += '$'+protseq+'*'   
                        d_sequence += '$'+decoyseq+'*'        
            
            del pep_prot_match
            
            #decoys
            pep_prot_match = [(z.strip().split('\t'), zi) for zi, z in enumerate(alllines) if zi in search_pepindx_decoy and z.strip()!='-1']
            
            for pep in pep_prot_match:
                pep_idx = pep[1]
                for entry in [int(i) for i in pep[0]]: 
                    protseq = str(records[entry].seq)
                    # if target prot seq already included, decoy already included, so skip
                    if protseq in secondlevelproteins:
                        continue
                    protseq = protseq[::-1]
                    if protseq in secondlevelproteins:
                        continue
                    secondlevelproteins.add(protseq)
                    num_proteins_per_peptide[pep_idx]['d']+=1
                           
                    sequence_split = [protseq[i:i+60] for i in range(0, len(protseq), 60)]
                    header = str(records[entry].name)
                    
                    write_db.write('>XXX_'+header+'\n'+'\n'.join(sequence_split)+'\n')
                    
                    if CREATE_DBPARTS_for_SEMITRYP_S2 == 1:
                        d_sequence += '$'+protseq+'*'
        
            del pep_prot_match
            del records
            del alllines
            
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
        
#        # write num of proteins per peptide to file
#        with open(os.path.join(S1output, 'SecondLevelDB_peplevel_FDR'+str(REFINED_DB_PEPLEVEL_FDR)+'_protdistribution.txt'),'w') as pepprotdist:
#            pepprotdist.write('#Peptide'+'\t'+'TargetProtNum'+'\t'+'DecoyProtNum'+'\n')
#            for item in num_proteins_per_peptide:
#                pepprotdist.write(unique_peplist[item]+'\t'+str(num_proteins_per_peptide[item]['t'])+'\t'+str(num_proteins_per_peptide[item]['d'])+'\n')
        shutil.rmtree(matchfilesdir)
        del secondlevelproteins
        write_db.close()
            
        if CREATE_DBPARTS_for_SEMITRYP_S2 == 1:
            print 'Finished refined database creation...'
            print 'Creating semi-tryptic peptide partitions...'
            CreateDBPartitions(ProteoStorm_dir, subfoldername,
                                     miscleavages, min_pep_len, max_pep_len,
                                     REMOVE_KRP, max_massdiff, 'S2', enzyme,
                                     RAMgb, parallel_n, cygwinpath)
    
    return time.time()-start_time
