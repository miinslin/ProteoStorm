# -*- coding: utf-8 -*-
import os, re
from ComputeFDR import Compute_FDR

def PeptideLvlFDRout(S2_OutputFiles, S2_peplvl_FDR, ProteoStormLOG):
    
    def checkfloat(number):
        try:
            float(number)
        except ValueError:
            return False
        return True
    
    #['#SpecFilename','Scan','Peptide','Protein','pvalue']
    specfilename_idx = 0
    scannum_idx = 1
    prepost_pepseq_idx = 2   
    protein_idx = 3
    pvalue_idx = 4 
    pepgroup_idx = 5
    
    combined_scans = {}
    with open(os.path.join(S2_OutputFiles,'ProteoStorm_output.txt'),'r') as infile:
        for line in infile:
            if line[0]!='#':
                sp = line.strip().split('\t')
                specfilename = sp[specfilename_idx]
                peptide = sp[prepost_pepseq_idx]
        
                pepsequence = ''.join([aa for aa in peptide[2:-2] if aa.isalpha()])            
                mods = [float(x.replace('+','')) for x in \
                        re.split(r'([-+]?\d+[.]?\d+)', peptide[2:-2]) if checkfloat(x)]
                mods = sum([x for x in mods if x!=57.021])
                pepgroup = pepsequence+'+'+str(mods)
                sp.append(pepgroup)
                
                if specfilename in combined_scans:
                    combined_scans[specfilename].append(sp)
                else:
                    combined_scans[specfilename] = [sp]
    
    #======= choose best psm for peptide ===================#
    decoy_prefix = 'XXX_'
    FDR_calc_score_direction = -1
    FDRcompute_targetpass_col = 0
    FDRcompute_decoypass_col = 1
    
    peplevel = {}
    for specf in combined_scans:
        for entry in combined_scans[specf]:
            pepgroup = entry[pepgroup_idx]
            PSMscore = float(entry[pvalue_idx])
            protein = entry[protein_idx]         
            
            # choose best psm for pepgroup
            if pepgroup in peplevel:
                PSMovw = PSMscore-float(peplevel[pepgroup][pvalue_idx])
                if PSMovw==0:
                    if decoy_prefix in peplevel[pepgroup][protein_idx] and decoy_prefix not in protein:
                        peplevel[pepgroup] = entry
                    continue
                if PSMovw/abs(PSMovw) == FDR_calc_score_direction:
                    peplevel[pepgroup] = entry
            else:
                peplevel[pepgroup] = entry
    
    #======= 1% pep level peptides (pooled) ===================#
    peplistall = list(peplevel.values())
    peplevelFDR = Compute_FDR(peplistall, S2_peplvl_FDR, pvalue_idx, FDR_calc_score_direction, protein_idx, decoy_prefix)
    
    ProteoStormLOG.write('Peptide groups at 1% peptide-level FDR: '+str(len(peplevelFDR[FDRcompute_targetpass_col]))+'\n')
    ProteoStormLOG.write('p-value threshold for 1% peptide-level FDR: '+str(peplevelFDR[2])+'\n')
    
    target_peptides_all = set([z[pepgroup_idx] for z in peplevelFDR[FDRcompute_targetpass_col]])
    
    header = ['#Specfilename', 'Scan','Peptide','Protein','pvalue']
    
    #write peptides passing pep level FDR to file
    with open(os.path.join(S2_OutputFiles,'Pooled_0.01_pepFDR.tsv'),'w') as outfile,\
    open(os.path.join(S2_OutputFiles,'peplevelFDR_0.01_peptides.txt'),'w') as EMoutfile,\
    open(os.path.join(S2_OutputFiles,'peplevelFDR_0.01_luckydecoys.txt'),'w') as withdecoy:
        outfile.write('#'+'\t'.join(header)+'\n')
        peplist = set()
        for targetpsm in peplevelFDR[FDRcompute_targetpass_col]:
            outfile.write('\t'.join(targetpsm)+'\n')
            pepseq_nomods = ''.join([aa for aa in targetpsm[prepost_pepseq_idx][2:-2] if aa.isalpha()])
            peplist.add(pepseq_nomods)
    
        EMoutfile.write('\n'.join(sorted(peplist)))
        
        peplist_decoy = set()
        for decoypsm in peplevelFDR[FDRcompute_decoypass_col]:
            pepseq_nomods = ''.join([aa for aa in decoypsm[prepost_pepseq_idx][2:-2] if aa.isalpha()])
            peplist_decoy.add(pepseq_nomods)
        
        withdecoy.write('\n'.join(sorted(peplist_decoy)))
    
    #=====calculate 1% PSM level fdr for each file (non-pooled) =========#
    with open(os.path.join(S2_OutputFiles,'Pooled_0.01_pepFDR_nonpooled_0.01_psmFDR.tsv'),'w') as outfile:
        outfile.write('#'+'\t'.join(header)+'\n')
        for specf in combined_scans:
            psmlistall = combined_scans[specf]
            PSMlevelFDR = Compute_FDR(psmlistall, S2_peplvl_FDR, pvalue_idx, FDR_calc_score_direction, protein_idx, decoy_prefix)
            if PSMlevelFDR=='cannot reach specified FDR':
                print specf, ' cannot reach specified fdr'
                continue
    
            for targetpsm in PSMlevelFDR[FDRcompute_targetpass_col]:
                if targetpsm[pepgroup_idx] in target_peptides_all:
                    outfile.write('\t'.join(targetpsm)+'\n')
    
    return