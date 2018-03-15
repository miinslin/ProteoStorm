# -*- coding: utf-8 -*-
import os, re, time

def S2_PSMs(S2output_dir, PSlogfile, subdir):
    
    PSlogfile.write('Writing PSMs...'+'\n')
    
    start_time = time.time()
    
    tsv_out = os.path.normpath(os.path.join(S2output_dir, 'PVAL_computations'))
    title_column = 1
    unrolled_prot_col = 4
    peptideseq_column = 3 # note should be peptide sequences only (can have mods, will remove)
    decoy_prefix = 'XXX_'
    betweenscans = '5,-1' # ex: 12,1  [col, 1 if greater is betterm -1 if smaller is better]

    decidePSMcolumn = int(betweenscans.split(',')[0]) # computed p-value
    decidePSMdirection = int(betweenscans.split(',')[1]) # 1 if greater is better, -1 if smaller is better hit
    
    protein_idx = 4
    rawscore_idx = 5

    def checkfloat(number):
        try:
            float(number)
        except ValueError:
            return False
        return True
    
    # pick best psm for same specfile scannum
    combined_scans = {}
    tsv_out_list = sorted(os.listdir(tsv_out))
    for tsv in tsv_out_list:        
        with open(os.path.join(tsv_out,tsv),'r') as infile:
            for line in infile:
                if line[0]!='#':
                    splitlines = line.strip().split('\t')
                    title = splitlines[title_column]
                    scannum = title.split('scan=')[1].split('"')[0]
                    specfilename = os.path.splitext(title.split('File:"')[1].split('"')[0])[0]
                    
                    full_pepseq = splitlines[peptideseq_column]
                    
                    if specfilename not in combined_scans:
                        combined_scans[specfilename] = {}
                        
                    psmdecidescore = splitlines[decidePSMcolumn]
                    protein = splitlines[unrolled_prot_col]
                    #strip pre and post aa
                    pepsequence = ''.join([aa for aa in full_pepseq[2:-2] if aa.isalpha()])
                    
                    mods = [float(x.replace('+','')) for x in \
                            re.split(r'([-+]?\d+[.]?\d+)', full_pepseq[2:-2]) if checkfloat(x)]
                    mods = sum([x for x in mods if x!=57.021])
                    
                    prepend = [specfilename, 
                               scannum, 
                               pepsequence+'+'+str(mods), 
                               full_pepseq, 
                               protein, 
                               psmdecidescore]
                        
                    if scannum in combined_scans[specfilename]:
                        PSM_score_diff = float(psmdecidescore)-float(combined_scans[specfilename][scannum][rawscore_idx])
                        if PSM_score_diff==0:
                            if decoy_prefix in combined_scans[specfilename][scannum][protein_idx] and decoy_prefix not in protein:
                                combined_scans[specfilename][scannum] = prepend[0:]
                            continue
                        if PSM_score_diff/abs(PSM_score_diff) == decidePSMdirection:
                            combined_scans[specfilename][scannum] = prepend[0:]
                    elif scannum not in combined_scans[specfilename]:
                        combined_scans[specfilename][scannum] = prepend[0:]
    
    with open(os.path.join(S2output_dir,'ProteoStorm_output.txt'),'w') as outfile:
        outfile.write('\t'.join(['#SpecFilename','Scan','Peptide','Protein','pvalue'])+'\n')
        for specfilename in combined_scans:
            for entry in combined_scans[specfilename].values():
                outfile.write('\t'.join([x for i,x in enumerate(entry) if i!=2])+'\n')
    
    return time.time()-start_time