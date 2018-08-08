# -*- coding: utf-8 -*-
import time
from operator import itemgetter
from os import path, listdir, stat
from bisect import bisect_left
from sys import getsizeof
import shutil

aa_monoisotopic_masses = {'A': 71.037113805,
 'C': 103.009184505 + 57.021463735,
 'D': 115.026943065,
 'E': 129.042593135,
 'F': 147.068413945,
 'G': 57.021463735,
 'H': 137.058911875,
 'I': 113.084064015,
 'K': 128.09496305,
 'L': 113.084064015,
 'M': 131.040484645,
 'N': 114.04292747,
 'P': 97.052763875,
 'Q': 128.05857754,
 'R': 156.10111105,
 'S': 87.032028435,
 'T': 101.047678505,
 'V': 99.068413945,
 'W': 186.07931298,
 'Y': 163.063328575,
 'X':0,
 'B':0,
 'J':0,
 'O':0,
 'U':0,
 'Z':0,
 '?':0,
 '-':0,
 '*':0}

HOH = 18.0105647
TMT_monoisotopicmass = 229.162932178

def calcmass_cmm(pepstring, TMT_mod):
    if TMT_mod ==0:
        return [aa_monoisotopic_masses[aa] for aa in pepstring]
    # TMT tags composed of amine-reactive NHS-ester group, a spacer arm and a mass reporter. 
    # Lysine side-chain NH2 or n-terminus NH2
    # if C at n-terminus is carbamidomethylated, TMT can still react with NH2
    elif TMT_mod ==1:
        TMTmodssum = []
        for i, aa in enumerate(pepstring):
            # if n-terminus of peptide, add 229.162932
            if i==0:
                # if n-terminus aa is also a K, add 229.162932 again (lysine side-chain)
                if aa == 'K':
                    TMTmodssum.append(aa_monoisotopic_masses[aa]+TMT_monoisotopicmass+TMT_monoisotopicmass)
                else:
                    TMTmodssum.append(aa_monoisotopic_masses[aa]+TMT_monoisotopicmass)
            # if not aa at n-terminus of peptide but a K, add 229.162932
            else:
                if aa == 'K':
                    TMTmodssum.append(aa_monoisotopic_masses[aa]+TMT_monoisotopicmass)
                else:
                    TMTmodssum.append(aa_monoisotopic_masses[aa])
        return TMTmodssum


def CopyPeptideS1(pepseq, info, M_ranges, M, C, perbin_bytes, output_directory, 
                  TMT_labeling):
    
    if any(aa in pepseq for aa in 'XBJOUZ?-*'):
        return True
    else:
        pepmass = round(sum(calcmass_cmm(pepseq, TMT_labeling))+HOH,4) # compute mass
        bin_i = bisect_left(M_ranges, pepmass)
        M[bin_i].append(pepseq.replace('L','I')+'\t'+str(pepmass)+'\t'+info+'\n')
        C[bin_i] += getsizeof(pepseq) + 8
        
        if C[bin_i]>perbin_bytes:
            outfile = open(path.join(output_directory, 'DBPart_'+str(bin_i)+'.txt'),'a')
            outfile.write(''.join(M[bin_i]))
            outfile.close()
            M[bin_i] = []
            C[bin_i] = getsizeof([])
            
        return False

def CopyPeptideS2(pepseq, info, M_ranges, M, C, perbin_bytes, output_directory,
                min_pep_len, max_pep_len, leftanchor, rightanchor, TMT_labeling):

    pepseq = pepseq.replace('L','I')
    peplength = len(pepseq)
    pepmass = calcmass_cmm(pepseq, TMT_labeling)

    for i in xrange(peplength):
        #if left anchor, and pepseq[0]=='M' and pre == '-', need to make semi tryptic peps with M removed
        if leftanchor == True:
            if pepseq[0] == 'M' and info[-2]=='-':
                if min_pep_len<=(peplength-1)-i<=max_pep_len:
                    l_spep = pepseq[1:peplength-i] 
                    if any(aa in l_spep for aa in 'XBJOUZ?-*'):
                        l_spep= ''
                    if l_spep:
                        l_spep_mass = round(sum(pepmass[1:peplength-i])+HOH, 4)
                        bin_i = bisect_left(M_ranges, l_spep_mass)
                        if i!=0: 
                            lspep_info = info[:-2]+'-' #preaa should always be -
                            lspep_info += pepseq[peplength-i]
                        else: #0 is the M-cleaved peptide, do not need to update info
                            lspep_info = info[0:]
                                
                        M[bin_i].append(l_spep+'\t'+str(l_spep_mass)+'\t'+lspep_info+'\n')
                        C[bin_i] += getsizeof(l_spep) + 8
                        
                        if C[bin_i]>perbin_bytes:
                            outfile = open(path.join(output_directory, 'DBPart_'+str(bin_i)+'.txt'),'a')
                            outfile.write(''.join(M[bin_i]))
                            outfile.close()
                            M[bin_i] = []
                            C[bin_i] = getsizeof([])
        
        if min_pep_len<=peplength-i<=max_pep_len:
            l_spep = ''
            r_spep = ''
        
            if leftanchor == True:
                l_spep = pepseq[0:peplength-i] 
                if any(aa in l_spep for aa in 'XBJOUZ?-*'):
                    l_spep= ''
                if l_spep:
                    l_spep_mass = round(sum(pepmass[0:peplength-i])+HOH, 4)
                    bin_i = bisect_left(M_ranges, l_spep_mass)
                    
                    if i!=0: #0 is the fully-tryptic peptide, do not need to update info
                        lspep_info = info[:-1] #preaa stays the same
                        lspep_info += pepseq[peplength-i]
                    else:
                        lspep_info = info[0:]
                        
                    M[bin_i].append(l_spep+'\t'+str(l_spep_mass)+'\t'+lspep_info+'\n')
                    C[bin_i] += getsizeof(l_spep) + 8
                    
                    if C[bin_i]>perbin_bytes:
                        outfile = open(path.join(output_directory, 'DBPart_'+str(bin_i)+'.txt'),'a')
                        outfile.write(''.join(M[bin_i]))
                        outfile.close()
                        M[bin_i] = []
                        C[bin_i] = getsizeof([])
                    
            if rightanchor == True:
                if i>0:
                    r_spep = pepseq[i:peplength]
                    if any(aa in r_spep for aa in 'XBJOUZ?-*'):
                        r_spep = ''
                    if r_spep:
                        r_spep_mass = round(sum(pepmass[i:peplength])+HOH, 4)
                        bin_i = bisect_left(M_ranges, r_spep_mass)
                        if pepseq[0] == 'M' and info[-2]=='-' and i==1:
                            rspep_info = info[:-2]+'-'+info[-1] 
                        else:
                            rspep_info = info[:-2]+pepseq[i-1]+info[-1] #postaa stays the same
                        M[bin_i].append(r_spep+'\t'+str(r_spep_mass)+'\t'+rspep_info+'\n')
                        C[bin_i] += getsizeof(r_spep) + 8           
                        if C[bin_i]>perbin_bytes:
                            outfile = open(path.join(output_directory, 'DBPart_'+str(bin_i)+'.txt'),'a')
                            outfile.write(''.join(M[bin_i]))
                            outfile.close()
                            M[bin_i] = []
                            C[bin_i] = getsizeof([])
            
    return 

def TrypticDigest(prot, fasta_protcount,
                  min_pep_len, max_pep_len, 
                  miscleavage, M_ranges, 
                  output_directory, massbin_n, 
                  perbin_bytes, td, reverseseq, stage, TMT_labeling):

    list_bytes = getsizeof([])
    # C[i] is the current number of entries for ith bin
    C = [list_bytes for x in xrange(massbin_n)]
    # M is peptide buffer
    M = [[] for x in xrange(massbin_n)]
    b=0 # previous cutsite
    a=0 # cutsite before b

    if stage == 'S1':
        protnum = -1
        f_idx = 0
        fasta_ID = fasta_protcount[f_idx][0]
        currentmax = fasta_protcount[f_idx][1]
        for i in xrange(len(prot)-1):
            if prot[i]=='$':
                protnum+=1
                if protnum>=currentmax:
                    f_idx+=1
                    fasta_ID = fasta_protcount[f_idx][0]
                    currentmax = fasta_protcount[f_idx][1]
                    protnum = 0
            if (prot[i] in ['K','R'] and prot[i+1]!='P') or prot[i+1]=='$':
                nonstandardaa = 'not assigned'                       
                
                # peptide @ prot[b+1:i+1]
                if min_pep_len<=i-b<=max_pep_len:
                    peptide = prot[b+1:i+1]
                    if td==True:
                        if reverseseq ==True:
                            d_protnum = (currentmax-1)-protnum
                        else:
                            d_protnum = protnum
                        add_info = fasta_ID+'|'+'d_'+str(d_protnum)+'|'+(prot[b]+prot[i+1]).replace('$','-')
                    else:
                        add_info = fasta_ID+'|'+str(protnum)+'|'+(prot[b]+prot[i+1]).replace('$','-')
                    nonstandardaa = CopyPeptideS1(peptide, add_info, M_ranges, 
                                            M, C, perbin_bytes, output_directory, TMT_labeling)
                
                if nonstandardaa!=True:
                    # N-term Methionine cleaved peptide @ prot[b+2:i+1]
                    if prot[b]=='$' and prot[b+1]=='M':
                        if min_pep_len<=i-(b+1)<=max_pep_len:
                            peptide = prot[b+2:i+1]
                            if td==True:
                                if reverseseq ==True:
                                    d_protnum = (currentmax-1)-protnum
                                else:
                                    d_protnum = protnum
                                add_info = fasta_ID+'|'+'d_'+str(d_protnum)+'|'+'-'+prot[i+1].replace('$','-')
                            else:
                                add_info = fasta_ID+'|'+str(protnum)+'|'+'-'+prot[i+1].replace('$','-')
                            CopyPeptideS1(peptide, add_info, M_ranges, 
                                            M, C, perbin_bytes, output_directory, TMT_labeling)
                            
                    if miscleavage ==True:                    
                        # when b!=a, miscleavage possible
                        if b!=a:
                            # Miscleavage peptide @ prot[a+1:i+1]
                            if min_pep_len<=i-a<=max_pep_len:
                                peptide = prot[a+1:i+1]
                                if td==True:
                                    if reverseseq ==True:
                                        d_protnum = (currentmax-1)-protnum
                                    else:
                                        d_protnum = protnum
                                    add_info = fasta_ID+'|'+'d_'+str(d_protnum)+'|'+(prot[a]+prot[i+1]).replace('$','-')
                                else:
                                    add_info = fasta_ID+'|'+str(protnum)+'|'+(prot[a]+prot[i+1]).replace('$','-')
                                CopyPeptideS1(peptide, add_info, M_ranges, 
                                            M, C, perbin_bytes, output_directory, TMT_labeling)
        
                            # N-term Methionine cleaved Miscleavage peptide @ prot[a+2:i+1]
                            if prot[a] == '$' and prot[a+1]=='M':
                                if min_pep_len<=i-a-1<=max_pep_len:
                                    peptide = prot[a+2:i+1]
                                    if td==True:
                                        if reverseseq ==True:
                                            d_protnum = (currentmax-1)-protnum
                                        else:
                                            d_protnum = protnum
                                        add_info = fasta_ID+'|'+'d_'+str(d_protnum)+'|'+'-'+prot[i+1].replace('$','-')
                                    else:
                                        add_info = fasta_ID+'|'+str(protnum)+'|'+'-'+prot[i+1].replace('$','-')
                                    CopyPeptideS1(peptide, add_info, M_ranges, 
                                            M, C, perbin_bytes, output_directory, TMT_labeling)
    
                #if nonstandardaa not assigned for original peptide (due to length)
                if nonstandardaa == 'not assigned':
                    if any(aa in prot[b+1:i+1] for aa in 'XBJOUZ?-*'):
                        nonstandardaa = True
                        
                if nonstandardaa == True:
                    a=i
                    if prot[i+1]=='$':
                        a+=1
                    b=a
                    continue
    
                # update a and b
                a = b
                b = i
                if prot[i+1]=='$':
                    a = i+1
                    b = a

    if stage == 'S2':
        protnum = -1
        f_idx = 0
        fasta_ID = fasta_protcount[f_idx][0]
        currentmax = fasta_protcount[f_idx][1]
        for i in xrange(len(prot)-1):
            if prot[i]=='$':
                protnum+=1
                if protnum>=currentmax:
                    f_idx+=1
                    fasta_ID = fasta_protcount[f_idx][0]
                    currentmax = fasta_protcount[f_idx][1]
                    protnum = 0
            if (prot[i] in ['K','R'] and prot[i+1]!='P') or prot[i+1]=='$':
                
                if miscleavage == False:
                    peptide = prot[b+1:i+1]
                    if td==True:
                        if reverseseq ==True:
                            d_protnum = (currentmax-1)-protnum
                        else:
                            d_protnum = protnum
                        add_info = fasta_ID+'|'+'d_'+str(d_protnum)+'|'+(prot[b]+prot[i+1]).replace('$','-')
                    else:
                        add_info = fasta_ID+'|'+str(protnum)+'|'+(prot[b]+prot[i+1]).replace('$','-')
                    CopyPeptideS2(peptide, add_info, M_ranges, M, C, perbin_bytes, 
                                    output_directory, min_pep_len, max_pep_len, True, True, TMT_labeling)   
    
                if miscleavage == True:
                    peptide = prot[a+1:i+1]
                    if td==True:
                        if reverseseq ==True:
                            d_protnum = (currentmax-1)-protnum
                        else:
                            d_protnum = protnum
                        add_info = fasta_ID+'|'+'d_'+str(d_protnum)+'|'+(prot[a]+prot[i+1]).replace('$','-')
                    else:
                        add_info = fasta_ID+'|'+str(protnum)+'|'+(prot[a]+prot[i+1]).replace('$','-')
                    CopyPeptideS2(peptide, add_info, M_ranges, M, C, perbin_bytes, 
                                    output_directory, min_pep_len, max_pep_len, True, True, TMT_labeling)
                    
                    if prot[a]=='$':
                        peptide = prot[a+1:b+1]
                        if td==True:
                            add_info = fasta_ID+'|'+'d_'+str(d_protnum)+'|'+(prot[a]+prot[b+1]).replace('$','-')
                        else:
                            add_info = fasta_ID+'|'+str(protnum)+'|'+(prot[a]+prot[b+1]).replace('$','-')         
                        CopyPeptideS2(peptide, add_info, M_ranges, M, C, perbin_bytes, 
                                    output_directory, min_pep_len, max_pep_len, False, True, TMT_labeling)                                           
                    
                    if prot[i+1]=='$':
                        peptide = prot[b+1:i+1]
                        if td==True:
                            add_info = fasta_ID+'|'+'d_'+str(d_protnum)+'|'+(prot[b]+prot[i+1]).replace('$','-')
                        else:
                            add_info = fasta_ID+'|'+str(protnum)+'|'+(prot[b]+prot[i+1]).replace('$','-')
                        CopyPeptideS2(peptide, add_info, M_ranges, M, C, perbin_bytes, 
                                    output_directory, min_pep_len, max_pep_len, True, False, TMT_labeling)
            
                a = b
                b = i
                if prot[i+1]=='$':
                    a = i+1
                    b = a

    for bin_i in xrange(massbin_n):
        if M[bin_i]:
            outfile = open(path.join(output_directory, 'DBPart_'+str(bin_i)+'.txt'),'a')
            outfile.write(''.join(M[bin_i]))
            outfile.close()
    return 

def MakePeptides(output_directory, fastadir, B1size, B2size, 
                 min_pep_len, max_pep_len, miscleavage, 
                 massbin_n, M_ranges, fastaidxmap, td, reverseseq, TMT_labeling, stage):

    begintime = time.time()
    maxbytes = B2size/massbin_n
    mfasta_sizes = 0
    mfasta_files = []
    for fastafile in listdir(fastadir):
        mfasta_sizes+=stat(path.join(fastadir, fastafile)).st_size 
        mfasta_files.append(fastafile)
        if mfasta_sizes>B1size:
            print '\treading ', len(mfasta_files), ' fasta files at ', mfasta_sizes/(1024*1024.0), ' bytes.'
            proteins = []
            fasta_protcount = []
            for fastafile in mfasta_files:
                with open(path.join(fastadir, fastafile),'r') as infile:
                    fasta_ID = fastaidxmap[path.basename(fastafile).replace('.fasta','')]
                    fasta_sequences = ''.join(['$' if line[0]=='>' else line.strip() for line in infile.readlines()])
                    fasta_protcount.append((fasta_ID, fasta_sequences.count('$')))
                    if td ==True:
                        if reverseseq ==True:
                            proteins.append('$'+fasta_sequences[::-1][:-1])
                        else:
                            proteins.append(fasta_sequences)
                    else:
                        proteins.append(fasta_sequences)
            
            proteins = ''.join(proteins)+'$'
            
            TrypticDigest(proteins,fasta_protcount,
                                min_pep_len, max_pep_len, 
                                miscleavage, M_ranges,
                                output_directory, massbin_n, 
                                maxbytes, td, reverseseq, stage, TMT_labeling)

            mfasta_sizes = 0
            mfasta_files = []

    # remaining fasta files
    if mfasta_files:
        print '\t#reading ', len(mfasta_files), ' fasta files at ', mfasta_sizes/(1024*1024.0), ' bytes.'
        proteins = []
        fasta_protcount = []
        for fastafile in mfasta_files:
            with open(path.join(fastadir, fastafile),'r') as infile:
                fasta_ID = fastaidxmap[path.basename(fastafile).replace('.fasta','')]
                fasta_sequences = ''.join(['$' if line[0]=='>' else line.strip() for line in infile.readlines()])                
                fasta_protcount.append((fasta_ID, fasta_sequences.count('$')))
                if td ==True:
                    if reverseseq ==True:
                        proteins.append('$'+fasta_sequences[::-1][:-1])
                    else:
                        proteins.append(fasta_sequences)
                else:
                    proteins.append(fasta_sequences)
        
        proteins = ''.join(proteins)+'$'
        
        TrypticDigest(proteins,fasta_protcount,
                            min_pep_len, max_pep_len, 
                            miscleavage, M_ranges,
                            output_directory, massbin_n, 
                            maxbytes, td, reverseseq, stage, TMT_labeling)

    return time.time()-begintime

def CreateDBPartitions(tempdir, finaldir, stage):
    begin_time= time.time()
    
    for f in listdir(tempdir):
        f1 = path.join(tempdir, f)
        f2 = path.join(finaldir, f)
        
        peptides = {}
        pepmasses =  {}
        with open(f1,'r') as infile:
            for line in infile:
                sp = line.strip().split('\t')
                pep = sp[0]            
                if pep in peptides:
                    peptides[pep].add(sp[2])
                else:
                    peptides[pep] = set([sp[2]])
                    pepmasses[pep] = float(sp[1])
    
        pepmasses = sorted(pepmasses.items(), key = itemgetter(1))

        with open(f2,'w') as outfile:
            for index, entry in enumerate(pepmasses):
                mass = str(entry[1])
                pepseq = entry[0]
                mappings = peptides[pepseq]
                prepost = [v.split('|')[2] for v in mappings]
                mappings_keep = []
                td = 0
                
                if stage== 'S1':
                    tryptic_target_idx = [v for v in mappings if 'd_' not in v]
                    if tryptic_target_idx:
                        pepseq_keep = tryptic_target_idx[0]
                        mappings_keep = tryptic_target_idx
                    else:
                        td = 1
                        tryptic_decoy_idx = set(mappings)-set(tryptic_target_idx)
                        mappings_keep = tryptic_decoy_idx
                        pepseq_keep = list(tryptic_decoy_idx)[0]  
    
                    prepost = pepseq_keep.split('|')[2]
        
                    if td ==1:
                        outfile.write(mass+'\t'+prepost[0]+'.'+pepseq+'.'+prepost[1]+'\t'+'XXX_'+str(index)+'\t'+';'.join(mappings_keep)+'\n')
                    if td==0:
                        outfile.write(mass+'\t'+prepost[0]+'.'+pepseq+'.'+prepost[1]+'\t'+str(index)+'\t'+';'.join(mappings_keep)+'\n')
    
                if stage== 'S2':
                    tryptic_idx = [zi for zi,z in enumerate(prepost) \
                        if (z[0] in ['K','R','-'] and (pepseq[-1] in ['K','R'] or z[1]=='-'))]
                
                    if tryptic_idx:
                        tryptic_target_idx = [v for vi,v in enumerate(mappings) \
                        if vi in set(tryptic_idx) and 'd_' not in v]
                            
                        if tryptic_target_idx:
                            pepseq_keep = tryptic_target_idx[0]
                        else:
                            td = 1
                            tryptic_decoy_idx = [v for vi, v in enumerate(mappings) \
                            if vi in set(tryptic_idx) and 'd_' in v]
                            pepseq_keep = tryptic_decoy_idx[0]
                            
                    elif not tryptic_idx:
                        semi_target_idx = [v for v in mappings if 'd_' not in v]
                        if semi_target_idx:
                            pepseq_keep = semi_target_idx[0]
                        else:
                            td=1
                            semi_decoy_idx = [v for v in mappings if 'd_' in v]
                            pepseq_keep = semi_decoy_idx[0]
                    
                    prepost = pepseq_keep.split('|')[2]
                    
                    if td ==1:
                        outfile.write(mass+'\t'+prepost[0]+'.'+pepseq+'.'+prepost[1]+'\t'+'XXX_'+str(index)+'\n')
                    if td==0:
                        outfile.write(mass+'\t'+prepost[0]+'.'+pepseq+'.'+prepost[1]+'\t'+str(index)+'\n')

    shutil.rmtree(tempdir)
    return time.time()-begin_time