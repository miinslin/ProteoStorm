# -*- coding: utf-8 -*-
import os, re, time, subprocess
import shutil
from Bio import SeqIO
import numpy as np
from InSilicoDigest import insilicodigest
from CalculatePeptideMass import calcmass_cmm

def CreateInputFiles(fasta_thmmass_file, stage, \
                     Proteostorm_dir, subdir):
       
    S1_preprocessing_dir = os.path.join(Proteostorm_dir, 'S1_PreprocessingOutput')    
    
    if stage == 'S1':    
        protein_mappingsdir = os.path.join(S1_preprocessing_dir, 'ProteinMappings')
        PS_dir = os.path.join(S1_preprocessing_dir, 'ProteoStorm_input')
        numsorted_deduplicated_combined = os.path.join(S1_preprocessing_dir, \
                                                       'temp_DBfiles', 'NumSorted_combined_dd.txt')
        if not os.path.exists(protein_mappingsdir):
            os.makedirs(protein_mappingsdir)
    
    if stage == 'S2':
        PS_dir = os.path.join(Proteostorm_dir, subdir, 'S2_InputFiles', 'ProteoStorm_input')
        numsorted_deduplicated_combined = os.path.join(Proteostorm_dir, subdir, \
                                                       'S2_InputFiles', 'temp_DBfiles', 'NumSorted_combined_dd.txt')

    for newd in [PS_dir]:
        if not os.path.exists(newd):
            os.makedirs(newd)

    with open(fasta_thmmass_file,'r') as fastapartmass:
        dbp_ranges = np.array([float(z.split('\t')[1]) for z in fastapartmass.readlines()])
    
    with open(numsorted_deduplicated_combined,'r') as infile:
        for index, line in enumerate(infile): 
            if line.strip() == '100000000000000000000000000000000':
                break
            sp = line.strip().split('\t')
            thmass = float(sp[0])
            pepseq = sp[1]
            
            if index==0:
                DBidx = np.searchsorted(dbp_ranges,thmass,side='left')
                current_max = dbp_ranges[DBidx]
                if stage =='S1':
                    mappings_f = open(os.path.join(protein_mappingsdir,'DBPart_'+str(DBidx)+'.mapping'),'w')
                PS_f = open(os.path.join(PS_dir,'DBPart_'+str(DBidx)+'.txt'),'w')

            if thmass >= current_max:
                if stage =='S1':
                    mappings_f.close()
                PS_f.close()               
                
                DBidx = np.searchsorted(dbp_ranges,thmass,side='left')
                current_max = dbp_ranges[DBidx]
                
                if stage =='S1':
                    mappings_f = open(os.path.join(protein_mappingsdir,'DBPart_'+str(DBidx)+'.mapping'),'w')
                PS_f = open(os.path.join(PS_dir,'DBPart_'+str(DBidx)+'.txt'),'w')

                
            if thmass<current_max:
                # chunk|startidx/ 'd_' or ''|pre+post
                mappings = sp[2].split(';')
                prepost = [v.split('|')[2] for v in mappings]
                mappings_keep = []
                
                #and z[1]!='P'
                td = 0
                tryptic_idx = [zi for zi,z in enumerate(prepost) \
                if (z[0] in ['K','R','-'] and (pepseq[-1] in ['K','R'] or z[1]=='-'))]

                # if fully-tryptic peptide exists
                if tryptic_idx:
                    # take target if possible
                    tryptic_target_idx = [v for vi,v in enumerate(mappings) \
                    if vi in set(tryptic_idx) and 'd_' not in v]
                        
                    if tryptic_target_idx:
                        pepseq_keep = tryptic_target_idx[0]
                        mappings_keep = tryptic_target_idx
                    else:
                        # elif no target, take any tryptic decoy
                        td = 1
                        tryptic_decoy_idx = [v for vi, v in enumerate(mappings) \
                        if vi in set(tryptic_idx) and 'd_' in v]
                        pepseq_keep = tryptic_decoy_idx[0]
                        mappings_keep = tryptic_decoy_idx
                        
                #if no fully-tryptic peptide exists
                elif not tryptic_idx:
                    # all semi-tryptic
                    # target
                    semi_target_idx = [v for v in mappings if 'd_' not in v]
                    if semi_target_idx:
                        pepseq_keep = semi_target_idx[0]
                        mappings_keep = semi_target_idx
                    else:
                        td=1
                        semi_decoy_idx = [v for v in mappings if 'd_' in v]
                        pepseq_keep = semi_decoy_idx[0]
                        mappings_keep = semi_decoy_idx
                    
                prepost = pepseq_keep.split('|')[2]
                
                if td ==1:
                    PS_f.write(sp[0]+'\t'+prepost[0]+'.'+pepseq+'.'+prepost[1]+'\t'+'XXX_'+str(index)+'\n')
                if td==0:
                    PS_f.write(sp[0]+'\t'+prepost[0]+'.'+pepseq+'.'+prepost[1]+'\t'+str(index)+'\n')

                if stage =='S1':
                    if not mappings_keep:
                        raise ValueError('Error')
                    mappings_f.write(pepseq+'\t'+';'.join(mappings_keep)+'\n')

        if stage =='S1':
            mappings_f.close()
        PS_f.close()   

    return '\tCreated database partitions...'

def CreateDBPartitions(Proteostorm_dir, subdir, miscleavages, \
                       min_pep_len, max_pep_len, REMOVE_KRP,\
                       massdiff, stage, enzyme, RAMgb, parallel_n, cygwinpath):
    
    print 'Database partitioning...'
    
    start_time = time.time()
    
    S1_preprocessing_dir = os.path.join(Proteostorm_dir, 'S1_PreprocessingOutput')
    fasta_thmmass_file = os.path.join(S1_preprocessing_dir, 'PartRanges_mc_1_md'+str(massdiff)+'.txt')
    S1_fastachunk_dir = os.path.join(S1_preprocessing_dir, 'FastaChunks')
    enzymerules = {'trypsin':r'([KR*](?=[^P]))'}
    regex_enz = re.compile(enzymerules[enzyme])
    
    if stage == 'S1':
        tempfiledir = os.path.join(S1_preprocessing_dir, 'temp_DBfiles')
        
    if stage == 'S2':
        tempfiledir = os.path.join(Proteostorm_dir, subdir, 'S2_InputFiles', 'temp_DBfiles')    

    if not os.path.exists(tempfiledir):
        os.makedirs(tempfiledir)

    combinedfname = os.path.join(tempfiledir,'Unsorted_combined.txt')

    if stage == 'S1' and not os.path.exists(combinedfname):
        
        with open(combinedfname,'w') as outfile:
            outfile.write('ZZZZEND'+'\n')        
        
        for f in [chunk for chunk in os.listdir(S1_fastachunk_dir) if chunk.endswith('.fasta') and 'revCat.fasta' not in chunk]:            
            Fastafile = os.path.join(S1_fastachunk_dir,f)
            chunkID = Fastafile[:-6].split('_')[-1]
                
            # $ indicates start, * indicates end of protein
            allprots = ''.join(['$'+str(record.seq)+'*' for record in SeqIO.parse(Fastafile, "fasta")])
            # result of insilico digest is {peptide key: (start index or ''/'d_', pre+post)}
            result = insilicodigest(allprots, regex_enz, min_pep_len, max_pep_len, miscleavages,'full','', REMOVE_KRP)
            del allprots
            allprots_decoy = ''.join(['$'+str(record.seq)[::-1]+'*' for record in SeqIO.parse(Fastafile, "fasta")])
            result_d = insilicodigest(allprots_decoy, regex_enz, min_pep_len, max_pep_len, miscleavages,'full','d_', REMOVE_KRP)
            del allprots_decoy
        
            #write to file
            if os.path.exists(combinedfname):
                outfile = open(combinedfname,'a')
            elif not os.path.exists(combinedfname):
                outfile = open(combinedfname,'w')
        
            for pep in result:
                outfile.write(pep+'\t'+';'.join([chunkID+'|'+k[0]+'|'+k[1] for k in result[pep]]))
                if pep in result_d:
                    outfile.write(';'+';'.join([chunkID+'|'+k[0]+'|'+k[1] for k in result_d[pep]]))
                outfile.write('\n')
            for pep in result_d:
                if pep not in result:
                    outfile.write(pep+'\t'+';'.join([chunkID+'|'+k[0]+'|'+k[1] for k in result_d[pep]])+'\n')
        
            outfile.close()
            del result
            del result_d
            
        print '\tFinished creating unsorted file A...'

    ####
    # sort combined file alphabetically    
    sortedcf = os.path.join(tempfiledir,'Sorted_combined.txt') 
    
    if not os.path.exists(sortedcf):
        if cygwinpath:
            cmd = [cygwinpath,
                   'sort',
                   '-S'+str(RAMgb)+'G',
                   '--parallel='+str(parallel_n),
                   '-T',tempfiledir,
                   combinedfname,
                   '-o', sortedcf]
        else:
            cmd = ['sort',
                   '-S'+str(RAMgb)+'G',
                   '--parallel='+str(parallel_n),
                   '-T',tempfiledir,
                   combinedfname,
                   '-o', sortedcf]
        
        subprocess.Popen(cmd, stdout = subprocess.PIPE)
        
        #check for 'ZZZZEND' at end
        stop = 0
        while True:
            if os.path.exists(sortedcf):
                with open(sortedcf, 'r') as infile:
                    for line in infile:
                        if line.strip()=='ZZZZEND':
                            print '\t Finished alphabetically sorting file...'
                            stop=1
                if stop ==1:
                    break
                elif stop ==0:
                    time.sleep(20)
            
            elif not os.path.exists(sortedcf):
                time.sleep(5)
    
        #remove unsorted file
        os.remove(combinedfname)
    
    ###    
    # deduplicate sorted file and calculate mass of peptides
    deduplicatedf = os.path.join(tempfiledir,'Sorted_combined_dd.txt')

    if not os.path.exists(deduplicatedf):
        def deduplicateTemp(sortedcf, deduplicatedf):    
            lastpep = ''
            lastmap = set()
            with open(sortedcf,'r') as infile, open(deduplicatedf,'w') as outfile:
                # write 100000000000000000000000000000000 to file
                outfile.write('100000000000000000000000000000000'+'\n')
                
                for i,line in enumerate(infile):
                    #if last line ZZZZEND
                    if line.strip() == 'ZZZZEND':
                        outfile.write(';'.join(lastmap)+'\n')
                        return '\t Finished deduplicating file...'
                    
                    sp = line.strip().split('\t')
                    currentpep = sp[0]
                    currentmap = set(sp[1].split(';'))
        
                    if currentpep==lastpep:
                        lastmap.update(currentmap)
                        
                    elif currentpep!= lastpep:
                        if lastmap:
                            outfile.write(';'.join(lastmap)+'\n')
                        thmass = str(calcmass_cmm(currentpep, TMT_mod = 0))
                        outfile.write(thmass+'\t'+currentpep+'\t')
                        lastpep = str(currentpep)
                        lastmap = set(currentmap)
                        
        print deduplicateTemp(sortedcf, deduplicatedf)
        # remove sorted combined file
        os.remove(sortedcf)

    ### 
    # sort numerically deduplicated file
    numsort_deduplicatedf = os.path.join(tempfiledir,'NumSorted_combined_dd.txt')

    if not os.path.exists(numsort_deduplicatedf):
        if cygwinpath:
            cmd = [cygwinpath,
                   'sort',
                   '-n',
                   '-S'+str(RAMgb)+'G',
                   '--parallel='+str(parallel_n),
                   '-T',tempfiledir,
                   deduplicatedf,
                   '-o', numsort_deduplicatedf]
        else:
            cmd = ['sort',
                   '-n',
                   '-S'+str(RAMgb)+'G',
                   '--parallel='+str(parallel_n),
                   '-T',tempfiledir,
                   deduplicatedf,
                   '-o', numsort_deduplicatedf]
        
        subprocess.Popen(cmd, stdout = subprocess.PIPE)
    
        #check for '100000000000000000000000000000000' at end
        stop = 0
        while True:
            if os.path.exists(numsort_deduplicatedf):
                with open(numsort_deduplicatedf, 'r') as infile:
                    for line in infile:
                        if line.strip()=='100000000000000000000000000000000':
                            print '\t Finished sorting numerically...'
                            stop=1
                if stop ==1:
                    break
                elif stop ==0:
                    time.sleep(20)
            
            elif not os.path.exists(numsort_deduplicatedf):
                time.sleep(5)

        #remove deduplicated file
        os.remove(deduplicatedf)
    
    # create input files
    print CreateInputFiles(fasta_thmmass_file, stage, Proteostorm_dir, subdir)
    
    # delete temporary directory 
    shutil.rmtree(tempfiledir)

    return time.time()-start_time

#def S2_Partitions_Mods(Proteostorm_dir, subdir):
#    start_time = time.time()
#    
#    maindir = os.path.join(Proteostorm_dir, subdir, 'S2_InputFiles')
#    fasta_thmmass_file = os.path.join(Proteostorm_dir, 'S1_PreprocessingOutput','PartRanges_mc_1_md15.txt')
#    msgf_modificatoins_file = r'D:\MSGF_top4mods.txt'
#    
#    PS_dir = os.path.join(maindir, 'ProteoStorm_input')
#    PS_mods_dir = os.path.join(maindir, 'ProteoStorm_input_mods')
#    
#    if not os.path.exists(PS_mods_dir):
#        os.makedirs(PS_mods_dir)
#    
##==============================================================================
##     mods_mass = {'Carbamidomethyl':57.021464,
##                  'Oxidation':15.994915,
##                  'Deamidated':0.984016,
##                  'Carbamyl':43.005814,
##                  'Glu->pyro-Glu':-18.010565,
##                  'Gln->pyro-Glu':-17.026549,
##                  'Acetyl':42.010565,
##                  'Methyl':14.015650,
##                  'Phospho':79.966331,
##                  'Formyl':27.994915
##                  }
##==============================================================================
#    
#    #msgfplus/src/main/java/edu/ucsd/msjava/msutil/Composition.java
#    element_mass = {'C':12.0,
#                    'H':1.007825035,
#                    'N':14.003074,
#                    'O':15.99491463,
#                    'S':31.9720707,
#                    'P':30.973762,
#                    'Br':78.9183361,
#                    'Cl':34.96885272,
#                    'Fe':55.9349393,
#                    'Se':79.9165196
#                    }
#    
#    
#    with open(fasta_thmmass_file,'r') as fastapartmass:
#        dbp_ranges = np.array([float(z.split('\t')[1]) for z in fastapartmass.readlines()])
#    
#    modslist = {}
#    #set modifications
#    with open(msgf_modificatoins_file, 'r') as infile:
#        modid = -1
#        for line in infile:
#            if line[0]!='#' and line.strip()!='':
#                if 'NumMods' in line.strip():
#                    NumMods = line.strip().split('NumMods=')[1]
#                #Mass or CompositionStr, Residues, ModType, Position, Name (all the five fields are required)
#                else:
#                    item = line.strip().split()[0].split(',')
#                    # CompositionStr (C[Num]H[Num]N[Num]O[Num]S[Num]P[Num]Br[Num]Cl[Num]Fe[Num])
#                    mass_composition = item[0]
#                    if any(aa in mass_composition for aa in ['C','H','N','O','S','P','Br','Cl','Fe']):
#                        mass_composition = [x for x in re.split(r'(-?\d+)', mass_composition) if x]
#                        mass_composition = [(mass_composition[i],mass_composition[i+1]) for i in range(0,len(mass_composition),2)]
#                        mass_composition = round(sum([element_mass[element[0]]*int(element[1]) for element in mass_composition]),6)
#                    else:
#                        mass_composition = round(float(mass_composition),6)
#                    
#                    # * for any residue. Can't be * and anywhere mod                        
#                    residues = item[1]
#                    if '*' in residues and len(residues)!=1:
#                        raise ValueError('residues cannot be * and amino acid...')
#                    if residues!='*' and not all(aa in 'ACEDGFIHKMLNQPSRTWVY' for aa in residues):
#                        raise ValueError('one or more residue not standard amino acid...')
#                    # fix or opt for variable                    
#                    modtype = item[2]
#                    #skip fixed mods
#                    if modtype == 'fix':
#                        continue
#                    # any (anywhere), N-term (peptide N-term), C-term (peptide C-term), 
#                    #   Prot-N-term (protein N-term), Prot-C-term (protein C-term)
#                    position = item[3]
#                    if position == 'any' and residues == '*':
#                        raise ValueError('position cannot be any (anywhere) if residues is *...')
#                    if position not in ['any', 'N-term', 'C-term', 'Prot-N-term', 'Prot-C-term']:
#                        raise ValueError('not acceptable position...')
#    #                        # Unimod PSI-MS name
#    #                        PSI_MSname = item[4]
#                    modid+=1
#                    # add to mods dictionary
#                    modslist[modid] = {'mass':mass_composition,'residues':residues, 'modtype':modtype, 'position':position}
#    
#    
#    
#    def findmodifiedpeptides(modslist, peptide, pepmass, preaa, postaa):
#        candidate_modified_peplist_all = {}
#        
#        for mod in modslist:  
#            position = modslist[mod]['position']
#            residues = modslist[mod]['residues']
#            mass = modslist[mod]['mass']
#    
#            modified_mass = pepmass+mass
#            
#            prefix = ''
#            if -1.0*mass<0:
#                prefix = '+'
#            
#            candidate_modified_peplist = []
#            # if modification can be on any residue at position...
#            if residues == '*':
#                if position == 'N-term':
#                    candidate_modified_peplist = [prefix+str(mass)+peptide]
#                elif position == 'C-term':
#                    candidate_modified_peplist = [peptide+prefix+str(mass)]
#                elif position == 'Prot-N-term':
#                    if preaa == '-':
#                        candidate_modified_peplist = [prefix+str(mass)+peptide]
#                elif position == 'Prot-C-term':
#                    if postaa == '-':
#                        candidate_modified_peplist = [peptide+prefix+str(mass)]
#            # if specific residues to be modified...
#            elif residues != '*':
#                if position == 'any':
#                    candidate_modified_peplist = [peptide[:ai+1]+prefix+str(mass)+peptide[ai+1:] for ai, aa in enumerate(peptide) if aa in residues]
#                elif position == 'N-term':
#                    if peptide[0] in residues:
#                        candidate_modified_peplist = [prefix+str(mass)+peptide]
#                elif position == 'C-term':
#                    if peptide[-1] in residues:
#                        candidate_modified_peplist = [peptide+prefix+str(mass)]
#                elif position == 'Prot-N-term':
#                    if preaa == '-' and peptide[0] in residues:
#                        candidate_modified_peplist = [prefix+str(mass)+peptide]
#                elif position == 'Prot-C-term':
#                    if postaa == '-' and peptide[-1] in residues:
#                        candidate_modified_peplist = [peptide+prefix+str(mass)]
#            
#            if candidate_modified_peplist:
#                if modified_mass not in candidate_modified_peplist_all:
#                    candidate_modified_peplist_all[modified_mass] = set() 
#                candidate_modified_peplist_all[modified_mass].update(candidate_modified_peplist)
#            
#        return candidate_modified_peplist_all
#    
#    ####################################################################################################
#    
#    print modslist
#    
#    targetmodindx = -1
#    decoymodindx = -1
#    
#    for DB_partition in range(len(dbp_ranges)):
#        PS_f = os.path.join(PS_dir, 'DBPart_'+str(DB_partition)+'.txt')
#        if os.path.exists(PS_f):
#            mods_toadd = {x:{} for x in range(len(dbp_ranges))}
#            
#            with open(PS_f,'r') as infile:
#                for line in infile:
#                    sp = line.strip().split('\t')
#                    #474.18227445	R.GGGGGGGG.R	0
#                    preaa, postaa = sp[1][0],sp[1][-1]
#                    peptide = sp[1][2:-2]
#                    td = 'target'
#                    if 'XXX_' in sp[2]:
#                        td = 'decoy'
#                    pepmass = float(sp[0])
#                    pepmodlist = findmodifiedpeptides(modslist, peptide, pepmass, preaa, postaa)
#                    for modmass in pepmodlist:
#                        # which partition to add to
#                        DBidx = np.searchsorted(dbp_ranges, modmass ,side='left')
#                        if modmass not in mods_toadd[DBidx]:
#                            mods_toadd[DBidx][modmass] = {'target':set(),'decoy':set()}
#                        mods_toadd[DBidx][modmass][td].update([preaa+'.'+cp+'.'+postaa for cp in pepmodlist[modmass]])
#            
#            # for each partition, write modified peptide to mod part.
#            for DBidx in mods_toadd:
#                if mods_toadd[DBidx]:
#                    if os.path.exists(os.path.join(PS_mods_dir,'DBPart_'+str(DBidx)+'_modpeps.txt')):
#                        PS_f_mods = open(os.path.join(PS_mods_dir,'DBPart_'+str(DBidx)+'_modpeps.txt'),'a')                    
#                    
#                    if not os.path.exists(os.path.join(PS_mods_dir,'DBPart_'+str(DBidx)+'_modpeps.txt')):
#                        PS_f_mods = open(os.path.join(PS_mods_dir,'DBPart_'+str(DBidx)+'_modpeps.txt'),'w')
#                    
#                    for modmass in mods_toadd[DBidx]:
#                        for mpepseq in mods_toadd[DBidx][modmass]['target']:
#                            targetmodindx+=1
#                            header = 'mod'+str(targetmodindx)
#                            PS_f_mods.write(str(modmass)+'\t'+mpepseq+'\t'+header+'\n')
#                            
#                        for mpepseq in mods_toadd[DBidx][modmass]['decoy']:
#                            decoymodindx+=1
#                            header = 'XXX_mod'+str(decoymodindx)
#                            PS_f_mods.write(str(modmass)+'\t'+mpepseq+'\t'+header+'\n')
#        
#    
#    #### sort
#    from shutil import copyfile
#    from operator import itemgetter
#    
#    for DB_partition in range(len(dbp_ranges)):
#        PS_f = os.path.join(PS_dir, 'DBPart_'+str(DB_partition)+'.txt')
#        PS_f_mods = os.path.join(PS_mods_dir,'DBPart_'+str(DB_partition)+'_modpeps.txt')
#        PS_f_mods_numsorted = os.path.join(PS_mods_dir,'Top4Mods_DBPart_'+str(DB_partition)+'.txt')
#        
#        if os.path.exists(PS_f): 
#            if os.path.exists(PS_f_mods):
#                with open(PS_f,'r') as a, open(PS_f_mods,'r') as b,\
#                open(PS_f_mods_numsorted,'w') as outfile:
#                    writelines = [x.strip() for x in a.readlines()]
#                    writelines.extend([x.strip() for x in b.readlines()])
#                    masses = sorted([(float(x.split('\t')[0]), xi) for xi,x in enumerate(writelines)], key = itemgetter(0))
#                    masses = [x[1] for x in masses]
#                    outfile.write('\n'.join([writelines[x] for x in masses]))
#                
#                #delete PS_f_mods
#                os.remove(PS_f_mods)
#        
#            elif not os.path.exists(PS_f_mods):
#                copyfile(PS_f, PS_f_mods_numsorted)
#    
#        else:
#            if os.path.exists(PS_f_mods):
#                os.rename(PS_f_mods, PS_f_mods_numsorted)
#    
#    return time.time()-start_time
