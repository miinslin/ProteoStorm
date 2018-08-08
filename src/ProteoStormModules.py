# -*- coding: utf-8 -*-
import time,  subprocess, shutil
import os
from CoreModule1_CreateDBPartitions import CreateDBPartitions, MakePeptides
from CoreModule1_CreateSpectralPartitions import CreateSpectralParts
from CoreModule2and3_PepFilter_pValueComputation import pepfilter, pvaluecomputation

def RunModules(stage, mods, filter_cutoff,
               ProteoStorm_dir,subdir,
               FASTA_dir, fastaIDmap,
               bufsize, bufsize2,
               spectral_dir, spectralparts_dir,
               spectra_remove, massbinranges,
               miscleavages, min_pep_len, max_pep_len,
               REMOVE_KRP, massbins,
               Pepfilter_exe, MSGFpvaluejar,
               aafreq_fasta, refinedDB_pepfdr,
               generaDB, PSlogfile, enzyme,
               precursormasstol, instrument, fragmentmet,
               RAMgb, parallel_n, cygwinpath, num_Spectra, 
               fragmasstol, save_space, modsfile, TMT_labeling):
    
    #log runtime
    RuntimeLog = {'spectralpartition':0,
                  'dbpartition':0,
                  'RefinedDBcreation':0,
                  'coremodule2':0,
                  'coremodule3':0}
    
    if stage == 'S1':
        print 'Beginning Core Module 1...'
        # Core module 1 - spectral partitioning
        if os.listdir(spectralparts_dir)!=[]:
            print '\tSpectral partitions already exist...skipping.'
            RuntimeLog['spectralpartition'] = 'Spectral partitioning skipped. Partitions already exsisted.'
            
        elif os.listdir(spectralparts_dir)==[]:
            sp_time = CreateSpectralParts(spectral_dir, spectralparts_dir, \
                                          spectra_remove, massbinranges,\
                                          num_Spectra, PSlogfile, precursormasstol, massbins)
            RuntimeLog['spectralpartition'] = sp_time
        
        # Core module 1 - database partitioning
        S1_preprocessing_dir = os.path.join(ProteoStorm_dir, 'S1_PreprocessingOutput')
        PS_dir = os.path.join(S1_preprocessing_dir, 'ProteoStorm_input')
        
        if os.path.exists(PS_dir) and os.listdir(PS_dir)!=[]:
            print '\tS1 database partitions already exist...skipping.'
            RuntimeLog['dbpartition'] = 'S1 database partitioning skipped. Partitions already existed.'
        
        else:
            print '\tCreating S1 database partitions...'
            PS_pre = os.path.join(S1_preprocessing_dir, 'PRE_ProteoStorm_input')
            if os.path.isdir(PS_pre):
                shutil.rmtree(PS_pre)
            os.mkdir(PS_pre)
            os.mkdir(PS_dir)
            
            bin_ends = [x[1] for x in massbinranges]
            #target peptides (reverse seq, False)
            t1 = MakePeptides(PS_pre, FASTA_dir, bufsize, bufsize2,
                                   min_pep_len, max_pep_len, miscleavages,
                                   massbins, bin_ends, 
                                   fastaIDmap, False, False, TMT_labeling, stage)
            #decoy peptides (reverse seq, True)
            t2 = MakePeptides(PS_pre, FASTA_dir, bufsize, bufsize2,
                                   min_pep_len, max_pep_len, miscleavages,
                                   massbins, bin_ends, 
                                   fastaIDmap, True, True, TMT_labeling, stage)
            #create partitions
            t3 = CreateDBPartitions(PS_pre, PS_dir, stage)
            RuntimeLog['dbpartition'] = t1+t2+t3
         
     
    if stage == 'S2':
        print 'Creating refined protein database...'
        # Core module 1 - database creation + partitioning
        if generaDB ==1:
            from CreateRefinedProtDB_genus import CreateRefinedDBGenus
            rdb_s2i_time = CreateRefinedDBGenus(ProteoStorm_dir, subdir,
                         refinedDB_pepfdr, PSlogfile, fastaIDmap, FASTA_dir, save_space)
            RuntimeLog['RefinedDBcreation'] = rdb_s2i_time
            
        else:
            from CreateRefinedProtDB_peptide import CreateRefinedDBPeptide
            rdb_s2i_time = CreateRefinedDBPeptide(ProteoStorm_dir, subdir,
                         refinedDB_pepfdr, PSlogfile, fastaIDmap, FASTA_dir, save_space)
            RuntimeLog['RefinedDBcreation'] = rdb_s2i_time

        #remove S1 pval out and log files
        if save_space ==1:
            pvalue_output_dir = os.path.join(ProteoStorm_dir, subdir, 'S1_OutputFiles', 'PVAL_computations')
            pvalue_logfile_dir = os.path.join(ProteoStorm_dir, subdir, 'S1_OutputFiles', 'PVAL_computation_logs')
            shutil.rmtree(pvalue_output_dir)
            shutil.rmtree(pvalue_logfile_dir)

        print 'Beginning Core Module 1...'
        print '\tCreating S2 database partitions...'
        S2_InputFiles = os.path.join(ProteoStorm_dir, subdir,'S2_InputFiles')
        PS_dir = os.path.join(S2_InputFiles, 'ProteoStorm_input')
        PS_pre = os.path.join(S2_InputFiles, 'PRE_ProteoStorm_input')
        os.makedirs(PS_dir)
        os.makedirs(PS_pre)
        refined_prot_DB_t = os.path.join(ProteoStorm_dir, subdir, 'S1_OutputFiles', 'RefinedProteinDB_t')
        refined_prot_DB_d = os.path.join(ProteoStorm_dir, subdir, 'S1_OutputFiles', 'RefinedProteinDB_d')
        
        s2_fastaIDmap = {'RefinedProteinDB_target':'0','RefinedProteinDB_decoy':'1'}
        
        bin_ends = [x[1] for x in massbinranges]
        #target peptides (reverse seq, False)
        t1 = MakePeptides(PS_pre, refined_prot_DB_t, bufsize, bufsize2,
                               min_pep_len, max_pep_len, miscleavages,
                               massbins, bin_ends, 
                               s2_fastaIDmap, False, False, TMT_labeling, stage)
        #decoy peptides (reverse seq, False)
        t2 = MakePeptides(PS_pre, refined_prot_DB_d, bufsize, bufsize2,
                               min_pep_len, max_pep_len, miscleavages,
                               massbins, bin_ends, 
                               s2_fastaIDmap, True, False, TMT_labeling, stage)
        #create partitions
        t3 = CreateDBPartitions(PS_pre, PS_dir, stage)
        #delete target decoy refined db directories
        shutil.rmtree(refined_prot_DB_t)
        shutil.rmtree(refined_prot_DB_d)
        RuntimeLog['dbpartition'] = t1+t2+t3

    #####
    # Core module 2 - PEPTIDE SPECTRUM PAIR FILTERING
    print 'Beginning Core Module 2...'
    pepfilter_cmds = pepfilter(Pepfilter_exe, ProteoStorm_dir,
                             subdir, spectralparts_dir,
                             stage, mods, filter_cutoff,
                             precursormasstol, fragmasstol, TMT_labeling)
    
    num_pepfilter_logfiles = pepfilter_cmds['num_commands']
    print '\texpecting ', num_pepfilter_logfiles, ' peptide filtering processes...'
    
    # commands are sorted by file size
    pepfilter_cmds = pepfilter_cmds['shcommand']
    # parallelize
    pepfilter_start_time = time.time()
    # list of sh files
    pepfilter_sh_files = []
    
    for i in range(0, parallel_n):
        task_idx = range(i,num_pepfilter_logfiles,parallel_n)
        print '\tTask',i+1,':', len(task_idx),' processes'
        pepfilter_sh = os.path.join(os.path.join(ProteoStorm_dir, subdir), stage+'_pepfilter_commands'+mods+'_task_'+str(i+1)+'.sh')
        pepfilter_sh_files.append(pepfilter_sh)
        
        with open(pepfilter_sh,'w') as outfile:
            outfile.write('&&'.join([pepfilter_cmds[x] for x in task_idx]))
        
        if cygwinpath:
            run_pepfilter_sh = [cygwinpath,'bash', pepfilter_sh]
        else:
            run_pepfilter_sh = ['bash', pepfilter_sh]        

        subprocess.Popen(run_pepfilter_sh, stdout=subprocess.PIPE)
    
    del pepfilter_cmds
    # CHECK LOGFILES
    logfile_dir = os.path.join(ProteoStorm_dir, subdir,\
                               stage+'_OutputFiles', 'ProteoStorm_filtering'+mods+'_logs')
    
    pass_write = int(num_pepfilter_logfiles/4)
    pass_files = set()
    while True:
        try:
            for log_fn in os.listdir(logfile_dir):
                if log_fn not in pass_files:
                    with open(os.path.join(logfile_dir,log_fn),'r') as infile:
                        lines = ' '.join(infile.readlines())
                        if 'Core Module 1 (peptide filtering) completed in' in lines:
                            pass_files.add(log_fn)
            if len(pass_files)==num_pepfilter_logfiles:
                pepfilter_end_time = time.time()-pepfilter_start_time
                break
            else:
                if len(pass_files)>=pass_write:
                    comment = "Peptide Filtering Processes Completed: "+str(len(pass_files))
                    PSlogfile.write(comment+'\n')
                    print '\t',comment
                    pass_write += int(num_pepfilter_logfiles/4) 
                time.sleep(30)
        except IOError:
            time.sleep(10)
            print 'wait...'
    
    RuntimeLog['coremodule2'] = pepfilter_end_time
    
    if save_space ==1:
        # delete S2 Input files
        if stage == 'S2':
            S2_InputFiles = os.path.join(ProteoStorm_dir, subdir,'S2_InputFiles')
            shutil.rmtree(S2_InputFiles)
        #delete command files
        for sh_file in pepfilter_sh_files:
            if os.path.exists(sh_file):
                os.remove(sh_file)
        
    # P-VALUE COMPUTATION FOR (P,S) PAIRS
    print 'Beginning Core Module 3...'
    compute_pvalue_cmds = pvaluecomputation(MSGFpvaluejar, ProteoStorm_dir,
                                          subdir, stage, spectralparts_dir,
                                          aafreq_fasta, precursormasstol,
                                          instrument, fragmentmet, RAMgb,
                                          min_pep_len, max_pep_len, modsfile)

    num_pvalue_logfiles = compute_pvalue_cmds['num_commands']
    print '\texpecting ', num_pvalue_logfiles, ' pvalue computation processes...'
    
    # commands are sorted by file size
    compute_pvalue_cmds = compute_pvalue_cmds['shcommand']
    # parallelize
    pvalue_start_time = time.time()    
    # list of sh files
    pvalue_sh_files = []
    
    for i in range(0, parallel_n):
        task_idx = range(i,num_pvalue_logfiles,parallel_n)
        print '\tTask',i+1,':', len(task_idx),' processes'
        compute_pvalue_sh = os.path.join(os.path.join(ProteoStorm_dir, subdir), stage+'_PVAL_computation_cmd'+'_task_'+str(i+1)+'.sh')
        pvalue_sh_files.append(compute_pvalue_sh)
        
        with open(compute_pvalue_sh,'w') as outfile:
            outfile.write('&&'.join([compute_pvalue_cmds[x] for x in task_idx]))
        
        if cygwinpath:
            run_compute_pvalue_sh = [cygwinpath,'bash', compute_pvalue_sh]
        else:
            run_compute_pvalue_sh = ['bash', compute_pvalue_sh]        

        subprocess.Popen(run_compute_pvalue_sh, stdout=subprocess.PIPE)
    
    del compute_pvalue_cmds
    # CHECK LOGFILES
    logfile_dir = os.path.join(ProteoStorm_dir, subdir,\
                          stage+'_OutputFiles','PVAL_computation_logs')
    
    pass_write = int(num_pvalue_logfiles/4)
    pass_files = set()
    while True:
        try:
            for log_fn in os.listdir(logfile_dir):
                if log_fn not in pass_files:
                    with open(os.path.join(logfile_dir,log_fn),'r') as infile:
                        lines = ' '.join(infile.readlines())
                        if 'p-value computation for peptide-spectrum pairs complete' in lines:
                            pass_files.add(log_fn)
            if len(pass_files)==num_pvalue_logfiles:
                pvalue_end_time = time.time()-pvalue_start_time
                break
            else:
                if len(pass_files)>=pass_write:
                    comment = "pValue Computation Processes Completed: "+str(len(pass_files))
                    PSlogfile.write(comment+'\n')
                    print '\t',comment
                    pass_write += int(num_pvalue_logfiles/4)
                time.sleep(30)
        except IOError:
            time.sleep(10)
            print 'wait...'    
    
    RuntimeLog['coremodule3'] = pvalue_end_time
    
    if save_space ==1:
        #delete command files
        for sh_file in pvalue_sh_files:
            if os.path.exists(sh_file):
                os.remove(sh_file)
        # delete pepfilter out and log files
        pepfilter_output_dir = os.path.join(ProteoStorm_dir, subdir, stage+'_OutputFiles', 'ProteoStorm_filtered'+mods)
        pepfilter_logfile_dir = os.path.join(ProteoStorm_dir, subdir, stage+'_OutputFiles', 'ProteoStorm_filtering'+mods+'_logs')
        shutil.rmtree(pepfilter_output_dir)
        shutil.rmtree(pepfilter_logfile_dir)
    
    return RuntimeLog
