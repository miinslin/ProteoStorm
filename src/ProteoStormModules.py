# -*- coding: utf-8 -*-
import time, os, subprocess
from CoreModule1_CreateDBPartitions import CreateDBPartitions
from CoreModule1_CreateSpectralPartitions import CreateSpectralParts
from CoreModule2and3_PepFilter_pValueComputation import pepfilter, pvaluecomputation

def RunModules(stage, mods, filter_cutoff,
               ProteoStorm_dir,subdir,
               spectral_dir, spectralparts_dir,
               spectra_remove, DBpart_massrange_file,
               miscleavages, min_pep_len, max_pep_len,
               REMOVE_KRP, max_massdiff,
               Pepfilter_exe, MSGFpvaluejar,
               aafreq_fasta, refinedDB_pepfdr,
               generaDB, PSlogfile, enzyme,
               precursormasstol, instrument, fragmentmet,
               RAMgb, parallel_n, cygwinpath, num_Spectra, fragmasstol):
    
    #log runtime
    RuntimeLog = {'spectralpartition':0,
                  'databasepartition':0,
                  'refineddb_s2input':0,
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
                                          spectra_remove, DBpart_massrange_file,\
                                          num_Spectra, PSlogfile, precursormasstol)
            RuntimeLog['spectralpartition'] = sp_time
        
        # Core module 1 - database partitioning
        PS_dir = os.path.join(ProteoStorm_dir, \
                              'S1_PreprocessingOutput', 'ProteoStorm_input')
        
        if os.path.exists(PS_dir) and os.listdir(PS_dir)!=[]:
            print '\tS1 database partitions already exist...skipping.'
            RuntimeLog['databasepartition'] = 'S1 database partitioning skipped. Partitions already existed.'

        else:
            dbp_time = CreateDBPartitions(ProteoStorm_dir, subdir,
                                     miscleavages,
                                     min_pep_len, max_pep_len, 
                                     REMOVE_KRP, max_massdiff,
                                     stage, enzyme, RAMgb,
                                     parallel_n, cygwinpath)
            RuntimeLog['databasepartition'] = dbp_time
            
    if stage == 'S2':
        print 'Beginning Core Module 1...'
        # Core module 1 - database creation + partitioning
        if generaDB ==1:
            from CreateRefinedProtDB_genus import CreateRefinedDBGenus
            rdb_s2i_time = CreateRefinedDBGenus(ProteoStorm_dir, subdir,
                         max_massdiff,
                         min_pep_len, max_pep_len,
                         miscleavages, refinedDB_pepfdr,
                         enzyme, PSlogfile,
                         RAMgb, parallel_n, cygwinpath)
            RuntimeLog['refineddb_s2input'] = rdb_s2i_time
            
        else:
            from CreateRefinedProtDB_peptide import CreateRefinedDBPeptide
            rdb_s2i_time = CreateRefinedDBPeptide(ProteoStorm_dir, subdir,
                         max_massdiff,
                         min_pep_len, max_pep_len,
                         miscleavages, refinedDB_pepfdr,
                         enzyme, PSlogfile,
                         RAMgb, parallel_n, cygwinpath)
            RuntimeLog['refineddb_s2input'] = rdb_s2i_time            

    # Core module 2 - PEPTIDE SPECTRUM PAIR FILTERING
    print 'Beginning Core Module 2...'
    pepfilter_cmds = pepfilter(Pepfilter_exe, ProteoStorm_dir,
                             subdir, spectralparts_dir,
                             stage, mods, filter_cutoff,
                             precursormasstol, fragmasstol)
    
    num_pepfilter_logfiles = pepfilter_cmds['num_commands']
    print '\texpecting ', num_pepfilter_logfiles, ' peptide filtering processes...'
    
    # commands are sorted by file size
    pepfilter_cmds = pepfilter_cmds['shcommand']
    # parallelize
    pepfilter_start_time = time.time()
    
    for i in range(0, parallel_n):
        task_idx = range(i,num_pepfilter_logfiles,parallel_n)
        print '\tTask',i+1,':', len(task_idx),' processes'
        pepfilter_sh = os.path.join(os.path.join(ProteoStorm_dir, subdir), stage+'_pepfilter_commands'+mods+'_task_'+str(i+1)+'.sh')
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
                    comment = "Peptide Filtering Processes Completed: "+str(len(pass_files))+\
                                    '\t'+time.strftime("%H:%M:%S")
                    PSlogfile.write(comment+'\n')
                    print '\t',comment
                    pass_write += int(num_pepfilter_logfiles/4) 
                time.sleep(30)
        except IOError:
            time.sleep(10)
            print 'wait...'
    
    RuntimeLog['coremodule2'] = pepfilter_end_time
    
    # P-VALUE COMPUTATION FOR (P,S) PAIRS
    print 'Beginning Core Module 3...'
    compute_pvalue_cmds = pvaluecomputation(MSGFpvaluejar, ProteoStorm_dir, \
                                          subdir, stage, spectralparts_dir,\
                                          aafreq_fasta, precursormasstol,\
                                          instrument, fragmentmet, RAMgb)

    num_pvalue_logfiles = compute_pvalue_cmds['num_commands']
    print '\texpecting ', num_pvalue_logfiles, ' pvalue computation processes...'
    
    # commands are sorted by file size
    compute_pvalue_cmds = compute_pvalue_cmds['shcommand']
    # parallelize
    pvalue_start_time = time.time()    

    for i in range(0, parallel_n):
        task_idx = range(i,num_pvalue_logfiles,parallel_n)
        print '\tTask',i+1,':', len(task_idx),' processes'
        compute_pvalue_sh = os.path.join(os.path.join(ProteoStorm_dir, subdir), stage+'_PVAL_computation_cmd'+'_task_'+str(i+1)+'.sh')
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
                        if 'MS-GF+ complete' in lines:
                            pass_files.add(log_fn)
            if len(pass_files)==num_pvalue_logfiles:
                pvalue_end_time = time.time()-pvalue_start_time
                break
            else:
                if len(pass_files)>=pass_write:
                    comment = "pValue Computation Processes Completed: "+str(len(pass_files))+\
                                '\t'+time.strftime("%H:%M:%S")
                    PSlogfile.write(comment+'\n')
                    print '\t',comment
                    pass_write += int(num_pvalue_logfiles/4)
                time.sleep(30)
        except IOError:
            time.sleep(10)
            print 'wait...'    
    
    RuntimeLog['coremodule3'] = pvalue_end_time
    
    return RuntimeLog
