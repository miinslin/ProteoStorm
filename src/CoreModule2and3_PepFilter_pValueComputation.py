# -*- coding: utf-8 -*-
import os

def pepfilter(proteostorm_exe, Proteostorm_dir, subdir, 
              spectraparts_dir, stage, mods, raw_score_cutoff,
              precursormasstol, fragmasstol, TMT_labeling):
    
    maindir = os.path.join(Proteostorm_dir, subdir)
    
    if stage =='S1':
        proteostorminput_dir = os.path.join(Proteostorm_dir, 'S1_PreprocessingOutput', 'ProteoStorm_input')
    
    if stage =='S2':
        proteostorminput_dir = os.path.join(maindir, 'S2_InputFiles', 'ProteoStorm_input'+mods)
    
    proteostormoutput_dir = os.path.join(maindir, stage+'_OutputFiles', 'ProteoStorm_filtered'+mods)
    proteostormlogfile_dir = os.path.join(maindir, stage+'_OutputFiles', 'ProteoStorm_filtering'+mods+'_logs')
    
    for newd in [proteostormoutput_dir, proteostormlogfile_dir]:
        if not os.path.exists(newd):
            os.makedirs(newd)
    
    number_of_commands = 0
    commandsoutput = []
    
    PS_inputfiles = os.listdir(proteostorminput_dir)
    # sort by filesize
    PS_inputfiles = [(x, os.path.getsize(os.path.join(proteostorminput_dir, x))) for x in PS_inputfiles]
    PS_inputfiles.sort(key=lambda s: s[1])
    PS_inputfiles = [x[0] for x in PS_inputfiles]
    
    tmtparam = ' -tmt 0'
    
    if TMT_labeling ==1:
        tmtparam = ' -tmt 1'
    
    for filenum, ps_input in enumerate(PS_inputfiles):
        indx = ps_input[:-4].split('_')[-1]
        specfiles = [s for s in os.listdir(spectraparts_dir) if 'SpecPart_'+indx+'_' in s]
        
        if not specfiles:
            continue
        
        number_of_commands+=len(specfiles)
        
        for s in specfiles:
            specfile = os.path.join(spectraparts_dir, s)
            logfilename = os.path.join(proteostormlogfile_dir, s[:-4]+'.log')
            outfilename = os.path.join(proteostormoutput_dir, s[:-4]+'.tsv')

            command = proteostorm_exe+' -i '+ '"'+specfile+'" -p "'+os.path.join(proteostorminput_dir, ps_input)+'"'+\
            ' -o "'+outfilename+'"'+\
            ' -ms1tol '+str(int(precursormasstol))+' -ms2tol '+str(fragmasstol)+tmtparam+\
            ' -s '+str(raw_score_cutoff)+' > "'+ logfilename+'"'
        
            commandsoutput.append(command)

    return {'shcommand':commandsoutput,'num_commands':number_of_commands}

def pvaluecomputation(MSGFpvaluejar, Proteostorm_dir, subdir,
                      stage, spectraparts_dir, aafreq_fasta,
                      precursormasstol, instrument, fragmentmet, RAMgb,
                      min_pep_len, max_pep_len, modsfile):
    
    num_threads = 1
    maindir = os.path.join(Proteostorm_dir, subdir)

    outputdir = os.path.join(maindir,stage+'_OutputFiles')
    ProteoStorm_filtered = os.path.join(outputdir, 'ProteoStorm_filtered')
    logdir = os.path.join(outputdir, 'PVAL_computation_logs')
    pval_output = os.path.join(outputdir,'PVAL_computations')
    ntt = '1' # 1 semi, 2 full
    number_of_commands = 0
    
    for newd in [pval_output, logdir]:
        if not os.path.exists(newd):
            os.makedirs(newd)
    
    pepfilter_files = os.listdir(ProteoStorm_filtered)
    # sort by filesize
    pepfilter_files = [(x, os.path.getsize(os.path.join(ProteoStorm_filtered, x))) for x in pepfilter_files]
    pepfilter_files.sort(key=lambda s: s[1])
    pepfilter_files = [x[0] for x in pepfilter_files]
    
    commandsoutput = []
    for ps_filter in pepfilter_files:
        #check if anything in output
        with open(os.path.join(ProteoStorm_filtered, ps_filter),'r') as infile:
            entries = len([x for x in infile.readlines() if x[0]!='#'])
            if entries ==0:
                continue
            
        specfile = os.path.join(spectraparts_dir, ps_filter[:-4]+'.mgf')
        logfilename = os.path.join(logdir, ps_filter[:-4]+'.log')
        outfilename = os.path.join(pval_output, ps_filter[:-4]+'.tsv')
        
        number_of_commands+=1
        command = 'java -Xmx'+RAMgb+'G -jar "'+MSGFpvaluejar+'" -s '+\
        '"'+specfile+'" -d "'+aafreq_fasta+'"'+' -x "'+os.path.join(ProteoStorm_filtered, ps_filter)+'"'+\
        ' -o "'+outfilename+'"'+\
        ' -t '+str(precursormasstol)+'ppm -tda 0 -m '+str(fragmentmet)+' -inst '+str(instrument)+\
        ' -minLength '+str(int(min_pep_len))+' -maxLength '+str(int(max_pep_len))+\
        ' -mod "'+modsfile+'"'+\
        ' -e 1 -ntt '+ntt+' -thread '+str(num_threads)+' > "'+logfilename+'"'
        
        commandsoutput.append(command)            
    
    return {'shcommand':commandsoutput,'num_commands':number_of_commands}
