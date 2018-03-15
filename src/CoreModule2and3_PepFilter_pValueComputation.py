# -*- coding: utf-8 -*-
import os

def pepfilter(proteostorm_exe, Proteostorm_dir, subdir, 
              spectraparts_dir, stage, mods, raw_score_cutoff,
              precursormasstol, fragmasstol):
    
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
    
    for filenum, ps_input in enumerate(PS_inputfiles):
        indx = ps_input[:-4].split('_')[-1]
        specfile = os.path.join(spectraparts_dir, 'SpecPart_'+indx+'.mgf')
        if not os.path.exists(specfile):
            continue
        number_of_commands+=1
        logfilename = os.path.join(proteostormlogfile_dir,'part_'+indx+'.log')
        
        command = proteostorm_exe+' -i '+ '"'+specfile+'" -p "'+os.path.join(proteostorminput_dir, ps_input)+'"'+\
        ' -o "'+os.path.join(proteostormoutput_dir, 'part_'+indx+'.tsv')+'"'+\
        ' -ms1tol '+str(int(precursormasstol))+' -ms2tol '+str(fragmasstol)+\
        ' -s '+str(raw_score_cutoff)+' > "'+ logfilename+'"'
        
        commandsoutput.append(command)

    return {'shcommand':commandsoutput,'num_commands':number_of_commands}

def pvaluecomputation(MSGFpvaluejar, Proteostorm_dir, subdir,
                      stage, spectraparts_dir, aafreq_fasta,
                      precursormasstol, instrument, fragmentmet, RAMgb):
    
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
            
        indx = ps_filter[:-4].split('_')[-1]
        specfile = os.path.join(spectraparts_dir, 'SpecPart_'+indx+'.mgf')
        logfilename = os.path.join(logdir,'SpecEvalues_part_'+indx+'.log')
        
        number_of_commands+=1
        command = 'java -Xmx'+RAMgb+'G -jar "'+MSGFpvaluejar+'" -s '+\
        '"'+specfile+'" -d "'+aafreq_fasta+'"'+' -x "'+os.path.join(ProteoStorm_filtered, ps_filter)+'"'+\
        ' -o "'+os.path.join(pval_output, 'SpecEvalues_part_'+indx+'.tsv')+'"'+\
        ' -t '+str(precursormasstol)+'ppm -tda 0 -m '+str(fragmentmet)+' -inst '+str(instrument)+\
        ' -e 1 -ntt '+ntt+' -thread 1 > "'+logfilename+'"'
        
        commandsoutput.append(command)            
    
    return {'shcommand':commandsoutput,'num_commands':number_of_commands}
