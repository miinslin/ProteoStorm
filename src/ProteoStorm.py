# -*- coding: utf-8 -*-
import os, time, argparse, sys, platform
from ProteoStormModules import RunModules
from PSM_output import S2_PSMs
from S2_peplevelFDR import PeptideLvlFDRout

## CHECK PLATFORM
platform_os = platform.system()
if platform_os not in ['Linux','Windows']:
    raise ValueError(platform_os, ' platform system not supported. Exiting.')

## PARSE ARGUMENTS
script_directory = os.path.abspath(os.path.dirname(sys.argv[0])) #src directory

Pepfilter_exe_OS = {'Linux':'CoreModule2_PeptideFiltering_Linux-x86_64.exe',\
                'Windows':'CoreModule2_PeptideFiltering_Windows-x86_64.exe'}

current_datetime = time.strftime("%b-%d-%Y_%H:%M:%S")
date_split = current_datetime.split('_')[0].split('-')
time_split = current_datetime.split('_')[1].split(':') 

parser = argparse.ArgumentParser()
# REQUIRED:
parser.add_argument("-O", "--output", dest = "output", \
                    help="Main directory for ProteoStorm output.")
parser.add_argument("-D", "--Database", dest = "db", \
                    help="Directory of protein database files in fasta format.")
parser.add_argument("-S", "--Spectra", dest = "spectra", \
                    help="Directory of spectral datasets in MGF format (must peak pick and convert from RAW format using MSConvert).")

parser.add_argument("-dist", "--PepMassDistribution", dest = "pmdist", default = os.path.join(script_directory,'DBmassDistributions','UniProt_1146874024.txt'),\
                    help="File containing peptide counts per theoretical mass.")

# REQUIRED for WINDOWS users:
parser.add_argument("-cygwin", "--CygwinPATH", dest = "cygwinpath", default = 'C:/cygwin64/bin/run.exe',\
                    help = "Path to cygwin for Windows users. Default = C:/cygwin64/bin/run.exe")

# OPTIONAL
parser.add_argument("-n", "--database_partitions", dest = "db_n", default = 200,\
                    help="Number of database partitions.")

parser.add_argument("-MSMS", "--SpectralDataset", dest = "spectralds", default = date_split[0]+date_split[1]+'_'+''.join(time_split), \
                    help="Name of metaproteomics dataset used as subdirectory name. Default: time date")
parser.add_argument("-r", "--RemoveSpectra",dest ="removespectra", default = 'na', \
                    help="File containing spectra filename and scan numbers to remove. See HS_matched_spectra.txt in example files. Default: 'na'")

parser.add_argument("-pepfilter", "--PepfilterEXE", dest = "pepfilterexe", default = os.path.join(script_directory,Pepfilter_exe_OS[platform_os]),\
                    help="Path to core module 2 (peptide filtering) executable. Default: determined by OS")
parser.add_argument("-pvaljar", "--pval_computation_jar", dest = "pvaljar", default = os.path.join(script_directory,'MSGFPlus_pvalue.jar'),\
                    help = "Path to p-value computation jar. Default: MSGFPlus_pvalue_withmods.jar")
parser.add_argument("-aafreq", "--aminoacid_freq", dest = "aafreq", default = os.path.join(script_directory,'364106_IL_transformed.fasta'),\
                    help = "Path to file specifying amino acid frequencies to use for core module 3. Do not recommend changing. Default: 364106_IL_transformed.fasta")

parser.add_argument("-ms1t", "--PrecursorMassTolerance", dest = "ms1pmt", default = '10',\
                    help = "Precursor mass tolerance in ppm (e.g., 10, 20) Default: 10")
parser.add_argument("-ms2t", "--FragmentMassTolerance", dest = "ms2mt", default = '0.015',\
                    help = "Fragment mass tolerance in Daltons (e.g., 0.015, 0.6). Default: 0.015")
parser.add_argument("-inst", "--InstrumentID", dest = "ms2dID", default = '3',\
                    help = "Identifier of the instrument to generate MS/MS spectra. Used to determine the scoring model in core module 3 (p-value computation). \
                    0: Low-res LCQ/LTQ, 1: Orbitrap/FTICR, 3: Q-Exactive (Default)")
parser.add_argument("-m", "--FragmentMethodID", dest = "fragmethod", default = '3',\
                    help = "Fragmentation method identifier. Used to determine the scoring model in core module 3 (p-value computation). \
                    1: CID, 3: HCD (Default)")
parser.add_argument("-S1spc", "--S1SharedPeakCount", dest = "s1spc", default = '7', \
                    help="Number of shared peak counts to use for S1 peptide filtering. Default: 7")
parser.add_argument("-S2spc", "--S2SharedPeakCount", dest = "s2spc", default = '7', \
                    help="Number of shared peak counts to use for S2 peptide filtering. Default: 7")

parser.add_argument("-genera", "--GeneraRestrictionApproach", dest = "gapp", default = '0',\
                    help = "0: Do not use genera-restriction approach (Default). \
                    1: Use genera-restriction approach. Place UniProt fasta files in /[Database_dir]/UniProt_fasta,\
                    and RefSeq fasta files in /[Database_dir]/RefSeq_fasta.")
parser.add_argument("-TMT", "--TMTlabeling", dest = "tmtlabel", default = '0',\
                    help = "Use TMT labeling. Default: 0")

# DO NOT CHANGE:
#parser.add_argument("-mem ", "--RAMusage_GB", dest = "memory", default = '8',\
#                    help = "Do not recommend changing. Default: 8GB")
parser.add_argument("-p", "--parallel", dest = "para_n", default = '1',\
                    help = "Do not recommend changing. Default: 1 thread")
parser.add_argument("-del", "--deletefiles", dest = "del_ss", default = '1',\
                    help = "0: Keep or 1: delete files for development purposes. Default: 1")
parser.add_argument("-s1fdr", "--refDBfdr", dest = "refinedDB_fdr",default = 0.05,\
                    help = "Do not recommend changing. Default: 0.05")

args = parser.parse_args()

## SET PARAMETERS
massbins = int(args.db_n)
S1_SPC = int(args.s1spc)
S2_SPC = int(args.s2spc)
parallel_n = int(args.para_n)
pepmass_dist_file = os.path.normpath(args.pmdist)
bufsize = 1000*1024*1024
bufsize2 = 600*1024*1024
miscleavages = 1  # max number of miscleavages allowed.
min_pep_len = 8  # minimum number of amino acids in peptide
max_pep_len = 40  # maximum number of amino acids in peptide
REFINEDDB_pepFDR = float(args.refinedDB_fdr) # peptide-level FDR for refined protein database
GENERA_RDB = int(args.gapp)

enzyme = 'trypsin' # does not cut after K/R if followed by P
precursormasstol = float(args.ms1pmt)
fragmasstol = float(args.ms2mt)
instrument = str(int(args.ms2dID))
fragmentmet = str(int(args.fragmethod))
TMT_labeling = int(args.tmtlabel)
RAMgb = '8' #str(int(args.memory))

num_Spectra = 100000*int(RAMgb)/8 # number of spectra to partition per iteration
save_space = int(args.del_ss)# save disk space
modsfile = os.path.normpath(os.path.join(script_directory,'MSGF_C57.txt'))

if TMT_labeling not in [0,1]:
    TMT_labeling = 0

if TMT_labeling == 1:
    modsfile = os.path.normpath(os.path.join(script_directory,'MSGF_C57_TMT.txt'))
    print 'using modfile: ', modsfile
    #pepmass_dist_file = os.path.normpath(os.path.join(script_directory, 'DBmassDistributions','RefUpTMT_2872778677.txt'))
    print 'using mass distribution file: ', pepmass_dist_file
    

## SET DIRECTORIES and FILES
FASTA_dir = os.path.normpath(args.db)
spectral_dir = os.path.normpath(args.spectra)
ProteoStorm_dir = os.path.normpath(args.output)
spectra_remove = os.path.normpath(args.removespectra)

subdir = str(args.spectralds)
Pepfilter_exe = os.path.normpath(args.pepfilterexe).split(os.sep)
Pepfilter_exe = '/'.join(Pepfilter_exe)
MSGFpvaluejar = os.path.normpath(args.pvaljar)
aafreq_fasta = os.path.normpath(args.aafreq)
main_dir = os.path.join(ProteoStorm_dir, subdir)
S1_preprocessing_dir = os.path.join(ProteoStorm_dir, 'S1_PreprocessingOutput')
spectralparts_dir = os.path.join(main_dir, 'S1_InputFiles', 'SpectralParts')

if platform_os=='Linux':
    cygwinpath = ''
    
DBpart_massrange_file = os.path.join(S1_preprocessing_dir,
                                     'PartRanges_mc_' + str(miscleavages) +\
                                     '_bins_' + str(massbins) + '.txt')

for directory in [spectralparts_dir, S1_preprocessing_dir]:
    if not os.path.exists(directory):
        os.makedirs(directory)

## CHECK DIRECTORIES and FILES
ErrorNum = 0
ErrorStr = {}
dne_dir_files = []

if platform_os=='Windows':
    cygwinpath = os.path.normpath(args.cygwinpath)
    if not os.path.isfile(cygwinpath):
        dne_dir_files.append(cygwinpath)
if spectra_remove != 'na':
    if not os.path.isfile(spectra_remove):
        dne_dir_files.append(spectra_remove)
for dirname in [FASTA_dir, spectral_dir, ProteoStorm_dir]:
    if not os.path.isdir(dirname):
        dne_dir_files.append(dirname)
for filename in [Pepfilter_exe, MSGFpvaluejar, aafreq_fasta, modsfile, pepmass_dist_file]:
    if not os.path.isfile(filename):
        dne_dir_files.append(filename)
if dne_dir_files:
    ErrorNum+=len(dne_dir_files)
    ErrorStr['Directory/File(s) dne'] = dne_dir_files

def checkSpecFileFormat(spectral_dir):
    spectra_error = 0
    for f in os.listdir(spectral_dir):
        if not f.endswith('.mgf'):
            spectra_error += 1  
    if spectra_error!=0:
        ErrorStr['Specfiles'] = str(spectra_error)+' spectrafiles not in .mgf format.'
        return 1

if ErrorStr:
    ErrorLog = os.path.join(script_directory,'Error.log')
    with open(ErrorLog,'w') as outfile:
        for error_key in ErrorStr:
            outfile.write(error_key+'\n'+'\n\t'.join(ErrorStr[error_key])+'\n')
    raise ValueError('Check Error.log')

## PROTEOSTORM LOGFILE
inst_type = {'0': 'Low-res LCQ/LTQ', '1': 'Orbitrap/FTICR', '3': 'Q-Exactive'}
fragment_type = {'1': 'CID', '3': 'HCD'}

ProteoStormLOG = open(os.path.join(ProteoStorm_dir,subdir,'ProteoStorm.log'),'w')
ProteoStormLOG.write('ProteoStorm'+'\t'+'-'.join(date_split)+'\t'+':'.join(time_split)+'\n')
ProteoStormLOG.write('Platform: '+ platform_os+'\n')
ProteoStormLOG.write('Cygwin path: '+cygwinpath+'\n')
ProteoStormLOG.write('#params'+'\n'+\
                     'Number of spectra to partition per iteration = '+str(num_Spectra)+'\n'+\
                     'Max miscleavages = '+str(miscleavages)+'\n'+\
                     'Number of database partitions = '+str(massbins)+'\n'+\
                     'Minimum peptide length = '+str(min_pep_len)+'\n'+\
                     'Maximum peptide length = '+str(max_pep_len)+'\n'+\
                     'Enzyme = '+str(enzyme)+'\n'+\
                     'Isotope Error <= '+'1'+'\n'+\
                     'Precursor mass tolerance (ppm) = '+str(precursormasstol)+'\n'+\
                     'Fragment mass tolerance (Da) = '+str(fragmasstol)+'\n'+\
                     'Instrument type = '+inst_type[instrument]+'\n'+\
                     'Fragment method = '+fragment_type[fragmentmet]+'\n'+\
                     'Refined DB peptide level FDR = '+str(REFINEDDB_pepFDR)+'\n'+\
                     'S1 peptide filter shared counts = '+str(S1_SPC)+'\n'+\
                     'S2 peptide filter shared counts = '+str(S2_SPC)+'\n'+\
                     'Genera restriction approach = '+str(GENERA_RDB)+'\n'+\
                     'Remove files to save space = '+str(save_space)+'\n'+\
                     'TMT labeling = '+str(TMT_labeling)+'\n'+\
                     'Core Module 2 (peptide filter exe) = '+Pepfilter_exe+'\n'+\
                     'Core Modeule 3 (pvalue computation) jar = '+MSGFpvaluejar+'\n'+\
                     'Modsfile = '+modsfile+'\n'+\
                     'Peptide mass distribution file = '+pepmass_dist_file+'\n')

miscleavages = True

###############################################################################
Total_start_time = time.time()
###############################################################################

####======= PREPROCESSING
def DefineMassRange(distributionfile, massbinnum):
    MassRangeList = []
    totalpeps = int(os.path.basename(distributionfile).split('_')[1].replace('.txt',''))
    entries_per_bin = int((totalpeps/massbinnum)+1)
    with open(distributionfile,'r') as infile:
        currentcount =0
        range_start = False
        for index,line in enumerate(infile):
            sp = line.strip().split('\t')
            pepmass = sp[0]
            if range_start == False:
                range_start = pepmass
            currentcount += int(sp[1])
            if currentcount>=entries_per_bin:
                MassRangeList.append((round(float(range_start),4),
                                      round(float(pepmass),4)))
                range_start = False
                currentcount = 0
        if currentcount<entries_per_bin:
            MassRangeList.append((round(float(range_start),4),
                                      round(float(pepmass),4)))
    return MassRangeList

print 'Defining partition mass ranges based on distribution file...'
if not os.path.exists(DBpart_massrange_file):
    massbinranges = DefineMassRange(pepmass_dist_file, massbins)
    with open(DBpart_massrange_file, 'w') as outfile:
        outfile.write('\n'.join(str(x[0])+'\t'+str(x[1]) for x in massbinranges))
else:
    print '\tmass range file already created...reading file...'
    massbinranges = []
    with open(DBpart_massrange_file, 'r') as infile:
        for line in infile:
            sp = line.strip().split('\t')
            massbinranges.append((float(sp[0]),float(sp[1])))

print 'Creating Fasta to Index mapping file...'
fasta_idx_mapping = os.path.join(S1_preprocessing_dir,'Fasta_index_mapping.txt')
if not os.path.isfile(fasta_idx_mapping):
    fastaIDmap = {}
    with open(fasta_idx_mapping,'w') as outfile:
        for fi, f in enumerate(os.listdir(FASTA_dir)):
            fastaIDmap[f.replace('.fasta','')] = str(fi)
            outfile.write(str(fi)+'\t'+f.replace('.fasta','')+'\n')
else:
    print '\t fasta index mapping file already created...reading file...'
    with open(fasta_idx_mapping,'r') as infile:
        fastaIDmap = {x.strip().split('\t')[1]:x.strip().split('\t')[0] for x in infile.readlines()}

   
###======= STAGE ONE
print 'Beginning stage one...'
STAGEONE_start = time.time()
stage = 'S1'
mods = ''
runtimes = RunModules(stage, mods, S1_SPC,
               ProteoStorm_dir,subdir,
               FASTA_dir, fastaIDmap,
               bufsize, bufsize2,
               spectral_dir, spectralparts_dir,
               spectra_remove, massbinranges,
               miscleavages, min_pep_len, max_pep_len,
               massbins,
               Pepfilter_exe, MSGFpvaluejar,
               aafreq_fasta, REFINEDDB_pepFDR,
               GENERA_RDB, ProteoStormLOG, enzyme,
               precursormasstol, instrument, fragmentmet,
               RAMgb, parallel_n, cygwinpath, num_Spectra, 
               fragmasstol, save_space, modsfile, TMT_labeling)

STAGEONE_end = time.time()-STAGEONE_start

ProteoStormLOG.write('Stage One Total:'+str(STAGEONE_end)+'\n'\
                     +'\t'+'CoreModule 1 (spectral partitioning):'+str(runtimes['spectralpartition'])+'\n'\
                     +'\t'+'CoreModule 1 (database partitioning):'+str(runtimes['dbpartition'])+'\n'\
                     +'\t'+'CoreModule 2 (peptide-spectrum pair filtering):'+str(runtimes['coremodule2'])+'\n'\
                     +'\t'+'CoreModule 3 (pvalue computation):'+str(runtimes['coremodule3'])+'\n')

print 'Stage one completed...\n'

####======= STAGE TWO ========#####
print 'Beginning stage two...'
STAGETWO_start = time.time()
stage = 'S2'
mods = ''
runtimes = RunModules(stage, mods, S2_SPC,
               ProteoStorm_dir,subdir,
               FASTA_dir, fastaIDmap,
               bufsize, bufsize2,
               spectral_dir, spectralparts_dir,
               spectra_remove, massbinranges,
               miscleavages, min_pep_len, max_pep_len,
               massbins,
               Pepfilter_exe, MSGFpvaluejar,
               aafreq_fasta, REFINEDDB_pepFDR,
               GENERA_RDB, ProteoStormLOG, enzyme,
               precursormasstol, instrument, fragmentmet,
               RAMgb, parallel_n, cygwinpath, num_Spectra, 
               fragmasstol, save_space, modsfile, TMT_labeling)

STAGETWO_end = time.time()-STAGETWO_start

ProteoStormLOG.write('Stage Two Total:'+str(STAGETWO_end)+'\n'\
                     +'\t'+'Refined Protein DB creation:'+str(runtimes['RefinedDBcreation'])+'\n'\
                     +'\t'+'CoreModule 1 (database partitioning):'+str(runtimes['dbpartition'])+'\n'\
                     +'\t'+'CoreModule 2 (peptide-spectrum pair filtering):'+str(runtimes['coremodule2'])+'\n'\
                     +'\t'+'CoreModule 3 (pvalue computation):'+str(runtimes['coremodule3'])+'\n')

print 'Stage two completed...\n'

####======= WRITE PROTEOSTORM OUTPUT
print '\nWriting PSMs to file...'
S2output_dir = os.path.join(ProteoStorm_dir, subdir, 'S2_OutputFiles')
S2_PSMs(S2output_dir, ProteoStormLOG, subdir)

Total_end_time = time.time()-Total_start_time
print 'ProteoStorm Complete:'+str(Total_end_time)+' seconds.'+'\n'
ProteoStormLOG.write('ProteoStorm Complete:'+str(Total_end_time)+'\n')
ProteoStormLOG.write(time.strftime("%H:%M:%S")+'\n')

####======= 1% peptide-level FDR
print 'Computing 1% peptide level FDR'
S2_peplvl_FDR = 0.01
PeptideLvlFDRout(S2output_dir, S2_peplvl_FDR, ProteoStormLOG)

ProteoStormLOG.close()

# delete pval output, logs, and S1 input files
if save_space ==1:
    import shutil
    pvalue_output_dir = os.path.join(ProteoStorm_dir, subdir, stage+'_OutputFiles', 'PVAL_computations')
    pvalue_logfile_dir = os.path.join(ProteoStorm_dir, subdir, stage+'_OutputFiles', 'PVAL_computation_logs')
    shutil.rmtree(pvalue_output_dir)
    shutil.rmtree(pvalue_logfile_dir)
    shutil.rmtree(os.path.join(main_dir, 'S1_InputFiles'))