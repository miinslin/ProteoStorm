# -*- coding: utf-8 -*-
import os, time, argparse
import platform
import numpy as np
from CalculatePeptideMass import calcmass_cmm
from Preprocessing import PanMicrobial_preprocessing, Genera_preprocessing
from ProteoStormModules import RunModules
from PSM_output import S2_PSMs
from S2_peplevelFDR import PeptideLvlFDRout

parser = argparse.ArgumentParser()
parser.add_argument("-O", "--output", dest = "output", \
                    help="Main directory for ProteoStorm output.")
parser.add_argument("-D", "--Database", dest = "db", \
                    help="Directory of protein database files in fasta format.")
parser.add_argument("-S", "--Spectra", dest = "spectra", \
                    help="Directory of spectral datasets in MGF format (converted from RAW format using MSConvert).")
parser.add_argument("-Sname", "--SpectralDataset", dest = "spectralds", \
                    help="Name of metaproteomics dataset. Will be used as subdirectory name.")
parser.add_argument("-r", "--RemoveSpectra",dest ="removespectra", default = 'na', \
                    help="File containing spectra filename and scan numbers to remove. See HS_matched_spectra.txt in example files.")
parser.add_argument("-pepfilter", "--PepfilterEXE", dest = "pepfilterexe", \
                    help="Path to core module 2 (peptide filtering) executable. \
                    Use pre-compiled binaries (CoreModule2_PeptideFiltering_Windows-x86_64.exe, CoreModule2_PeptideFiltering_Linux-x86_64.exe) \
                    or compile from source.")
parser.add_argument("-masswin", "--PartitionMassWindow", dest = "pmw", default = '15', \
                    help="Mass range of a partition in Daltons.")
parser.add_argument("-S1spc", "--S1SharedPeakCount", dest = "s1spc", default = '7', \
                    help="Number of shared peak counts to use for S1 peptide filtering.")
parser.add_argument("-S2spc", "--S2SharedPeakCount", dest = "s2spc", default = '6', \
                    help="Number of shared peak counts to use for S2 peptide filtering.")
parser.add_argument("-RefSeq", "--RefSeqCatalog", dest = "refcat", default = 'na',\
                    help = "Path to RefSeq Catalog.")
parser.add_argument("-ms1t", "--PrecursorMassTolerance", dest = "ms1pmt", default = '10',\
                    help = "Precursor mass tolerance in ppm (e.g., 10 for HCD).")
parser.add_argument("-ms2t", "--FragmentMassTolerance", dest = "ms2mt", default = '0.015',\
                    help = "Fragment mass tolerance in Daltons (e.g., 0.015 for HCD).")
parser.add_argument("-inst", "--InstrumentID", dest = "ms2dID", default = '3',\
                    help = "Identifier of the instrument to generate MS/MS spectra. Used to determine the scoring model in core module 3 (p-value computation). \
                    0: Low-res LCQ/LTQ, 1: Orbitrap/FTICR, 3: Q-Exactive (Default)")
parser.add_argument("-m", "--FragmentMethodID", dest = "fragmethod", default = '3',\
                    help = "Fragmentation method identifier. Used to determine the scoring model in core module 3 (p-value computation). \
                    1: CID, 3: HCD (Default)")
parser.add_argument("-mem ", "--RAMusage_GB", dest = "memory", default = '8',\
                    help = "Under development. Do not recommend changing.")
parser.add_argument("-genera", "--GeneraRestrictionApproach", dest = "gapp", default = '0',\
                    help = "0: Do not use genera-restriction approach. \
                    1: Use genera-restriction approach. Place UniProt fasta files in /[Database_dir]/UniProt_fasta,\
                    and RefSeq fasta files in /[Database_dir]/RefSeq_fasta.")
parser.add_argument("-p", "--parallel", dest = "para_n", default = '1',\
                    help = "Under development. Do not recommend changing.")
parser.add_argument("-pvaljar", "--pval_computation_jar", dest = "pvaljar", default = 'MSGFPlus_pvalue.jar',\
                    help = "Path to p-value computation jar.")
parser.add_argument("-aafreq", "--aminoacid_freq", dest = "aafreq", default = '364106_IL_transformed.fasta',\
                    help = "Path to file specifying amino acid frequencies to use for core module 3. Do not recommend changing.")
parser.add_argument("-cygwin", "--CygwinPATH", dest = "cygwinpath", default = '',\
                    help = "Path to cygwin for Windows users.")
args = parser.parse_args()

# PARAMETERS
masswindow = int(args.pmw)
S1_SPC = int(args.s1spc)
S2_SPC = int(args.s2spc)
parallel_n = int(args.para_n)
miscleavages = 1  # max number of miscleavages allowed.
min_pep_len = 8  # minimum number of amino acids in peptide
max_pep_len = 40  # maximum number of amino acids in peptide
REFINEDDB_pepFDR = 0.05 # peptide-level FDR for refined protein database
GENERA_RDB = int(args.gapp)
REMOVE_KRP = 1 # do not include semi-tryptic peptides with K.P or R.P for S2 input files
enzyme = 'trypsin' # does not cut after K/R if followed by P
precursormasstol = float(args.ms1pmt)
fragmasstol = float(args.ms2mt)
instrument = str(int(args.ms2dID))
fragmentmet = str(int(args.fragmethod))
RAMgb = str(int(args.memory))
MBsize = 60*int(RAMgb)/8 # split large pan-microbial database into chunks of size MBsize
num_Spectra = 200000*int(RAMgb)/8 # number of spectra to partition per iteration

platform_os = platform.platform() # determine platform
if platform_os.split('-')[0]=='Windows':
    cygwinpath = str(args.cygwinpath)
    if not os.path.exists(cygwinpath):
        raise ValueError(cygwinpath, ' does not exist. Should be path to run.exe for Cygwin.')
    cygwinpath = os.path.normpath(cygwinpath)
    
if platform_os.split('-')[0]=='Linux':
    cygwinpath = ''
    
if platform_os.split('-')[0] not in ['Linux','Windows']:
    raise ValueError('platform not supported!')

##### DEFINE DIRECTORIES AND FILENAMES
# INPUT FILES
FASTA_dir = os.path.normpath(args.db)
spectral_dir = os.path.normpath(args.spectra)
spectra_remove = str(args.removespectra)
if spectra_remove != 'na':
    if not os.path.exists(spectra_remove):
        raise ValueError('File does not exist: ',spectra_remove,'\nCheck usage using --help')
    spectra_remove = os.path.normpath(spectra_remove)
RefSeq_catalog = str(args.refcat)
if RefSeq_catalog != 'na':
    if not os.path.exists(RefSeq_catalog):
        raise ValueError('File does not exist: ',RefSeq_catalog,'\nCheck usage using --help')
    RefSeq_catalog = os.path.normpath(RefSeq_catalog)
    
# CORE MODULES 2 and 3 EXECUTABLES
Pepfilter_exe = os.path.normpath(args.pepfilterexe)
Pepfilter_exe = os.path.normpath(Pepfilter_exe).split(os.sep)
Pepfilter_exe = '/'.join(Pepfilter_exe)
MSGFpvaluejar = os.path.normpath(args.pvaljar)
aafreq_fasta = os.path.normpath(args.aafreq)

# PROTEOSTORM OUTPUT DIRECTORIES
ProteoStorm_dir = os.path.normpath(args.output)
subdir = str(args.spectralds)
S1_preprocessing_dir = os.path.join(ProteoStorm_dir, 'S1_PreprocessingOutput')
fastachunks_dir = os.path.join(S1_preprocessing_dir, 'FastaChunks')
main_dir = os.path.join(ProteoStorm_dir, subdir)
spectralparts_dir = os.path.join(main_dir, 'S1_InputFiles', 'SpectralParts')

# MASS RANGES FOR PARTITIONS
DBpart_massrange_file = os.path.join(S1_preprocessing_dir,
                                     'PartRanges_mc_' + str(miscleavages) +\
                                     '_md' + str(masswindow) + '.txt')

# CREATE DIRECTORIES
for directory in [fastachunks_dir, spectralparts_dir]:
    if not os.path.exists(directory):
        os.makedirs(directory)

#check if files exist
for filename in [FASTA_dir,spectral_dir, Pepfilter_exe,
                 MSGFpvaluejar, aafreq_fasta]:
    if not os.path.exists(filename):
        raise ValueError('File does not exist:', filename)

inst_type = {'0': 'Low-res LCQ/LTQ', '1': 'Orbitrap/FTICR', '3': 'Q-Exactive'}
fragment_type = {'1': 'CID', '3': 'HCD'}

# PROTEOSTORM LOGFILE
ProteoStormLOG = open(os.path.join(ProteoStorm_dir,subdir,'ProteoStorm.log'),'w')
ProteoStormLOG.write('ProteoStorm'+'\t'+time.strftime("%m/%d/%Y")+'\t'+time.strftime("%H:%M:%S")+'\n')
ProteoStormLOG.write('Platform: '+ platform_os+'\n')
ProteoStormLOG.write('Cygwin path: '+cygwinpath+'\n')
ProteoStormLOG.write('#params'+'\n'+\
                     'fasta chunk size ='+str(MBsize)+'\n'+\
                     'number of spectra to partition per iteration ='+str(num_Spectra)+'\n'+\
                     'max miscleavages ='+str(miscleavages)+'\n'+\
                     'partition mass window ='+str(masswindow)+'\n'+\
                     'minimum peptide length='+str(min_pep_len)+'\n'+\
                     'maximum peptide length='+str(max_pep_len)+'\n'+\
                     'enzyme='+str(enzyme)+'\n'+\
                     'precursor mass tolerance (ppm) ='+str(precursormasstol)+'\n'+\
                     'fragment mass tolerance (Da)'+str(fragmasstol)+'\n'+\
                     'instrument type ='+inst_type[instrument]+'\n'+\
                     'fragment method ='+fragment_type[fragmentmet]+'\n'+\
                     'refined DB peptide level FDR='+str(REFINEDDB_pepFDR)+'\n'+\
                     'S1 peptide filter shared counts='+str(S1_SPC)+'\n'+\
                     'S2 peptide filter shared counts='+str(S2_SPC)+'\n'+\
                     'Genera restriction approach ='+str(GENERA_RDB)+'\n')

##########################################################################################
Total_start_time = time.time()

####======= PREPROCESSING
if not os.path.exists(DBpart_massrange_file):
    print 'Creating partition mass range file...'
    # Create mass-range file
    minpepmass = int(np.floor(calcmass_cmm(min_pep_len * 'G', TMT_mod = 0)))
    maxpepmass = int(np.ceil(calcmass_cmm(max_pep_len * 'W', TMT_mod = 0)))
    set_masslist = range(minpepmass, maxpepmass + masswindow, masswindow)
    num_partitions = len(set_masslist)
    
    with open(DBpart_massrange_file, 'w') as outfile:
        outfile.write(
            '\n'.join(['\t'.join([str(x) for x in set_masslist[i:i + 2]])\
                     for i in range(num_partitions)\
                     if i != num_partitions - 1]))

if os.listdir(fastachunks_dir) == []:
    if GENERA_RDB ==0:
        print 'Splitting fasta files into chunks of size', MBsize, 'MB...'   
        ppfs_time= PanMicrobial_preprocessing(FASTA_dir,MBsize,fastachunks_dir)
        ProteoStormLOG.write('Preprocessing (split fasta): '+str(ppfs_time)+'\n')
    
    if GENERA_RDB ==1:
        print 'Retaining entries with genus...splitting fasta files into chunks of size', MBsize, 'MB...'
        ppfs_time = Genera_preprocessing(FASTA_dir, S1_preprocessing_dir, \
                                fastachunks_dir, MBsize,\
                                RefSeq_catalog, ProteoStormLOG)
        ProteoStormLOG.write('Preprocessing (identify genera + split fasta): '+str(ppfs_time)+'\n')

####======= STAGE ONE
print 'Beginning stage one...'
STAGEONE_start = time.time()
stage = 'S1'
mods = ''
runtimes = RunModules(stage, mods, S1_SPC,
               ProteoStorm_dir,subdir,
               spectral_dir, spectralparts_dir,
               spectra_remove, DBpart_massrange_file,
               miscleavages, min_pep_len, max_pep_len,
               REMOVE_KRP, masswindow,
               Pepfilter_exe, MSGFpvaluejar,
               aafreq_fasta, REFINEDDB_pepFDR,
               GENERA_RDB, ProteoStormLOG, enzyme,
               precursormasstol, instrument, fragmentmet,
               RAMgb, parallel_n, cygwinpath, num_Spectra, fragmasstol)

STAGEONE_end = time.time()-STAGEONE_start

ProteoStormLOG.write('Stage One Total:'+str(STAGEONE_end)+'\n'\
                     +'\t'+'CoreModule 1 (spectral partitioning):'+str(runtimes['spectralpartition'])+'\n'\
                     +'\t'+'CoreModule 1 (database partitioning):'+str(runtimes['databasepartition'])+'\n'\
                     +'\t'+'CoreModule 2 (peptide-spectrum pair filtering):'+str(runtimes['coremodule2'])+'\n'\
                     +'\t'+'CoreModule 3 (pvalue computation):'+str(runtimes['coremodule3'])+'\n')

print 'Stage one completed...\n'

####======= STAGE TWO
print 'Beginning stage two...'
STAGETWO_start = time.time()
stage = 'S2'
mods = ''
runtimes = RunModules(stage, mods, S2_SPC,
               ProteoStorm_dir,subdir,
               spectral_dir, spectralparts_dir,
               spectra_remove, DBpart_massrange_file,
               miscleavages, min_pep_len, max_pep_len,
               REMOVE_KRP, masswindow,
               Pepfilter_exe, MSGFpvaluejar,
               aafreq_fasta, REFINEDDB_pepFDR,
               GENERA_RDB, ProteoStormLOG, enzyme,
               precursormasstol, instrument, fragmentmet,
               RAMgb, parallel_n, cygwinpath, num_Spectra, fragmasstol)

STAGETWO_end = time.time()-STAGETWO_start

ProteoStormLOG.write('Stage Two Total:'+str(STAGETWO_end)+'\n'\
                     +'\t'+'CoreModule 1 (database partitioning):'+str(runtimes['refineddb_s2input'])+'\n'\
                     +'\t'+'CoreModule 2 (peptide-spectrum pair filtering):'+str(runtimes['coremodule2'])+'\n'\
                     +'\t'+'CoreModule 3 (pvalue computation):'+str(runtimes['coremodule3'])+'\n')

print 'Stage one completed...\n'

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