# ProteoStorm v 1.0

Software Requirements
==========
1. ```Python2.7``` ([Anaconda 2.7](https://www.anaconda.com/download/?lang=en-us)):
	- numpy
	- biopython
2. ```Java```
3. ```Cygwin``` (required for Windows users)

**Linux (ubuntu-16.04.4-desktop-amd64)**

* ```Anaconda```
	```sh
	$ wget https://repo.continuum.io/archive/Anaconda2-5.1.0-Linux-x86_64.sh
	$ bash Anaconda2-5.1.0-Linux-x86_64.sh
	$ conda install biopython
	```

* ```Java```
	```sh
	$ sudo apt-get update
	$ sudo apt-get install default-jre
	```

**Windows**

* ```Anaconda```
	
	* Download from https://www.anaconda.com/download/
	* Install biopython
	```sh
	$ conda install biopython
	```

* ```Cygwin```

	* Add "C:\cygwin64\bin\" to environment variables-->system variables-->PATH
	* Add "export PATH=/cygdrive/c/Users/user/Anaconda2:$PATH" to .bashrc

* ```Java```

	* Download from https://java.com/en/download/

Execution
==========
Use the following command to see all available options in ProteoStorm
```sh
python -u ./src/ProteoStorm.py --help
```

**Input:**
```
--Database # Directory of protein database files in fasta format.
--Spectra # Directory of spectral datasets in MGF format, peak-picked and converted using [MSConvert](http://proteowizard.sourceforge.net/tools.shtml)
--RemoveSpectra # File containing spectra filename and scan numbers to remove. See ./example/HS_matched_spectra.txt.
--output # Main directory for ProteoStorm output
--SpectralDataset # Name of metaproteomics dataset. Will be used as subdirectory name.
--PepfilterEXE # Path to core module 2 (peptide filtering) executable. Use pre-compiled binaries (CoreModule2_PeptideFiltering_Windows-x86_64.exe, CoreModule2_PeptideFiltering_Linux-x86_64.exe).
--PrecursorMassTolerance # Precursor mass tolerance in ppm (e.g., 10 for HCD).
--FragmentMassTolerance # Fragment mass tolerance in Daltons (e.g., 0.015 for HCD).
--InstrumentID # Identifier of the instrument to generate MS/MS spectra. Used to determine the scoring model in core module 3 (p-value computation). 0: Low-res LCQ/LTQ, 1: Orbitrap/FTICR, 3: Q-Exactive (Default).
--FragmentMethodID # Fragmentation method identifier. Used to determine the scoring model in core module 3 (p-value computation). 1: CID, 3: HCD (Default).
--pval_computation_jar # Path to p-value computation jar.
--CygwinPATH # Required for Windows users (e.g., C:/cygwin64/bin/run.exe)
--GeneraRestrictionApproach # 0: Do not use genera-restriction approach. 1: Use genera-restriction approach. Place UniProt fasta files in /[Database_dir]/UniProt_fasta, and RefSeq fasta files in /[Database_dir]/RefSeq_fasta.
--RefSeqCatalog # Path to RefSeq Catalog.
```

**Output:**
Files are located in [output]/[SpectralDataset]/S2_OutputFiles/
```
- ProteoStorm_output.txt: Peptide-spectrum matches (PSMs) identified after stage two of ProteoStorm, with p-values computed using the MS-GF+ generating function.
- Pooled_0.01_pepFDR.tsv: Peptides passing 1% peptide-level FDR (pooled)
- Pooled_0.01_pepFDR_nonpooled_0.01_psmFDR.tsv: PSMs passing both a 1% PSM-level FDR (per MS/MS experiment) and a 1% peptide-level FDR (pooled)
```

Demo
==========
**1.** Download and extract [demo files](https://drive.google.com/file/d/13k0VANfTPdeLEQ2beZ6Uu5DNaGyJEwiS/view?usp=sharing) into the ProteoStorm main directory.<br />
**2.** Run either demo command 1 or 2 (genera-restriction approach) as provided below. For Windows users, refer to NOTE below.<br />
**3.** The expected output files for the test runs are located at /example/ProteoStorm_Out/prerun_demo and /example/ProteoStorm_Out_GeneraRestrictionApproach/prerun_demo.

**NOTE for Windows users:** <br />

* Append Cygwin path parameter to command (e.g., --CygwinPATH "C:/cygwin64/bin/run.exe")<br />
* Specify --PepfilterEXE as ./src/CoreModule2_PeptideFiltering_Windows-x86_64.exe


**command 1**
```
python -u ./src/ProteoStorm.py --Database ./example/fasta --PartitionMassWindow 15 --Spectra ./example/mgf --SpectralDataset "demo" --RemoveSpectra ./example/HS_matched_spectra.txt --PepfilterEXE ./src/CoreModule2_PeptideFiltering_Linux-x86_64.exe --S1SharedPeakCount 7 --S2SharedPeakCount 6 --output ./example/ProteoStorm_Out --PrecursorMassTolerance 10 --FragmentMassTolerance 0.015 --InstrumentID 3 --FragmentMethodID 3 --pval_computation_jar ./src/MSGFPlus_pvalue.jar --aminoacid_freq ./src/364106_IL_transformed.fasta
```

**command 2 (genera-restriction approach)**
```
python -u ./src/ProteoStorm.py --Database ./example/fasta_genera_restriction_approach --PartitionMassWindow 15 --Spectra ./example/mgf --SpectralDataset "demo" --RemoveSpectra ./example/HS_matched_spectra.txt --PepfilterEXE ./src/CoreModule2_PeptideFiltering_Linux-x86_64.exe --S1SharedPeakCount 7 --S2SharedPeakCount 6 --output ./example/ProteoStorm_Out_GeneraRestrictionApproach --PrecursorMassTolerance 10 --FragmentMassTolerance 0.015 --InstrumentID 3 --FragmentMethodID 3 --GeneraRestrictionApproach 1 --RefSeqCatalog ./example/fasta_genera_restriction_approach/RefSeq-release85_SUB.catalog --pval_computation_jar ./src/MSGFPlus_pvalue.jar --aminoacid_freq ./src/364106_IL_transformed.fasta
```
