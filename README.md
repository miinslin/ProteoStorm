# ProteoStorm v06222018 #

Software Requirements
---------------

[MSConvert](http://proteowizard.sourceforge.net/tools.shtml) required for peak-picking and converting RAW files to MGF format.

### Linux ###
1. ```Anaconda 2.7``` ***or*** ```Python2.7``` with numpy <br />
```sh
$ wget https://repo.continuum.io/archive/Anaconda2-5.1.0-Linux-x86_64.sh
$ bash Anaconda2-5.1.0-Linux-x86_64.sh
```

2. ```Java (1.8 or above)```<br />
```sh
$ sudo apt-get update
$ sudo apt-get install default-jre
```

### Windows ### 
1. ```Anaconda 2.7```<br />
```sh
download from https://www.anaconda.com/download/?lang=en-us
```

2. ```Java (1.8 or above)```<br />
```sh
download from https://java.com/en/download/
```

3. ```Cygwin```<br />
```sh
download from https://cygwin.com/install.html
Add "C:\cygwin64\bin\" to environment variables-->system variables-->PATH
Add "export PATH=/cygdrive/c/Users/user/Anaconda2:$PATH" to .bashrc
```

Usage
---------------
```sh
python -u ./src/ProteoStorm.py
	-D DatabaseDirectory (Directory containing database files in fasta format)
	-S SpectralDirectory (Directory containing spectral datasets in MGF format, peak-picked and converted from RAW using MSConvert)
	-O OutputDirectory
	[-MSMS SubdirectoryName] (Name of metaproteomics dataset, Default: date_time)
	[-ms1t PrecursorMassTolerance] (e.g., 10, 20, 50, Default: 10)
	[-ms2t FragmentMassTolerance] (e.g., 0.015, 0.6 Default: 0.015)
	[-inst MS2DetectorID] (0: Low-res LCQ/LTQ, 1: Orbitrap/FTICR, 2: TOF, 3: Q-Exactive(Default))
	[-m FragmentMethodID] (1: CID, 3: HCD)
	[-S1spc S1SharedPeaksCount] (Default: 7)
	[-S2spc S2SharedPeaksCount] (Default: 7)
	[-genera 0/1] (0: Create refined protein DB using peptide-level FDR (Default), 1: genera-restriction approach)
	
Output: ProteoStorm_output.txt (Peptide-spectrum matches (PSMs) with p-values computed using the MS-GF+ generating function.)
```
	
Demo
---------------
**1.** Download and extract [demo files](https://drive.google.com/open?id=1LGNSdR4n0mNlztmIryi9SeLoy9420ugy) into ./ProteoStorm/example<br />
**2.** Download [Mass distribution files](https://drive.google.com/open?id=1FfQFeheN2BUSZ9BaA-t0PyspKpYIJrwn) into ./ProteoStorm/src/DBmassDistributions<br />
**3.** Run either Command 1 or 2 (genera-restriction approach) as provided below.<br />
**4.** The expected output files for the test runs are located at ./ProteoStorm/example/ProteoStorm_Out/prerun_demo and ./ProteoStorm/example/ProteoStorm_Out_GeneraRestrictionApproach/prerun_demo.

### Linux ### 

Command 1
```
python -u ./src/ProteoStorm.py --Database ./example/fasta --Spectra ./example/mgf --RemoveSpectra ./example/HS_matched_spectra.txt --SpectralDataset "demo" --output ./example/ProteoStorm_Out --PrecursorMassTolerance 10 --FragmentMassTolerance 0.015 --InstrumentID 3 --FragmentMethodID 3
```
Command 2 (genera-restriction approach)
```
python -u ./src/ProteoStorm.py --Database ./example/fasta_genera_restriction_approach --Spectra ./example/mgf --RemoveSpectra ./example/HS_matched_spectra.txt --SpectralDataset "demo" --output ./example/ProteoStorm_Out_GeneraRestrictionApproach --PrecursorMassTolerance 10 --FragmentMassTolerance 0.015 --InstrumentID 3 --FragmentMethodID 3 --GeneraRestrictionApproach 1 --refDBfdr 0.01
```

### Windows ### 
***Replace "C:/cygwin64/bin/run.exe" with corresponding path in system.***

Command 1
```
python -u ./src/ProteoStorm.py --Database ./example/fasta --Spectra ./example/mgf --RemoveSpectra ./example/HS_matched_spectra.txt --SpectralDataset "demo" --output ./example/ProteoStorm_Out --PrecursorMassTolerance 10 --FragmentMassTolerance 0.015 --InstrumentID 3 --FragmentMethodID 3 --CygwinPATH "C:/cygwin64/bin/run.exe"
```
Command 2 (genera-restriction approach)
```
python -u ./src/ProteoStorm.py --Database ./example/fasta_genera_restriction_approach --Spectra ./example/mgf --RemoveSpectra ./example/HS_matched_spectra.txt --SpectralDataset "demo" --output ./example/ProteoStorm_Out_GeneraRestrictionApproach --PrecursorMassTolerance 10 --FragmentMassTolerance 0.015 --InstrumentID 3 --FragmentMethodID 3 --GeneraRestrictionApproach 1 --refDBfdr 0.01 --CygwinPATH "C:/cygwin64/bin/run.exe"
```

Alternative configurations for ProteoStorm
---------------
[Alternative configurations for ProteoStorm PDF download](https://drive.google.com/file/d/17y7JxIy7duakiXHIFQ9P1IlMX1gMtM-i/view?usp=sharing)