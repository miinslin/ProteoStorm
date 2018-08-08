# -*- coding: utf-8 -*-
import os, re, time
import numpy as np
#Input mgf files have to be converted from RAW using MSConvert

def CreateSpectralParts(spectra_dir, spectralparts_dir, \
spectra_remove, dbp_ranges, num_Spectra, PSlogfile, precursormasstol,\
massbin_n):
    
    idxcount = [0 for x in xrange(massbin_n)]
    
    Buffsize = 209715200
    
    # CRC_Handbook_of_Chemistry_and_Physics
    # mass of proton 1.00727646688
    # msgfplus/src/main/java/edu/ucsd/msjava/msutil/Composition.java
    proton_mass = 1.00727649
    # mass difference between C12 and C13 is 1.003355
    isotope = 1.003355
#    multiplier = 0.9995
    
    print 'Spectral partitioning...'
    
    start_time = time.time()
    
    HS_spectra = {}
    if not os.path.exists(spectra_remove) and spectra_remove!='na':
        PSlogfile.write('Error: Spectra to remove should either be a file or "na".'+'\n')
        raise ValueError('spectra to remove parameter should either be a file or na')
    if os.path.exists(spectra_remove):
        with open(spectra_remove,'r') as infile:
            for line in infile:
                HS_list = line.strip().split()
                HS_spectra[HS_list[0]] = set([int(x) for x in HS_list[1:]])

    # uses 10 ppm mass tolerance
    def mmrange_given_tmrange(tmrange):
        s = tmrange[0]*(1-(precursormasstol/(10**6)))
        e = (tmrange[1]*(1+(precursormasstol/(10**6))))+1*(isotope)
        return (s,e)

    dbp_ranges = [mmrange_given_tmrange((x[0],x[1])) for x in dbp_ranges]
    
    def numpeaks(mzlist):
        return sum([1 for x in mzlist[1:-1] if '=' not in x])
    
    def writespectratofile(sorted_pepmass, spectra, idxcount):
        # list of (spectralindex,precursormass)
        sorted_pepmass.sort(key=lambda x: x[1])
        pmass_only = np.array([pmass[1] for pmass in sorted_pepmass])

        for fastapart, fp_mmrange in enumerate(dbp_ranges):
            spectraoutfile = os.path.join(spectralparts_dir, 'SpecPart_'+str(fastapart)+'_'+str(idxcount[fastapart])+'.mgf')
            leftindex = np.searchsorted(pmass_only,fp_mmrange[0], side='left')
            rightindex = np.searchsorted(pmass_only,fp_mmrange[1], side='right')
            
            if not sorted_pepmass[leftindex:rightindex]:
                continue
            
            if os.path.exists(spectraoutfile):
                #check file size
                current_filesize = os.stat(spectraoutfile).st_size
                if current_filesize>Buffsize:
                    idxcount[fastapart]+=1
                    spectraoutfile = os.path.join(spectralparts_dir, 'SpecPart_'+str(fastapart)+'_'+str(idxcount[fastapart])+'.mgf')
                else:
                    outfile = open(spectraoutfile,'a')
            if not os.path.exists(spectraoutfile):
                outfile = open(spectraoutfile,'w')
            
            for sp in sorted_pepmass[leftindex:rightindex]:
                outfile.write('\n'.join(spectra[sp[0]])+'\n')
            
            outfile.close()
        return ''
    
    sorted_pepmass = []
    spectra = []
    spectralindex = -1
    totalnumspecparsed = 0

    for SpecFileName in os.listdir(spectra_dir):
        if SpecFileName not in HS_spectra:
            matched_spec = []
        elif SpecFileName in HS_spectra:
            matched_spec = HS_spectra[SpecFileName]           
        with open(os.path.join(spectra_dir,SpecFileName),'r') as specfile:
            begin = 1
            for line in specfile:
                if line.strip() =='BEGIN IONS':
                    begin = 0
                    charge = 0
                    spectrum = []
                    spectrum.append(line.strip())
                    continue
                if begin==0:
                    if line.strip()[:6] == 'TITLE=':
                        scanNumber = int(line.split('scan=')[1].split('"')[0])
                        if scanNumber in matched_spec:
                            begin = 1
                            continue
                        spectrum.append(line.strip())
                    elif line.strip()[:8] == 'PEPMASS=':
                        precursorMassCharge = float(line.strip().split()[0].split('PEPMASS=')[1])
                        spectrum.append(line.strip())
                    elif line.strip()[:7] == 'CHARGE=':
                        charge = float(re.findall('\d+', line.strip())[0])
                        precursormass = round((precursorMassCharge*charge)-(proton_mass*charge), 4)
                        if precursormass<dbp_ranges[0][0] or precursormass>dbp_ranges[-1][1]:
                            begin =1
                            continue
                        spectrum.append(line.strip())
                    elif line.strip() == 'END IONS':
                        spectrum.append(line.strip())
                        if charge == 0:
                            begin =1
                            spectrum = []
                            continue
                        # if less than 10 peaks, ignore.
                        if numpeaks(spectrum)<10:
                            begin = 1
                            spectrum = []
                            continue
                        spectra.append(spectrum)
                        spectralindex +=1
                        totalnumspecparsed+=1
                        sorted_pepmass.append((spectralindex,precursormass))
                        precursormass = 0
                        charge = 0
                        if spectralindex >= num_Spectra:
                            if spectra:
                                writespectratofile(sorted_pepmass, spectra, idxcount)
                            sorted_pepmass = []
                            spectra = []
                            spectralindex = -1
                    else:
                        spectrum.append(line.strip())
    if spectra:
        writespectratofile(sorted_pepmass, spectra, idxcount)
        sorted_pepmass = []
        spectra = []        

    print '\t',str(totalnumspecparsed), ' spectra partitioned...'
    PSlogfile.write(str(totalnumspecparsed)+ ' spectra partitioned.'+'\n')
    
    if totalnumspecparsed==0:
        PSlogfile.write('Exiting. No spectra partitioned. Check that mgf contains charge.')
        raise ValueError('No spectra partitioned. Check that mgf contains charge.')
    
    return time.time()-start_time
