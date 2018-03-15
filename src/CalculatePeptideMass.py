import re 

#msgfplus/src/main/java/edu/ucsd/msjava/msutil/Composition.java
element_mass = {'C':12.0,
                'H':1.007825035,
                'N':14.003074,
                'O':15.99491463,
                'S':31.9720707,
                'P':30.973762,
                'Br':78.9183361,
                'Cl':34.96885272,
                'Fe':55.9349393,
                'Se':79.9165196
                }

def FormulaToMass(mass_composition):
    mass_composition = [x for x in re.split(r'(-?\d+)', mass_composition) if x]
    mass_composition = [(mass_composition[i],mass_composition[i+1]) for i in range(0,len(mass_composition),2)]
    return round(sum([element_mass[element[0]]*int(element[1]) for element in mass_composition]),10)

#http://www.unimod.org/modifications_list.php
mods_chemical_formula = {'Carbamidomethyl':'C2H3N1O1',
             'Oxidation':'O1',
             'Deamidated':'H-1N-1O1',
             'Carbamyl':'H1C1N1O1',
             'Glu->pyro-Glu':'H-2O-1',
             'Gln->pyro-Glu':'H-3N-1',
             'Acetyl':'C2H2O1',
             'Methyl':'H2C1',
             'Phospho':'H1O3P1',
             'Formyl':'C1O1'
             }

#mods_monoisotopic_masses = {}
#for x in mods_chemical_formula:
#    mods_monoisotopic_masses[x] = FormulaToMass(mods_chemical_formula[x])

mods_monoisotopic_masses = {'Acetyl': 42.0105647,
 'Carbamidomethyl': 57.021463735,
 'Carbamyl': 43.005813665,
 'Deamidated': 0.984015595,
 'Formyl': 27.99491463,
 'Gln->pyro-Glu': -17.026549105,
 'Glu->pyro-Glu': -18.0105647,
 'Methyl': 14.01565007,
 'Oxidation': 15.99491463,
 'Phospho': 79.966330925}

aa_chemical_formula = {'A': 'C3H5N1O1',
'C': 'C3H5N1O1S1', 
'E': 'C5H7N1O3', 
'D': 'C4H5N1O3', 
'G': 'C2H3N1O1', 
'F': 'C9H9N1O1', 
'I': 'C6H11N1O1', 
'H': 'C6H7N3O1', 
'K': 'C6H12N2O1', 
'M': 'C5H9N1O1S1', 
'L': 'C6H11N1O1', 
'N': 'C4H6N2O2', 
'Q': 'C5H8N2O2', 
'P': 'C5H7N1O1', 
'S': 'C3H5N1O2', 
'R': 'C6H12N4O1', 
'T': 'C4H7N1O2', 
'W': 'C11H10N2O1', 
'V': 'C5H9N1O1', 
'Y': 'C9H9N1O2'}

#aa_monoisotopic_masses = {}
#for x in aa_chemical_formula:
#    aa_monoisotopic_masses[x] = FormulaToMass(aa_chemical_formula[x])

aa_monoisotopic_masses = {'A': 71.037113805,
 'C': 103.009184505,
 'D': 115.026943065,
 'E': 129.042593135,
 'F': 147.068413945,
 'G': 57.021463735,
 'H': 137.058911875,
 'I': 113.084064015,
 'K': 128.09496305,
 'L': 113.084064015,
 'M': 131.040484645,
 'N': 114.04292747,
 'P': 97.052763875,
 'Q': 128.05857754,
 'R': 156.10111105,
 'S': 87.032028435,
 'T': 101.047678505,
 'V': 99.068413945,
 'W': 186.07931298,
 'Y': 163.063328575}

def calcmass_cmm(pepstring, TMT_mod = 0):
    if TMT_mod ==0:
        pepmass = element_mass['H']+element_mass['O']+element_mass['H']
        for aa in pepstring:
            # add mass of amino acid
            pepmass += aa_monoisotopic_masses[aa]
            # default have fixed C Carbamidomethyl modification
            # C has side chain CH2-SH
            # iodoacetamide reacts with reduced cysteine, prevents the formation of disulfide bonds
            if aa == 'C':
                pepmass += mods_monoisotopic_masses['Carbamidomethyl']                
        return pepmass
    
    # TMT tags composed of amine-reactive NHS-ester group, a spacer arm and a mass reporter. 
    # Lysine side-chain NH2 or n-terminus NH2
    # if C at n-terminus is carbamidomethylated, TMT can still react with NH2
    elif TMT_mod ==1:
        TMTmodssum = element_mass['H']+element_mass['O']+element_mass['H']
        for i, aa in enumerate(pepstring):
            # add mass of amino acid
            TMTmodssum += aa_monoisotopic_masses[aa]
            # if n-terminus of peptide, add 229.162932
            if i==0:
                TMTmodssum += 229.162932
                # if n-terminus aa is also a K, add 229.162932 again (lysine side-chain)
                if aa == 'K':
                    TMTmodssum += 229.162932
                # if n-terminus aa is a C, side-chain can still be modified
                if aa == 'C':
                    TMTmodssum += mods_monoisotopic_masses['Carbamidomethyl']
            # if not aa at n-terminus of peptide but a K, add 229.162932
            if i!=0:
                if aa == 'K':
                    TMTmodssum += 229.162932
                # fixed C Carbamidomethyl modification
                if aa == 'C':
                    TMTmodssum += mods_monoisotopic_masses['Carbamidomethyl']
        return TMTmodssum


#def parse_modfile(modfilename):
#    # parse modificatoins file
#    modslist = {}
#    with open(modfilename, 'r') as infile:
#        modid = -1
#        for line in infile:
#            if line[0]!='#' and line.strip()!='':
#                if 'NumMods' in line.strip():
#                    NumMods = line.strip().split('NumMods=')[1]
#                #Mass or CompositionStr, Residues, ModType, Position, Name (all the five fields are required)
#                else:
#                    item = line.strip().split()[0].split(',')
#                    # CompositionStr (C[Num]H[Num]N[Num]O[Num]S[Num]P[Num]Br[Num]Cl[Num]Fe[Num])
#                    CompositionStr = item[0]
#                    if any(aa in CompositionStr for aa in ['C','H','N','O','S','P','Br','Cl','Fe']):
#                        CompositionMass = FormulaToMass(CompositionStr)
#                    else:
#                        CompositionMass = round(float(CompositionStr),6)
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
#                    modslist[modid] = {'mass':CompositionMass,'residues':residues, 'modtype':modtype, 'position':position}
#    return modslist
#
## function to determine modified peptide sequence
#def findmodifiedpeptides(modslist, peptide, pepmass, preaa, postaa):
#    candidate_modified_peplist_all = {}
#    
#    for mod in modslist:  
#        position = modslist[mod]['position']
#        residues = modslist[mod]['residues']
#        mass = modslist[mod]['mass']
#
#        modified_mass = pepmass+mass
#        
#        prefix = ''
#        if -1.0*mass<0:
#            prefix = '+'
#        
#        candidate_modified_peplist = []
#        # if modification can be on any residue at position...
#        if residues == '*':
#            if position == 'N-term':
#                candidate_modified_peplist = [prefix+str(mass)+peptide]
#            elif position == 'C-term':
#                candidate_modified_peplist = [peptide+prefix+str(mass)]
#            elif position == 'Prot-N-term':
#                if preaa == '-':
#                    candidate_modified_peplist = [prefix+str(mass)+peptide]
#            elif position == 'Prot-C-term':
#                if postaa == '-':
#                    candidate_modified_peplist = [peptide+prefix+str(mass)]
#        # if specific residues to be modified...
#        elif residues != '*':
#            if position == 'any':
#                candidate_modified_peplist = [peptide[:ai+1]+prefix+str(mass)+peptide[ai+1:] for ai, aa in enumerate(peptide) if aa in residues]
#            elif position == 'N-term':
#                if peptide[0] in residues:
#                    candidate_modified_peplist = [prefix+str(mass)+peptide]
#            elif position == 'C-term':
#                if peptide[-1] in residues:
#                    candidate_modified_peplist = [peptide+prefix+str(mass)]
#            elif position == 'Prot-N-term':
#                if preaa == '-' and peptide[0] in residues:
#                    candidate_modified_peplist = [prefix+str(mass)+peptide]
#            elif position == 'Prot-C-term':
#                if postaa == '-' and peptide[-1] in residues:
#                    candidate_modified_peplist = [peptide+prefix+str(mass)]
#        
#        if candidate_modified_peplist:
#            if modified_mass not in candidate_modified_peplist_all:
#                candidate_modified_peplist_all[modified_mass] = set() 
#            candidate_modified_peplist_all[modified_mass].update(candidate_modified_peplist)
#        
#    return candidate_modified_peplist_all