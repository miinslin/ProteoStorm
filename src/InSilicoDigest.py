# -*- coding: utf-8 -*-

def insilicodigest(sequence, regex_enz, min_pep_len, max_pep_len, miscleavages, ntt, td, REMOVE_KRP):
    #determine if KRP semipeptide
    if ntt=='semi':
        def nonKRPpep(pepstring, preaastr, postaastr, anchortype):
            if anchortype == 'right':
                # -.P is fully tryptic
                if preaastr in ['K','R'] and pepstring[0] == 'P':
                    return False
                else:
                    return True
                #
            elif anchortype == 'left':
                if pepstring[-1] in ['K','R'] and postaastr == 'P':
                    return False
                else:
                    return True
    
    # rules for length and non-standard amino acids
    def passrules(pepseq):
        if min_pep_len<=len(pepseq)<=max_pep_len:
            if any(aa in pepseq for aa in 'BJOUXZ?-'):
                return False
            else:
                return True
        else:
            return False
    
    # given concatenated protein sequence $XXXXXXXXX*$XXXXXXXXXXX*$XXXXXXXXXXXXXXXXXXXX*,
    # where $ indicates start of protein, and * indicates end of protein, find cut sites
    # given trypsin rules
    cutsites = regex_enz.finditer(sequence)
    trypsites = [0]
    trypsites.extend([pep.end() for pep in cutsites])
    trypsites.append(len(sequence))

    # trypsites[i] is where peptide sequence starts in allprot. Required to determine which protein in fasta chunk later.
    # (pepsequence, trypsite, pre+post)
    pep_list = {}

    # given start sites, regenerate peptide, and determine pre and post amino acids.
    # If semi-tryptic, also generate semi-tryptic peptides from tryptic peptide. 
    for i in range(len(trypsites)-1):
        mcleaved_peptide = 'notmiscleavage'
        # set of peptides         
        add_p = set()
        # peptide sequence
        pepseq = sequence[trypsites[i]:trypsites[i+1]].replace('*','').strip()

        if not pepseq:
            continue
        
        # ==== determine pre post amino acids ===== #
        # if beginning of protein, pre aa is -
        if sequence[trypsites[i]]=='$':
            if sequence[trypsites[i+1]-1]=='*' or sequence[trypsites[i+1]]=='*':
                prepost = '--'
            else:
                prepost = '-'+sequence[trypsites[i+1]]
        # if end of protein, post aa is -
        elif sequence[trypsites[i+1]-1]=='*' or sequence[trypsites[i+1]]=='*':
            prepost = sequence[trypsites[i]-1]+'-'
        # else grab pre post aa
        else:
            prepost = sequence[trypsites[i]-1]+sequence[trypsites[i+1]]
        
        # if start of protein, and if first amino acid is Met
        if pepseq[0]=='$':
            #remove $ sign
            pepseq = pepseq[1:]
            if pepseq[0]=='M':
                #add m-cleaved version of peptide
                add_p.add((pepseq[1:], prepost))
                mcleaved_peptide = pepseq
        
        # add peptide to set
        add_p.add((pepseq,prepost))

        # if miscleavage >0
        if miscleavages!=0:            
            for mi in range(1,miscleavages+1):           
                # if miscleavage peptide breaches total concatenated protein string
                try:
                    trypsites[i+1+mi]
                except IndexError:
                    break
                # miscleavage peptide sequence
                mc_pep = sequence[trypsites[i]:trypsites[i+1+mi]].strip()   
                
                # if miscleavage peptide includes sequences from >1 protein ,move on
                if '*$' in mc_pep:
                    break                

                # if miscleavage peptide only from single protein                
                elif '*$' not in mc_pep:
                    # == determine pre post aa
                    # if peptide at the end of concatenated protein string, post aa is -
                    if trypsites[i+1+mi] == len(sequence):
                        mc_pep_prepost = prepost[0] + '-'
                    # if peptide not at end of concatenated protein string....
                    if trypsites[i+1+mi] != len(sequence):
                        # if last peptide of protein, post aa is -
                        if sequence[trypsites[i+1+mi]-1]=='*' or sequence[trypsites[i+1+mi]]=='*':
                            mc_pep_prepost = prepost[0] +'-'
                        # else take pre aa from miscleavage == 0 peptide, and find post aa
                        else:
                            mc_pep_prepost = prepost[0] + sequence[trypsites[i+1+mi]]                    
                    
                    # remove '*' 
                    mc_pep = mc_pep.replace('*','')
                    # if start of protein, and if first amino acid is Met
                    if mc_pep[0]=='$':
                        #remove $ sign
                        mc_pep = mc_pep[1:]
                        if mc_pep[0] == 'M':
                            #add m-cleaved version of peptide
                            add_p.add((mc_pep[1:], mc_pep_prepost))

                    # add peptide
                    add_p.add((mc_pep, mc_pep_prepost))


        #==== add peptides, cutsite information, and pre post aa to pep_list ===#
        # if fully-tryptic digest
        if ntt =='full':
            # check if peptide passes rule              
            for n in add_p:
                if passrules(n[0]):
                    if n[0] in pep_list:
                        # peptide key: (start index, pre+post)
                        pep_list[n[0]].add((td+str(trypsites[i]), n[1]))
                    else:
                        pep_list[n[0]] = set([(td+str(trypsites[i]), n[1])]) 
            del add_p
            continue

        # if semi-tryptic digest
        if ntt=='semi':
            # create semi-tryptic versions of peptide, and check if peptide passes rule               
            for n in add_p:
                for z in range(len(n[0])):     
                    # RIGHT ANCHOR
                    semipepseq = n[0][z:]
                    postaa = n[1][1]
                    preaa = n[0][z-1]
                    # first is always fully-tryptic peptide
                    z_indx = [0]
                    # if peptide starts with $M, should consider M cleaved as fully-tryptic
                    if n[0][:len(mcleaved_peptide)] == mcleaved_peptide:
                        z_indx = [0,1]
                    # if fully-tryptic on left side, preaa should =='-'
                    if z in z_indx:
                        preaa = n[1][0]
                    # check if peptide is a K/R.P peptide 
                    # in which case we will ignore because p value/specEvalue calculation 
                    # cannot differentiate K/R.P peptides as semi-tryptic  
                    keepKRP = True
                    if REMOVE_KRP ==1:
                        keepKRP = nonKRPpep(semipepseq, preaa, postaa, 'right')
                    if keepKRP:
                        # only include standard amino acid peptides wtihin length requirements
                        if passrules(semipepseq):
                            if semipepseq in pep_list:
                                pep_list[semipepseq].add((td, preaa+postaa))
                            else:
                                pep_list[semipepseq] = set([(td, preaa+postaa)])

                    # LEFT ANCHOR
                    if z in range(len(n[0])-1):
                        semipepseq = n[0][:z+1]
                        preaa = n[1][0]
                        postaa = n[0][z+1]
                        keepKRP = True
                        if REMOVE_KRP ==1:
                            keepKRP = nonKRPpep(semipepseq, preaa, postaa, 'left')
                        if keepKRP:
                            #check if peptides pass rules
                            if passrules(semipepseq):
                                # td is 'd_' or ''
                                if semipepseq in pep_list:
                                    pep_list[semipepseq].add((td, preaa+postaa))
                                else:
                                    pep_list[semipepseq] = set([(td, preaa+postaa)])
            del add_p
            
    return pep_list