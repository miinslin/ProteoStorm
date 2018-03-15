# -*- coding: utf-8 -*-

# Input should only have one PSM per scan.
def Compute_FDR(psmlist, user_FDR,PSMscore_col,FDR_calc_score_direction, prot_col,decoy_prefix):
    FDR_cutoffscore = 'na'
    target_peptides = []
    decoy_PSMs = []
    target_scores = []
    decoy_scores = []
    
    #Store target and decoy scores
    for PSM in psmlist:
        proteins = PSM[prot_col]
        PSMscore = float(PSM[PSMscore_col])
        if decoy_prefix not in proteins:
            target_scores.append(PSMscore)
            target_peptides.append((PSMscore,PSM))
        else:
            decoy_PSMs.append((PSMscore,PSM))
            decoy_scores.append(PSMscore)
    # Dictionary with scores (either target or decoy) as keys, and values as 't', 'd', or 'dt'
    scores = {}
    # scores should be sorted in order from worst to best
    if FDR_calc_score_direction ==1:
        decoy_scores = sorted(decoy_scores)
    elif FDR_calc_score_direction ==-1:
        decoy_scores = sorted(decoy_scores, reverse=True)
    # Start calculating #decoy/#target from first decoy score
    firstdecoyscore = decoy_scores[0]
    
    #=======pre-calculate the number of PSMs passing score x======#
    totalnumdecoy = len(decoy_scores)
    # number of decoy PSMs passing score x    
    numdecoy = {}
    # number of decoy PSMs with the same score
    numsamescoredecoy = {}
    for index, decoy_s in enumerate(decoy_scores):
        if decoy_s in numdecoy:
            numsamescoredecoy[decoy_s]+=1
        elif decoy_s not in numdecoy:
            numdecoy[decoy_s] = totalnumdecoy-index
            scores[decoy_s] = 'd'
            numsamescoredecoy[decoy_s]=1
    # scores should be sorted in order from worst to best
    if FDR_calc_score_direction ==1:
        target_scores = sorted(target_scores)
    elif FDR_calc_score_direction ==-1:
        target_scores = sorted(target_scores, reverse=True)
    firsttargetscore = target_scores[0]
    totalnumtarget = len(target_scores)
    # number of targets PSMs passing score x  
    numtarget = {}
    # number of target PSMs with the same score
    numsamescoretarget = {}
    for index, target_s in enumerate(target_scores):
        if target_s in numtarget:
            numsamescoretarget[target_s]+=1
        elif target_s not in numtarget:
            numtarget[target_s] = totalnumtarget-index
            numsamescoretarget[target_s]=1
            if target_s in scores:
                scores[target_s]='dt'
            elif target_s not in scores:
                scores[target_s]='t'
    scores = scores.items()
    if FDR_calc_score_direction ==1:
        scores = sorted(scores, key=lambda x: x[0])
    elif FDR_calc_score_direction==-1:
        scores = sorted(scores, key=lambda x: x[0], reverse=True)
    #keep score for previous 't' and 'd'. Necessary when current score does not include both a target and decoy
    prevd = firstdecoyscore*1.0
    prevt = firsttargetscore*1.0 
    
    skipfile=0
    begin=-1
    for s in scores:
        if begin ==-1:
            if s[0]==firstdecoyscore:
                begin=1
        if begin==1:
            if s[1]=='t':
                prevt = s[0]
                pass_target_scores = numtarget[s[0]]
                pass_decoy_scores = numdecoy[prevd]-numsamescoredecoy[prevd]
            elif s[1]=='d':
                prevd = s[0]
                pass_decoy_scores = numdecoy[s[0]]
                pass_target_scores = numtarget[prevt]-numsamescoretarget[prevt]
            elif s[1] =='dt':
                prevt = s[0]
                prevd = s[0]
                pass_decoy_scores = numdecoy[s[0]]
                pass_target_scores = numtarget[s[0]]
            
            # if number of passing target scores =0, already reached last target score
            if pass_target_scores==0 and pass_decoy_scores>0:
                skipfile =1
                break   

            currentrate = (1.0*pass_decoy_scores)/pass_target_scores
            
            if currentrate<=user_FDR:
                FDR_cutoffscore = s[0] 
                break

    if skipfile ==1 or FDR_cutoffscore=='na':
        return 'cannot reach specified FDR'
#        raise ValueError('Could not reach specified FDR...')
    
    if FDR_calc_score_direction ==1:
        decoy_PSMs = [z[1] for z in decoy_PSMs if z[0]>=FDR_cutoffscore]
        target_peptides = [z[1] for z in target_peptides if z[0]>=FDR_cutoffscore]
    elif FDR_calc_score_direction==-1:
        decoy_PSMs = [z[1] for z in decoy_PSMs if z[0]<=FDR_cutoffscore]
        target_peptides = [z[1] for z in target_peptides if z[0]<=FDR_cutoffscore]
    
    return (target_peptides, decoy_PSMs, \
    repr(FDR_cutoffscore), \
    str(len(decoy_PSMs))+'|'+str(len(target_peptides)),\
    repr((1.0*len(decoy_PSMs))/len(target_peptides)))