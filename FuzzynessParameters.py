def FuzzyParameters():
    #We ensure that number of big mistakes allowed BM is the maximum number of mistakes 
    #that can be tolerated in the matching (it could be less in the result).
    #Big mistakes contain non-canonical interaction missing 
    #and labels/non-canonical interactions distant from more than eleminatory_isostericity_threshold
    number_BM_allowed = 0#5
    
    #Fixed the limit of isostericity that eliminate the matching if overcomes.
    elim_iso_threshold = 10 
    
    #Fixed the limit of isostericity in which it is just a big mistake.
    #Big mistakes in isostericity are the interval: 
    #[BM_iso_threshold, elim_iso_threshold[
    BM_iso_threshold = 4
    
    #If this boolean is 1 it is eliminatory, otherwise it is just a big mistake.
    B53_missing_elim = 1
    
    #If small mistakes (nonBM) are allowed, isostericity below big_mistake_isostericity_threshold
    #will also be softly taken into account into the global error count.
    allow_iso_nonBM = 1

    #We define if edges missing are consider as a big mistakes or are eliminatory -1 means eliminatory
    #Others gives the number of bg mistakes that are counted for such a mistake.
    BM_by_missing_edge = -1

    #Gaps are counted as multiple big mistakes whose number can be specified here.
    BM_by_gap = 1

    #If it is 1, number_BM_allowed is for a pattern of around 20 edges.
    # #In this case, scaling is done according to the size of each pattern.
    scale_BM_with_P_size=0
    return number_BM_allowed, elim_iso_threshold, BM_iso_threshold, B53_missing_elim, allow_iso_nonBM, BM_by_missing_edge, BM_by_gap, scale_BM_with_P_size