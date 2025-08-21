"""
gen2_ngrams.py

This is a Python script for generation of normal (i.e., continuous) n-grams and skippy (i.e., discontinuous) n-grams, regular or extended.
This is a re-implementation of gen_ngrams.py and dispenses with gen_extended_skippy_ngrams(..). Now, gen_skippy_ngrams(..) generates extended skippy n-grams by setting extended = True.

Created on 2025/08/20 by Kow Kuroda (kow.kuroda@gmail.com)

Modifications
2025/08/21 i) revised the algorithm to avoid unwanted removal of duplicates; ii) simplified the processing;

"""

##
def gen_ngrams (S: list, n: int, sep: str = " ", as_list: bool = False, check: bool = False):
    """
    takes a list S of segments and returns a list R of n-grams out of them.
    """
    
    assert n > 0
    ##
    if check:
        print(f"#S: {S}")
    
    ##
    segs = [ seg for seg in S if len(seg) > 0 ]
    if len(segs) <= n:
        if as_list:
            return [ segs ]
        else:
            return [ sep.join(segs) ]
    
    ## main
    R = [ ]
    for i in range(len(segs)):
        try:
            gram = segs[i : i + n] # get an n-gram
            if len(gram) == n:
                R.append(gram)
        except IndexError:
            pass
    ##
    if as_list:
        return R
    else:
        return [ sep.join(r) for r in R ]

##
def segment(t: str, pattern: str = r"", as_tuple: bool = False):

    import re
    if as_tuple:
        return ( x for x in re.split(pattern, t) if len(x) > 0 )
    else:
        return [ x for x in re.split(pattern, t) if len(x) > 0 ]

##
def gen_source(L: list, gap_mark: str = "…", as_tuple: bool = False):
    
    if as_tuple:
        return ( (x, y) for x, y in zip(L, gap_mark * len(L)) )
    else:
        return [ (x, y) for x, y in zip(L, gap_mark * len(L)) ]

##
def count_elements(L: list, gap_mark: str = "…"):

    """
    counts the number of non-gap elements in a segment sequence and returns it
    """

    #return len(x for x in L if x != gap_mark) # fails due to generator mishandling
    return len([ x for x in L if x != gap_mark ])

##
def simplify_gaps(segs: list, gap_mark: str, verbose: bool = False, check: bool = False):

    """
    simplifies repeated gap_marks in compliance to the definition of extended skippiness
    """

    ##
    if check:
        print(f"input segs: {segs}")
    
    ## main
    R = [ ]
    R_size = len(R)
    for i, seg in enumerate(segs):
        if check and verbose:
            print(f"seg{i}: {seg}")
        if seg == gap_mark:
            try:
                if R[-1] == gap_mark:
                    pass
                else:
                    R.append(seg)
            except IndexError:
                R.append(seg)
        else:
            R.append(seg)
        ##
        if check:
            if len(R) > R_size:
                print(f"added seg{i}: {seg}")
        R_size = len(R) # update
    
    ##
    if check:
        print(f"R: {R}")
    return R

##
def remove_gaps(segs: list, gap_mark: str):

    """simply removes all gaps in a segment list and returns the result"""
    
    return [ seg for seg in segs if seg != gap_mark ]

##
def gen_skippy_ngrams(L: list, n: int, max_gap_size: int = None, extended: bool = True, inclusive: bool = False, sep: str = " ", gap_mark: str = "…", as_list: bool = False, verbose: bool = False, check: bool = True):

    """
    general generator function that can be called
    """
    
    assert n > 0

    ## filter out empty elements
    segs = [ seg for seg in L if len(seg) > 0 ]
    if len(segs) <= n:
        if as_list:
            return [segs]
        else:
            return [sep.join(segs)]

    ## set search_span
    divided = True
    if max_gap_size is None or len(segs) <= max_gap_size:
        search_span = len(segs)
        divided = False
        if check:
            print(f"search is not divided with span: {search_span} = len(segs): {len(segs)}")
    else:
        search_span = max_gap_size
        if check:
            print(f"search is divided with span: {search_span} < len(segs): {len(segs)}")
    
    ## generate a lattice of segs, with or without subdivision
    import itertools
    if not divided:
        #P = itertools.product(gen_source(segs)) # not work!
        P = list(itertools.product(*gen_source(segs))) # list(..) is needed to retain result
    else:
        subsegs_set = gen_ngrams(segs, search_span, sep = sep, as_list = True)
        set_size = len(subsegs_set)
        P = [ ]
        Px = [ ] # checker of iso-forms
        if check:
            print(f"set_size: {set_size}")
        for i, subsegs in enumerate(subsegs_set):
            if check:
                print(f"{i} subsegs0: {subsegs}")
            ## adjust edges
            if i == 0:
                subsegs = subsegs + [gap_mark]
            elif i > 0 and i < (set_size - 1):
                subsegs = [gap_mark] + subsegs + [gap_mark]
            elif i == (set_size - 1):
                subsegs = [gap_mark] + subsegs
            else:
                pass
            ##
            if check:
                print(f"{i} subsegs1: {subsegs}")
            
            ## filter out iso-forms
            for p in list(itertools.product(*gen_source(subsegs))):
                px = simplify_gaps(p, gap_mark)
                if check:
                    print(f"px: {px}")
                if not p in P and not px in Px:
                    P.append(p)
                    Px.append(px) # Crucial to avoid seemingly repeated occurrences
    ##
    if check:
        print(f"P0: {P}")
    
    ## make result fit to unextended version
    if not extended:
        X = [ ]
        for p in P:
            if p[0] != gap_mark or p[-1] != gap_mark:
                X.append(p)
            elif p[0] == gap_mark and p[-1] == gap_mark:
                X.append(p)
            else:
                pass
        P = X
    if check:
        print(f"P: {P}")
    
    ## regulate gaps
    Q = [ ]
    for i, p in enumerate(P):
        if check:
            print(f"p{i}: {p}")
        ##
        segs_size = count_elements(p)
        if check:
            print(f"segs_size: {segs_size}")
        if not inclusive:
            if segs_size == n:
                q = simplify_gaps(p, gap_mark = gap_mark, check = check)
                Q.append(q)
        else:
            if segs_size <= n:
                q = simplify_gaps(p, gap_mark = gap_mark, check = check)
                Q.append(q)
    
    ## remove gap_mark singleton
    Q = [ p for p in Q if count_elements(p) > 0 ]
    
    ## handle 1-grams
    if not extended:
        Q = [ remove_gaps(p, gap_mark) if count_elements(p) == 1 and len(p) > 1 else p for p in Q ]
    
    ## return
    if check:
        print(f"Q: {Q}")
    if as_list:
        return Q
    else:
        return [ sep.join(x) for x in Q ]

### end of file