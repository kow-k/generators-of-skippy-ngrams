#!/usr/bin/env python3
# -*- encoding: utf-8 -*-

"""
gen2_ngrams.py

This is a Python script for generation of normal (i.e., continuous) n-grams and skippy (i.e., discontinuous) n-grams, regular or extended.

This script is a re-implementation of gen_ngrams.py and dispenses with gen_extended_skippy_ngrams(..). There are two main differences from its predecessor. First, gen_skippy_ngrams(..) generates extended skippy n-grams with option extended = True. Second, gen_skippy_ngrams(..) optionally generates inclusive n-grams by "inclusive = True" option, thereby dispensing with incremental generation of 1-gram up to n-grams afterwards.

Created on 2025/08/20 by Kow Kuroda (kow.kuroda@gmail.com)

Modifications
2025/08/21 1) revised the algorithm to avoid unwanted removal of duplicates; 2) simplified the processing;
2025/08/22 made minor changes;
2025/08/26 fixed bugs in filtering of p from P, where non-extended cases are treated as extended ones;
2025/08/28 implemented a better solution to mishandling of n_for_ngram in recursion;

"""

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
        return ( [x, y] for x, y in zip(L, gap_mark * len(L)) )
    else:
        return [ [x, y] for x, y in zip(L, gap_mark * len(L)) ]

##
def count_elements(L: list, gap_mark: str = "…"):

    """
    counts the number of non-gap elements in a segment sequence and returns it
    """

    #return len(x for x in L if x != gap_mark) # fails due to generator mishandling
    return len([ x for x in L if x != gap_mark ])

##
def count_gaps(L: list, gap_mark: str = "…"):

    """
    counts the number of gaps in a segment sequence and returns it
    """

    return len([ x for x in L if x == gap_mark ])

##
def drop_gap_at_end(segs: list, gap_mark: str):
    
    F = []
    for i, seg in enumerate(segs):
        if seg != gap_mark:
            F.append(seg)
        else:
            if i == 0 or i == (len(segs) - 1):
                pass
            else:
                F.append(seg)
    ##
    return F

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
            print(f"{i} seg: {seg}")
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
        if check and verbose:
            if len(R) > R_size:
                print(f"added seg{i}: {seg}")
        R_size = len(R) # update

    ##
    if check:
        print(f"#R: {R}")
    return R

##
def remove_gaps(segs: list, gap_mark: str):

    """simply removes all gaps in a segment list and returns the result"""

    return [ seg for seg in segs if seg != gap_mark ]

##
## Beware to make recursively = True. it procudes extra strings;
def gen_ngrams (S: list, n_for_ngram: int, inclusive: bool = False, recursively: bool = False, sep: str = " ", as_list: bool = False, check: bool = False):
    """
    takes a list S of segments and returns a list R of n-grams out of them.
    """

    assert n_for_ngram > 0
    if check:
        print(f"#S: {S}")

    ##
    segs = [ seg for seg in S if len(seg) > 0 ]
    if len(segs) <= n_for_ngram:
        if recursively:
            G = gen_ngrams(segs, n_for_ngram - 1, inclusive = inclusive, recursively = recursively, sep = sep, as_list = as_list, check = check)
            if as_list:
                if not segs in G:
                    G = G + segs
                return G
            else:
                g = sep.join(segs)
                if not g in G:
                    G.append(segs)
                return [ sep.join(segs) for segs in G ]
        else:
            if as_list:
                return [ segs ]
            else:
                return [ sep.join(segs) ]

    ## main
    R = [ ]
    if inclusive:
        for j in range(1, n_for_ngram + 1):
            for i in range(len(segs)):
                try:
                    gram = segs[i : i + j] # get an n-gram
                    if len(gram) == j:
                        R.append(gram)
                except IndexError:
                    pass
    else:
        for i in range(len(segs)):
            try:
                gram = segs[i : i + n_for_ngram] # get an n-gram
                if len(gram) == n_for_ngram:
                    R.append(gram)
            except IndexError:
                pass
    ##
    if as_list:
        return R
    else:
        return [ sep.join(r) for r in R ]

##
def filter_segs(subsegs_pool: list, n_for_ngram: int, max_gap_size: int, extended: bool = True, inclusive: bool = True, gap_mark: str = "…", verbose: bool = False, check: bool = False):
    
    if check and verbose:
        print(f"#max_gap_size: {max_gap_size}")
    import itertools
    Q = [ ]; xQ = [ ] # checker of iso-forms
    for i, subsegs in enumerate(subsegs_pool):
        if check:
            print(f"#{i} subsegs: {subsegs}")
        
        ## process over the segs of a given subsegs
        #for segs in list(itertools.product(*gen_source(subsegs))):
        ## The following code increases Cython-compatibility
        for i, segs in enumerate([ list(x) for x in itertools.product(*gen_source(subsegs)) ]):
            if check:
                print(f"#{i} segs: {segs}")
            
            ## define xsegs for later reference
            n_segs = len(segs)
            n_elements = count_elements(segs, gap_mark)
            n_gaps = count_gaps(segs, gap_mark)
            xsegs = simplify_gaps(segs, gap_mark = gap_mark, check = check)
            
            ## exclude sequences of gap_markers
            if n_elements == 0:
                if check:
                    print(f"#ignored: {segs} [n_elements: {n_elements} == 0]\n...")
                continue
            
            ## excludes segs longer then max_gap_size
            if n_segs > max_gap_size:
                if check:
                    print(f"#ignored: {segs} [n_segs: {n_segs} > max_gap_size: {max_gap_size}]\n...")
                continue
            
            ## excludes if count_elements(p) > n_for_ngram
            if n_elements > n_for_ngram:
                if check:
                    print(f"#ignored: {segs} [n_elements: {n_elements} > n_for_ngram: {n_for_ngram}]\n...")
                continue
            
            ## includes if and only if count_elements(p) == n_for_ngram
            if inclusive:
                pass
            else:
                if n_elements < n_for_ngram:
                    if check:
                        print(f"#ignored: {segs} [n_elements: {n_elements} < n_for_ngram: {n_for_ngram}]\n...")
                    continue
            
            ## select by extendedness
            if extended:
                if n_elements == 1 and n_gaps == 0:
                    if check:
                        print(f"#ignored: {segs} [n_elements == 1 or n_gaps == 0]\n...")
                    continue
            else: # complicated selection for segs
                if segs[0] != gap_mark or segs[-1] != gap_mark:
                    if segs not in Q and xsegs not in xQ:
                        Q.append(segs)
                        xQ.append(xsegs)
                elif segs[0] == gap_mark or segs[-1] == gap_mark:
                    if n_elements == 1:
                        Q.append(segs)
                        xQ.append(xsegs)
                    else:
                        if check:
                            print(f"#ignored: {segs} [segs[0] or segs[-1] == gap_mark]\n...")
                        continue
            ##
            if check:
                print(f"#xsegs: {xsegs}")
            if not segs in Q and not xsegs in xQ:
                Q.append(segs)
                xQ.append(xsegs)
    ##
    if check and verbose:
        print(f"#xQ [size: {len(xQ)}]: {xQ}")
    if check:
        print(f"#Q [size: {len(Q)}]: {Q}")
    return Q

##
def gen_skippy_ngrams(L: list, n_for_ngram: int, max_gap_size: int = None, extended: bool = True, inclusive: bool = True, recursively: bool = False, sep: str = " ", gap_mark: str = "…", as_list: bool = False, recursion_level: int = 0, verbose: bool = False, check: bool = False):

    """
    general generator function that can be called
    """

    ##
    if check:
        print(f"#recursion_level: {recursion_level}")
    
    ## confirm assumption
    assert n_for_ngram > 0
    
    ## filter out empty elements
    base_segs = [ seg for seg in L if len(seg) > 0 ]
    
    ##
    n_base_segs = len(base_segs)
    if n_base_segs <= n_for_ngram:
        if recursively:
            recursion_level += 1
            G = gen_skippy_ngrams(base_segs, n_for_ngram - recursion_level, max_gap_size = max_gap_size, extended = extended, inclusive = inclusive, recursively = recursively, sep = sep, gap_mark = gap_mark, as_list = as_list, recursion_level = recursion_level, verbose = verbose, check = check)
            if as_list:
                if not segs in G:
                    G = G + [ base_segs ]
                return G
            else:
                g = sep.join(base_segs)
                if not g in G:
                    G.append(base_segs)
                return [ sep.join(base_segs) for segs in G ] 
        else:
            if as_list:
                return [ base_segs ]
            else:
                return [ sep.join(base_segs) ]

    ## define divided or not
    if max_gap_size is None or n_base_segs <= max_gap_size:
        divided = False
    else:
        divided = True
    if check:
        print(f"#search divided = {divided} [len(segs) = {n_base_segs}]")
    ##
    subsegs_full = gen_ngrams(base_segs, n_base_segs, inclusive = True, sep = sep, as_list = True)
    ## prevent irrelevant elements generated wrongly
    subsegs_full = [ x for x in subsegs_full if type(x) is list or type(x) is tuple ]
    if check:
        print(f"#subsegs_full (size: {len(subsegs_full)}): {subsegs_full}")
    
    ## create subsegs_base with considration of max_gap_size
    if divided:
        subsegs_base = []
        for G in gen_ngrams(subsegs_full, max_gap_size, inclusive = True, as_list = True, check = check):
            for g in G:
                if not g in subsegs_base:
                    subsegs_base.append(g)
    else:
        subsegs_base = subsegs_full
    ##
    base_size = len(subsegs_base)
    if check:
        print(f"#subsegs_base (size: {base_size}): {subsegs_base}")
    
    ## generate a lattice of segs, with or without subdivision
    import itertools
    ##P = itertools.product(gen_source(segs)) # not work!
    #P = list(itertools.product(*gen_source(subsegs))) # list(..) is needed to retain result
    P = filter_segs(subsegs_base, n_for_ngram, max_gap_size, extended = extended, inclusive = inclusive, gap_mark = gap_mark, verbose = verbose, check = check)
    if check:
        print(f"#P0 [size: {len(P)}]: {P}")
        
    ## The following is incompatible with recursively = True
    if not inclusive:
        P = [ p for p in P if count_elements(p, gap_mark = gap_mark) == n_for_ngram ]
    if check:
        print(f"#P1 [size: {len(P)}]: {P}")

    ## regulate gaps
    Q = [ ]
    for i, p in enumerate(P):
        if check and verbose:
            print(f"#p{i}: {p}")
        
        ## simplify a series of gaps
        q = simplify_gaps(p, gap_mark = gap_mark, check = check)
        if check and verbose:
            print(f"#q1: {q}")
        
        ## remove gaps at ends
        if not extended:
            q = drop_gap_at_end(q, gap_mark = gap_mark)
        if check:
            print(f"#q2: {q}")
        Q.append(q)
    ##
    if check:
        print(f"#Q: {Q}")
    
    ## return
    if check:
        print(f"#Q: {Q}")
    if as_list:
        return Q
    else:
        return [ sep.join(x) for x in Q ]

## aliases
gen_sk_ngrams = gen_skippy_ngrams

##
def test_gen_ngrams(docs, max_n_for_ngram: int, inclusive: bool = True, as_list: bool = False, verbose: bool = False, reordered: bool = True, check: bool = False):

    ## generate skippy n-grams recursively
    import pprint as pp
    for i, doc in enumerate(docs):
        print("#==========")
        print(f"#generating normal n-grams from word {i}: {doc}")
        doc_segs = segment(doc)
        print(f"#doc_segs: {doc_segs}")
        for i in range(1, max_n_for_ngram + 1):
            print("#----------")
            print(f"#normal {i}-grams with inclusive = {inclusive} ...")
            O = gen_ngrams(doc_segs, i, inclusive = inclusive, sep = "", as_list = as_list, check = check)
            ##
            if reordered:
                O = sorted(O, key = lambda x: count_elements(x), reverse = True)
            if verbose:
                pp.pprint(O)
            else:
                print(O)

##
def test_gen_skippy_ngrams(docs, max_n_for_ngram: int, max_gap_size: int, inclusive: bool = True, extended: bool = True, as_list: bool = False, verbose: bool = False, reordered: bool = True, check: bool = False):

    ## generate skippy n-grams recursively
    import pprint as pp
    for i, doc in enumerate(docs):
        print("#==========")
        print(f"#generating skippy n-grams from word {i}: {doc}")
        doc_segs = segment(doc)
        print(f"#doc_segs: {doc_segs}")
        for i in range(1, max_n_for_ngram + 1):
            #print("#----------")
            print(f"#skippy {i}-grams [max_gap_size = {max_gap_size}, extended = {extended}, inclusive = {inclusive}]")
            O = gen_skippy_ngrams(doc_segs, i, max_gap_size, extended = extended, inclusive = inclusive, sep = "", as_list = as_list, check = check)
            ##
            if reordered:
                O = sorted(O, key = lambda x: count_elements(x), reverse = True)
            if verbose:
                pp.pprint(O)
            else:
                print(O)
    ##
    print(f"#parameters: extended = {extended}; inclusive = {inclusive}; max_n_for_ngram = {max_n_for_ngram}; max_gap_size = {max_gap_size}")

##
def main():

    """
    test code
    """
    
    t = "abcde"
    u = "abc"
    docs = [ t, u ]
    #max_seg_n = 5
    #docs = [ doc for doc in docs if len(doc) <= max_seg_n ]
    print(f"docs: {docs}")
    as_list = False
    extended = True
    inclusive = True
    max_n_for_ngram = 4
    max_gap_size = max(len(doc) for doc in docs)
    verbose = False
    check = False
    
    ## test gen_skippy_ngrams
    show_normal_ngrams = False
    if show_normal_ngrams:
        test_gen_ngrams(docs, max_n_for_ngram = max_n_for_ngram, inclusive = inclusive, as_list = as_list, verbose = verbose, check = check)
    
    ## test gen_ngrams
    test_gen_skippy_ngrams(docs, max_n_for_ngram = max_n_for_ngram, max_gap_size = max_gap_size, extended = extended, inclusive = inclusive, as_list = as_list, verbose = verbose, check = check)

##
if __name__ == "__main__":
    main()

### end of file