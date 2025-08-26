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
        if check:
            if len(R) > R_size:
                print(f"added seg{i}: {seg}")
        R_size = len(R) # update

    ##
    if check:
        print(f"R: {R}")
    return R

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
def remove_gaps(segs: list, gap_mark: str):

    """simply removes all gaps in a segment list and returns the result"""

    return [ seg for seg in segs if seg != gap_mark ]

##
def gen_ngrams (S: list, n_for_ngram: int, inclusive: bool = False, recursively: bool = True, sep: str = " ", as_list: bool = False, check: bool = False):
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
            G = gen_ngrams(list(segs), n_for_ngram - 1, inclusive = inclusive, recursively = recursively, sep = sep, as_list = as_list, check = check)
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
def gen_skippy_ngrams(L: list, n_for_ngram: int, base_n: int, max_gap_size: int = None, extended: bool = True, inclusive: bool = False, sep: str = " ", gap_mark: str = "…", as_list: bool = False, recursively: bool = True, verbose: bool = False, check: bool = False):

    """
    general generator function that can be called
    """

    assert n_for_ngram > 0
    
    ## filter out empty elements
    segs = [ seg for seg in L if len(seg) > 0 ]
    if len(segs) <= n_for_ngram:
        if recursively:
            G = gen_skippy_ngrams(list(segs), n_for_ngram - 1, base_n = base_n, max_gap_size = max_gap_size, extended = extended, inclusive = inclusive, sep = sep, gap_mark = gap_mark, as_list = as_list, recursively = recursively, verbose = verbose, check = check)
            if as_list:
                if not segs in G:
                    G = G + [segs]
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

    ## set search_span
    if max_gap_size is None or len(segs) <= max_gap_size:
        divided = False
        search_span = len(segs) + 1
    else:
        divided = True
        search_span = len(segs) - max_gap_size
    if check:
        print(f"#search divided = {divided} with search_span = {search_span} [len(segs) = {len(segs)}]")
    ##
    subsegs_full = gen_ngrams(list(segs), len(segs), inclusive = True, sep = sep, as_list = True)
    print(f"#subsegs_full: {subsegs_full}")
    
    ## create subsegs_base with considration of max_gap_size
    if divided:
        subsegs_base = []
        for G in gen_ngrams(list(subsegs_full), max_gap_size, inclusive = True, as_list = True, check = check):
            for g in G:
                if not g in subsegs_base:
                    subsegs_base.append(g)
    else:
        subsegs_base = subsegs_full
    base_size = len(subsegs_base)
    if check:
        print(f"#subsegs_base [size: {base_size}]: {subsegs_base}")
    
    ## generate a lattice of segs, with or without subdivision
    import itertools
    ##P = itertools.product(gen_source(segs)) # not work!
    #P = list(itertools.product(*gen_source(subsegs))) # list(..) is needed to retain result
    P = [ ]; xP = [ ] # checker of iso-forms
    for i, subsegs in enumerate(subsegs_base):
        if check:
            print(f"#{i} subsegs0: {subsegs}")
        ##
        for segs in list(itertools.product(*gen_source(subsegs))):
            if check:
                print(f"#segs: {segs}")
            
            ## filter by n_for_ngram
            if not count_elements(list(segs), gap_mark) <= base_n: # n_for_ngram:
                continue
            
            ## filter by inclusiveness
            if inclusive:
                if not extended:
                    pass
            else:
                if not extended:            
                    if segs[0] == gap_mark or segs[-1] == gap_mark:
                        continue
            ## handle gaps
            xsegs = simplify_gaps(list(segs), gap_mark = gap_mark, check = check)
            if check:
                print(f"#xsegs: {xsegs}")
            if not segs in P and not xsegs in xP:
                P.append(segs)
                xP.append(xsegs)
    if check:
        print(f"#P0: {P}")
    
    ## The following is incompatible with recursively = True
    if not inclusive:
        P = [ p for p in P if count_elements(list(p), gap_mark = gap_mark) == base_n ]
    if check:
        print(f"#P1: {P}")

    ## regulate gaps
    Q = [ ]
    for i, p in enumerate(P):
        if check:
            print(f"#p{i}: {p}")
        
        ## simplify a series of gaps
        q = simplify_gaps(list(p), gap_mark = gap_mark, check = check)
        if check:
            print(f"#q1: {q}")
        
        ## remove gaps at ends
        if not extended:
            q = drop_gap_at_end(list(q), gap_mark = gap_mark)
        if check:
            print(f"#q2: {q}")
        Q.append(q)
    if check:
        print(f"#Q: {Q}")
    
    ## remove gap_mark singleton
    Q = [ p for p in Q if count_elements(list(p)) > 0 ]

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
            print("#----------")
            print(f"#skippy {i}-grams with max_gap_size = {max_gap_size}, extended = {extended}, inclusive = {inclusive} ...")
            O = gen_skippy_ngrams(doc_segs, i, i, max_gap_size, extended = extended, inclusive = inclusive, sep = "", as_list = as_list, check = check)
            ##
            if reordered:
                O = sorted(O, key = lambda x: count_elements(x), reverse = True)
            if verbose:
                pp.pprint(O)
            else:
                print(O)
    ##
    print(f"#parameters: max_n_for_ngram = {max_n_for_ngram}; max_gap_size = {max_gap_size}; extended = {extended}; inclusive = {inclusive};")

##
if __name__ == "__main__":

    s = "you are my true woman"
    docs = s.split()
    #docs = [ d for d in docs if len(d) < 4 ]
    print(f"docs: {docs}")
    max_n_for_ngram = 3
    max_gap_size = 4
    extended = True
    inclusive = False
    as_list = False
    check = False
    
    ## test gen_skippy_ngrams
    show_norma_ngrams = False
    if show_norma_ngrams:
        test_gen_ngrams(docs, max_n_for_ngram = max_n_for_ngram, inclusive = inclusive, as_list = as_list, check = check)
    
    ## test gen_ngrams
    test_gen_skippy_ngrams(docs, max_n_for_ngram = max_n_for_ngram, max_gap_size = max_gap_size, extended = extended, inclusive = inclusive, as_list = as_list, check = False)
    
### end of file