#!/usr/bin/env python3
# -*- encoding: utf-8 -*-

"""
gen2_ngrams_cy.pyx

This is a Cython-compatible Python script for generation of normal (i.e., continuous) n-grams and skippy (i.e., discontinuous) n-grams, regular or extended.

This script is a re-implementation of gen_ngrams.py and dispenses with gen_extended_skippy_ngrams(..). There are two main differences from its predecessor. First, gen_skippy_ngrams(..) generates extended skippy n-grams with option extended = True. Second, gen_skippy_ngrams(..) optionally generates inclusive n-grams by "inclusive = True" option, thereby dispensing with incremental generation of 1-gram up to n-grams afterwards.

Creation
2025/08/22 by Kow Kuroda (kow.kuroda@gmail.com) by converting gen2_ngrams.py to gen2_ngrams.pyx using encrypt.py [https://github.com/Pranjalab/pyencrypt] with renaming it to gen2_ngrams_cy.pyx for disambiguation with gen2_ngrams.py.

Modifications

"""

##
def gen_ngrams (S: list, n: int, inclusive: bool = False, sep: str = " ", as_list: bool = False, check: bool = False, recursively: bool = True):
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
        if recursively:
            G = gen_ngrams(segs, n - 1, inclusive = inclusive, sep = sep, as_list = as_list, check = check, recursively = recursively)
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
        for n in range(1, n + 1):
            for i in range(len(segs)):
                try:
                    gram = segs[i : i + n] # get an n-gram
                    if len(gram) == n:
                        R.append(gram)
                except IndexError:
                    pass
    else:
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
def gen_source(L: list, gap_mark: str = "…"):

    ## for Cython compatibility
    #return [ [x, y] for x, y in zip(L, gap_mark * len(L)) ]
    return [ (x, y) for x, y in zip(L, gap_mark * len(L)) ]

##
def gen_source_prev(L: list, gap_mark: str = "…", as_tuple: bool = False):

    if as_tuple:
        return ( [x, y] for x, y in zip(L, gap_mark * len(L)) )
    else:
        #return [ (x, y) for x, y in zip(L, gap_mark * len(L)) ]
        return [ [x, y] for x, y in zip(L, gap_mark * len(L)) ] # for Cython compatibility
##
def count_elements(segs: list, gap_mark: str = "…"):

    """
    counts the number of non-gap elements in a segment sequence and returns it
    """

    #return len(seg for seg in segs if x != gap_mark) # fails due to generator mishandling
    return len([ seg for seg in segs if seg != gap_mark ])

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
def gen_skippy_ngrams(L: list, n: int, max_gap_size: int = None, extended: bool = True, inclusive: bool = False, sep: str = " ", gap_mark: str = "…", as_list: bool = False, verbose: bool = False, check: bool = True, recursively: bool = True):

    """
    general generator function that can be called
    """

    assert n > 0

    ## filter out empty elements
    segs = [ seg for seg in L if len(seg) > 0 ]
    if len(segs) <= n:
        if recursively:
            G = gen_skippy_ngrams(segs, n - 1, max_gap_size = max_gap_size, extended = extended, inclusive = inclusive, sep = sep, gap_mark = gap_mark, as_list = as_list, verbose = verbose, check = check, recursively = recursively)
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
        subsegs_set = gen_ngrams(segs, search_span, inclusive = False, sep = sep, as_list = True)
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
        segs_size = count_elements(list(p))
        if check:
            print(f"segs_size: {segs_size}")
        if not inclusive:
            if segs_size == n:
                q = simplify_gaps(list(p), gap_mark = gap_mark, check = check)
                Q.append(q)
        else:
            if segs_size <= n:
                q = simplify_gaps(list(p), gap_mark = gap_mark, check = check)
                Q.append(q)

    ## remove gap_mark singleton
    Q = [ p for p in Q if count_elements(list(p)) > 0 ]

    ## handle 1-grams
    if not extended:
        Q = [ remove_gaps(p, gap_mark) if count_elements(list(p)) == 1 and len(p) > 1 else p for p in Q ]

    ## return
    if check:
        print(f"Q: {Q}")
    if as_list:
        return Q
    else:
        return [ sep.join(x) for x in Q ]

##
def test_gen_ngrams(docs, max_n_for_ngram: int, inclusive: bool = True, as_list: bool = False, check: bool = False, verbose: bool = False, reordered: bool = True):

    ## generate skippy n-grams recursively
    import pprint as pp
    for i, doc in enumerate(docs):
        print("------------------------------------")
        print(f"generating normal n-grams from word {i}: {doc}")
        doc_segs = segment(doc)
        print(f"doc_segs: {doc_segs}")
        for i in range(1, max_n_for_ngram + 1):
            print(f"normal {i}-grams with inclusive = {inclusive} ...")
            O = gen_ngrams(doc_segs, i, inclusive = inclusive, sep = "", as_list = as_list, check = False)
            if reordered:
                O = sorted(O, key = lambda x: count_elements(list(x)), reverse = True)
            pp.pprint(O)

##
def test_gen_skippy_ngrams(docs, max_n_for_ngram: int, max_gap_size: int, inclusive: bool = True, extended: bool = True, as_list: bool = False, check: bool = False, verbose: bool = False, reordered: bool = True):

    ## generate skippy n-grams recursively
    import pprint as pp
    for i, doc in enumerate(docs):
        print("------------------------------------")
        print(f"generating skippy n-grams from word {i}: {doc}")
        doc_segs = segment(doc)
        print(f"doc_segs: {doc_segs}")
        for i in range(1, max_n_for_ngram + 1):
            print(f"skippy {i}-grams with max_gap_size = {max_gap_size}, extended = {extended}, inclusive = {inclusive} ...")
            O = gen_skippy_ngrams(doc_segs, i, max_gap_size, extended = extended, inclusive = inclusive, sep = "", as_list = as_list, check = False)
            if reordered:
                O = sorted(O, key = lambda x: count_elements(list(x)), reverse = True)
            pp.pprint(O)

##
if __name__ == "__main__":

    s = "you are truly my true woman"
    docs = s.split()
    print(f"docs: {docs}")
    inclusive = True
    as_list = False
    check = False
    
    ## test gen_skippy_ngrams
    test_gen_ngrams(docs, max_n_for_ngram = 3, inclusive = inclusive, as_list = as_list, check = check)
    
    ## test gen_ngrams
    test_gen_skippy_ngrams(docs, max_n_for_ngram = 3, max_gap_size = 3, inclusive = inclusive, as_list = as_list, check = check)
    
### end of file
