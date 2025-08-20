"""
gen2_ngrams.py

This is a Python script for generation of n-grams and skippy n-grams, regular or extended. This is a re-implementation of gen_ngrams.py and dispenses with gen_extended_skippy_ngrams(..). Now, gen_skippy_ngrams(..) generates extended skippy n-grams by setting extended = True.

Created on 2020/08/20 by Kow Kuroda (kow.kuroda@gmail.com)

Modifications

"""

##
def segment(t: str, pattern: str = r"", check: bool = False):
    import re
    return [ x for x in re.split(pattern, t) if len(x) > 0 ]

##
def gen_source(L: list, gap_mark: str = "…"):
    return ((x, y) for x, y in zip(L, gap_mark * len(L)))
    
##
def count_elements(L: list, gap_mark: str = "…"):
    #return len(x for x in L if x != gap_mark) # fails due to generator mishandling
    return len([ x for x in L if x != gap_mark ])

##
def simplify_gaps(segs: list, gap_mark: str, check: bool = False):

    """
    simplify repeated gap_marks in row
    """

    ##
    if check:
        print(f"input: {segs}")
    ## main
    R = [ ]
    for i, seg in enumerate(segs):
        if check:
            print(f"seg{i}: {seg}")
        ## main
        if i == 0:
            R.append(seg)
        else:
            if seg == gap_mark and R[-1] == gap_mark:
                pass
            else:
                R.append(seg)
                if check:
                    print(f"added seg{i}: {seg}")
    ##
    if check:
        print(f"R: {R}")
    return R

##
def regularize_gaps(segs: list, gap_mark: str, restrictive: bool = False, check: bool = False):

    """
    regularize occurrences of gap_marks due to standard (= unextended) definition
    """
    
    if check:
        print(f"input: {segs}")
    ##
    R = [ ]
    for i, seg in enumerate(segs):
        if i == 0:
            if seg != gap_mark:
                R.append(seg)
        elif i == (len(segs) - 1):
            if seg != gap_mark:
                R.append(seg)
        else:
            try:
                if R[-1] != gap_mark:
                    R.append(seg)
            except IndexError:
                if seg != gap_mark:
                    R.append(seg)
    
    ## remove remaining gap_marks at ends
    if gap_mark in R:
        if R[0] == gap_mark and R[-1] == gap_mark:
            R = R[1:-1]
        elif R[0] == gap_mark:
            R = R[1:]
        elif R[-1] == gap_mark:
            R = R[:-1]
    ##
    return R

##
def classify_segs(Q: list, n: int, inclusive: bool, gap_mark: str, restrictive: bool = True, check: bool = False):
    """
    sort out segments according to condition for n
    """
    R = [ ]
    for r in Q:
        n_elements = count_elements(r, gap_mark = gap_mark)
        if n_elements == 0 or n_elements > n:
            if check:
                print(f"ignored {r} [{n} for n_elements: {n_elements}]")
        else:
            if n_elements == n:
                if not restrictive or not r in R:
                    R.append(r)
                    if check:
                        print(f"added {r} [{n} for n_elements: {n_elements}]")
            else:
                if inclusive:
                    if not restrictive or r not in R:
                        R.append(r)
                        if check:
                            print(f"added {r} [{n} for n_elements: {n_elements}]")
                else:
                    if check:
                        print(f"ignored {r} [{n} for n_elements: {n_elements}]")
    ##
    if check:
        print(f"R: {R}")
    ##
    return R

##
def gen_skippy_ngrams_with_fixed_size(segs: list, n: int, extended: bool, gap_mark: str = "…", inclusive: bool = False, restrictive: bool = True, check: bool = True): 
    
    """
    core processing of skippy n-gram without consideration of max_gap_size
    """

    import itertools
    #P = itertools.product(gen_source(segs)) # not work!
    P = list(itertools.product(*gen_source(segs))) # result vanishes unless list(..) is applied
    if check:
        print(f"P: {P}")
    
    ## extended or regular
    if extended:
        Q = [ simplify_gaps(p, gap_mark = gap_mark, check = check) for p in P ]
    else:
        Q = [ regularize_gaps(p, gap_mark = gap_mark, check = check) for p in P ]
    
    ## classify by n_elements
    R = classify_segs(Q, n = n, inclusive = inclusive, restrictive = restrictive, gap_mark = gap_mark, check = check)
    
    ##
    return R


## skippy n-gram generator, both extended or not
def gen_skippy_ngrams(L: list, n: int, extended: bool, inclusive: bool = False, restrictive: bool = True, max_gap_size: int = None, sep: str = " ", gap_mark: str = "…", as_list: bool = False, check: bool = True):

    """general function to be accessed"""
    
    ## filter out empty elements
    segs = [ seg for seg in L if len(seg) > 0 ]
    if len(segs) <= n:
        if as_list:
            return [segs]
        else:
            return [sep.join(segs)]
    
    ## process subsegments identified by max_gap_size
    if max_gap_size is None or len(segs) <= max_gap_size:
        segs_pool = [segs]
    else:
        segs_pool = [ ]
        d = len(segs) - max_gap_size
        if extended:
            for i in range(d + 1):
                if i == 0:
                    segs_pool.append(segs[i: i + max_gap_size] + [gap_mark])
                elif i == d:
                    segs_pool.append([gap_mark] + segs[i: i + max_gap_size])
                else:
                    segs_pool.append([gap_mark] + segs[i: i + max_gap_size] + [gap_mark])
        else:
            for i in range(d + 1):
                segs_pool.append(segs[i: i + max_gap_size])
    ##
    if check:
        print(f"segs_pool: {segs_pool}")
    ##
    R = [ ]
    for segs in segs_pool:
        ##
        G = gen_skippy_ngrams_with_fixed_size(segs, n, extended = extended, inclusive = inclusive, restrictive = restrictive, gap_mark = gap_mark, check = check)
        if check:
            print(f"G: {G}")
        ##
        O = [ segs for segs in G if not segs in R ]
        if check:
            print(f"O: {O}")
        ##
        R.extend(O)
    ##
    if as_list:
        return R
    else:
        return [ sep.join(x) for x in R ]

## normal, continus n-gram generator
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
    ##
    if len(segs) <= n:
        if as_list:
            return [ segs ]
        else:
            return [ sep.join(segs) ]
    
    ## main
    R = [ ]
    for i, seg in enumerate(segs):
        try:
            gram = segs[ i : i + n] # get an n-gram
            if len(gram) == n: # check its length
                R.append(gram)
        except IndexError:
            pass
    ##
    if as_list:
        return R
    else:
        return [ sep.join(r) for r in R ]

### end of file