# Generators of (extended) (skippy) n-grams from words or sentences

## Overview

This repository provides "gen2_ngrams.py" (or "gen2_ngrams_cy.pyx") I developed to enhance the usability of its predecessor "gen_ngrams.py", which is available at my repositories like [LDA-mixed-terms](https://github.com/kow-k/LDA-mixed-terms), [self-supervised-word-classification](https://github.com/kow-k/self-supervised-word-classification), [HDP-spell-sound-typology](https://github.com/kow-k/HDP-spell-sound-analyzer) and some others.

There are two main differences from its predecessor, "gen_ngrams.py". First, functinality is integrated and gen_skippy_ngrams(..) now generates extended skippy n-grams with "extended = True" option.

Second, gen_skippy_ngrams(..) now can generate "inclusive" n-grams, thereby dispensing with incremental generation of n-grams from 1-grams afterwards.

## Limitations
Availablity of Cython-enhancement is limited. Apple Silicons like M1 and M2 (M3 is not tested yet) will not accept it, though only Python 3.10 running on M1 does.

## Files

1. [gen_ngram_runner.ipynb (Jupyter notebook)](gen_ngrams-runner.ipynb) is a Jupyter notebook that demonstrates how to use "gen2_ngrams.py" (or "gen2_ngrams_cy.pyx").

2. [gen2_ngrams.py (Python file)](gen2_ngrams.py) is a Python script/module that can be imporeted from a Python program.

3. [gen2_ngrams_cy.pyx (Cython file)](gen2_ngrams_cy.pyx) is a Cython script/module that can be imporeted from a Python program.
