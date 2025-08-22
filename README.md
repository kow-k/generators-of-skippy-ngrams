# Generators of (extended) (skippy) n-grams from words or sentences

This Jupyter notebook demonstrates how to use "gen2_ngrams.py" (or "gen2_ngrams_cy.pyx") I developed to enhance the usability of its predecessor "gen_ngrams.py".

There are two main differences from its predecessor. First, gen_skippy_ngrams(..) generates extended skippy n-grams with "extended = True" option. Second, gen_skippy_ngrams(..) cann generate inclusive n-grams, thereby dispensing with incremental generation of n-grams from 1-grams.

## Limitations
Availablity of Cython-enhancement is limited. Apple Silicons like M1 and M2 (M3 is not tested yet) will not accept it, though only Python 3.10 running on M1 does.

## Files

1. [gen_ngram_runner.ipynb (Jupyter notebook)](gen_ngram_runner.ipynb)
2. [gen2_ngrams.py (Python file)](gen2_ngrams.py)
3. [gen2_ngrams_cy.pyx (Cython file)](gen2_ngrams_cy.pyx)