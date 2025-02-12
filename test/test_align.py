# Importing Dependencies
import pytest
from align import NeedlemanWunsch, read_fasta
import numpy as np

def test_nw_alignment():
    """
    TODO: Write your unit test for NW alignment
    using test_seq1.fa and test_seq2.fa by
    asserting that you have correctly filled out
    the your 3 alignment matrices.
    Use the BLOSUM62 matrix and a gap open penalty
    of -10 and a gap extension penalty of -1.
    """
    seq1, _ = read_fasta("./data/test_seq1.fa")
    seq2, _ = read_fasta("./data/test_seq2.fa")
    
    # Initalise NW algorithm
    # Using the BLOSUM62 matrix and a gap open penalty
    # of -10 and a gap extension penalty of -1.
    NW_alg = NeedlemanWunsch(sub_matrix_file = 'substitution_matrices/BLOSUM62.mat',
                             gap_open = -10,
                             gap_extend = -1
                            )
                            
    align_score, seq1_align, seq2_align = NW_alg.align(seq1, seq2)                    
    
    pass
    

def test_nw_backtrace():
    """
    TODO: Write your unit test for NW backtracing
    using test_seq3.fa and test_seq4.fa by
    asserting that the backtrace is correct.
    Use the BLOSUM62 matrix. Use a gap open
    penalty of -10 and a gap extension penalty of -1.
    """
    seq3, _ = read_fasta("./data/test_seq3.fa")
    seq4, _ = read_fasta("./data/test_seq4.fa")
    
    # Initalise NW algorithm
    # Using the BLOSUM62 matrix and a gap open penalty
    # of -10 and a gap extension penalty of -1.
    NW_alg = NeedlemanWunsch(sub_matrix_file = 'substitution_matrices/BLOSUM62.mat',
                             gap_open = -10,
                             gap_extend = -1
                            )
                            
    align_score, seq3_align, seq4_align = NW_alg.align(seq3, seq4)  
    
    # print(align_score)
    # print(seq3_align)
    # print(seq4_align)
    # raise ValueError("????")
   
    pass




