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
    
    # 1. Check align_matrix: 
    
    # align matrix should look like this: 
    #[[  0. -inf -inf -inf]
    # [-inf   5. -11. -13.]
    # [-inf -12.   4.  -8.]
    # [-inf -12.  -1.   5.]
    # [-inf -14.  -6.   4.]]
   
    true_align_mat = np.array([
                                [0,   -float('inf'), -float('inf'), -float('inf')],
                                [-float('inf'),   5,  -11,  -13], #????13],
                                [-float('inf'), -12,    4,   -8],
                                [-float('inf'), -12,   -1,    5],
                                [-float('inf'), -14,   -6,    4]
                              ])
    
    
    assert np.allclose(true_align_mat, NW_alg._align_matrix), "Error: Estimated alignment matrix is incorrect (for seq1 and seq2)" 


    # 2. Check gapA_matrix
    # Basic structure of gapA_matrix should look like this: 
    #[[-10 -inf. -inf. -inf.]
    # [-11 -inf. -inf. -inf.]
    # [-12  ?.  ?.  ?.]
    # [-13  ?.  ?.  ?.]
    # [-14  ?.  ?.  ?.]]
    
    # Using this basic structure, and the alignment score matrix above,
    # by hand, gapA_matrix should look like this:
    #[
    # [-10 -inf. -inf. -inf.]
    # [-11 -inf. -inf. -inf.]
    # [-12,  5 - 11 = -6,  -11 -11 = -22., -13 - 11 = -24]
    # [-13  max(-6 -1, -12 - 11) = -7, max(-22 - 1, 4 - 11) = -7,   max(-24 -1, -8 - 11) = -19]
    # [-14  max(-7 -1, -14 - 11) = -8,  max(-7 -1, -1 - 11) = -8, max(-19 - 1, 5 - 11) = -6]
    # ]
    true_gapA_matrix = np.array([
                                  [-10, -float('inf'), -float('inf'), -float('inf')],
                                  [-11, -float('inf'), -float('inf'), -float('inf')],
                                  [-12, -6, -22, -24],
                                  [-13, -7, -7, -19],
                                  [-14, -8, -8, -6]
                                ])
                                
    assert np.allclose(NW_alg._gapA_matrix, true_gapA_matrix), "Error: Estimated gapA matrix is incorrect"
    
    # 3. Check gapB_matrix
    # Basic structure of gapB_matrix should look like this: 
    #[[-10 -11. -12. -13.]
    # [-inf  -inf. ?. ?.]
    # [-inf  -inf.  ?.  ?.]
    # [-inf -inf.  ?.  ?.]
    # [-inf -inf.  ?.  ?.]]
    
    # Using this basic structure, and the alignment score matrix above,
    # by hand, gapB_matrix should look like this:
    #[[-10 -11. -12. -13.]
    # [-inf  -inf. 5-11 = -6, max(-11 - 11, -6 - 1) = -7]
    # [-inf  -inf.  -12 -11 = -23,   max(4 - 11,-23 - 1) = -7]
    # [-inf -inf.  -12 - 11 = -23,  max(-1 -11 , -23 -1) = -12]
    # [-inf -inf.  -14 - 11 = -25.  max(-6-11, 25 - 1) = -17]]
    
    true_gapB_matrix = np.array([
                                  [-10, -11, -12, -13],
                                  [-float('inf'), -float('inf'), -6, -7],
                                  [-float('inf'), -float('inf'), -23, -7],
                                  [-float('inf'), -float('inf'), -23, -12],
                                  [-float('inf'), -float('inf'), -25, -17]
                               ])

    assert np.allclose(true_gapB_matrix, NW_alg._gapB_matrix), "Error: Estimated gapB matrix is incorrect"
    

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
    
    # Check alignment score is correct: 
    assert align_score == 17, 'Error in NW backtracking: alignment score was not calculated corretly'
    
    # Check that aligned sequences are correct:
    assert seq3_align == 'MAVHQLIRRP', 'Error in NW backtracking: sequence was not aligned correctly (seq3)'
    assert seq4_align == 'M---QLIRHP', 'Error in NW backtracking: sequence was not aligned correctly (seq4)'
   
    



