# Importing Dependencies
import numpy as np
from typing import Tuple

# Defining class for Needleman-Wunsch Algorithm for Global pairwise alignment
class NeedlemanWunsch:
    """ Class for NeedlemanWunsch Alignment

    Parameters:
        sub_matrix_file: str
            Path/filename of substitution matrix
        gap_open: float
            Gap opening penalty
        gap_extend: float
            Gap extension penalty

    Attributes:
        seqA_align: str
            seqA alignment
        seqB_align: str
            seqB alignment
        alignment_score: float
            Score of alignment from algorithm
        gap_open: float
            Gap opening penalty
        gap_extend: float
            Gap extension penalty
    """
    def __init__(self, sub_matrix_file: str, gap_open: float, gap_extend: float):
        # Init alignment and gap matrices
        self._align_matrix = None
        self._gapA_matrix = None
        self._gapB_matrix = None

        # Init matrices for backtrace procedure
        self._back = None
        self._back_A = None
        self._back_B = None

        # Init alignment_score
        self.alignment_score = 0

        # Init empty alignment attributes
        self.seqA_align = ""
        self.seqB_align = ""

        # Init empty sequences
        self._seqA = ""
        self._seqB = ""

        # Setting gap open and gap extension penalties
        self.gap_open = gap_open
        assert gap_open < 0, "Gap opening penalty must be negative."
        self.gap_extend = gap_extend
        assert gap_extend < 0, "Gap extension penalty must be negative."

        # Generating substitution matrix
        self.sub_dict = self._read_sub_matrix(sub_matrix_file) # substitution dictionary

    def _read_sub_matrix(self, sub_matrix_file):
        """
        DO NOT MODIFY THIS METHOD! IT IS ALREADY COMPLETE!

        This function reads in a scoring matrix from any matrix like file.
        Where there is a line of the residues followed by substitution matrix.
        This file also saves the alphabet list attribute.

        Parameters:
            sub_matrix_file: str
                Name (and associated path if not in current working directory)
                of the matrix file that contains the scoring matrix.

        Returns:
            dict_sub: dict
                Substitution matrix dictionary with tuple of the two residues as
                the key and score as value e.g. {('A', 'A'): 4} or {('A', 'D'): -8}
        """
        with open(sub_matrix_file, 'r') as f:
            dict_sub = {}  # Dictionary for storing scores from sub matrix
            residue_list = []  # For storing residue list
            start = False  # trigger for reading in score values
            res_2 = 0  # used for generating substitution matrix
            # reading file line by line
            for line_num, line in enumerate(f):
                # Reading in residue list
                if '#' not in line.strip() and start is False:
                    residue_list = [k for k in line.strip().upper().split(' ') if k != '']
                    start = True
                # Generating substitution scoring dictionary
                elif start is True and res_2 < len(residue_list):
                    line = [k for k in line.strip().split(' ') if k != '']
                    # reading in line by line to create substitution dictionary
                    assert len(residue_list) == len(line), "Score line should be same length as residue list"
                    for res_1 in range(len(line)):
                        dict_sub[(residue_list[res_1], residue_list[res_2])] = float(line[res_1])
                    res_2 += 1
                elif start is True and res_2 == len(residue_list):
                    break
        return dict_sub

    def align(self, seqA: str, seqB: str) -> Tuple[float, str, str]:
        """
        TODO
        
        This function performs global sequence alignment of two strings
        using the Needleman-Wunsch Algorithm
        
        Parameters:
        	seqA: str
         		the first string to be aligned
         	seqB: str
         		the second string to be aligned with seqA
         
        Returns:
         	(alignment score, seqA alignment, seqB alignment) : Tuple[float, str, str]
         		the score and corresponding strings for the alignment of seqA and seqB
        """
        # Resetting alignment in case method is called more than once
        self.seqA_align = ""
        self.seqB_align = ""

        # Resetting alignment score in case method is called more than once
        self.alignment_score = 0

        # Initializing sequences for use in backtrace method
        self._seqA = seqA
        self._seqB = seqB
        
        # TODO: Initialize matrix private attributes for use in alignment
        # create matrices for alignment scores, gaps, and backtracing
        
        # Length of each sequence
        n = len(self._seqA)
        m = len(self._seqB)
        
        # Initialise alignment scores matrix
        # Default penalty should be -inf
        self._align_matrix = np.full((n + 1, m + 1), -float('inf'))
        self._align_matrix[0][0] = 0
        
        # Initialise gap matrices 
        self._gapA_matrix = np.full((n + 1, m + 1), -float('inf'))
        self._gapA_matrix = np.full((n + 1, m + 1), -float('inf'))
        
        # Initialise backtracing matrix
        self._back = np.full((n + 1, m + 1), -float('inf'))
        
        # Now calculate gap penalty scores
        # for fully gapped sequences 
        for i in range(0, n + 1):
            # penalty for opening gap + penalty for extending gap * length of extension
            self._gapA_matrix[i,0] = self.gap_open + self.gap_extend * i
        #  M[i,0] = gap penalty score * 
        
        for j in range(0, m + 1):
            # penalty for opening gap + penalty for extending gap * length of extension
            self._gapB_matrix[0,j] = self.gap_open + self.gap_extend * j
        # M[0,j] = 
        
        # Build Local Table
        
        # Start with M(0, j) and M(i, 0)
        
        # tuple of the two residues as the key and score as value e.g. {('A', 'A'): 4} or {('A', 'D'): -8}
        
        # for every combination of bases in seqA and seqB 
        # calculate local alignment score row by row 
        # (row by row calculation is why i/seqA indexes are the outer for loop)
        for i in range(1, n + 1):
            for j in range(1, m + 1):
        
            
            # M(i,j) = max(M(i-1, j) + 0, M(i, j-1) + 0, M(i-1, j-1) + t(i, j)
            
            # calculate local alignment score row by row
            # M(i,j) = max(0,
            #              M(i, j-1) + 0, 
            #              M(i-1, j) + 0,
            #              M(i-1,j-1) + sim(s1[i], s2[j])
            #              )
            
            # get relevant subsitition score (sim(s1[i], s2[j]))
            sub_score = self.sub_dict[(seqA[i], seqB[j])]
            
            # self._align_matrix[i-1,j-1] + sim_score
            # 
            # self._align_matrix[i-1,j] #+ gap penalty
            
            
        
        



        # return (self.alignment_score, self.seqA_align, self.seqB_align)
        # pass
        # 
        # 
        # # TODO: Implement global alignment here
        # pass      		
        		    
        return self._backtrace()

    def _backtrace(self) -> Tuple[float, str, str]:
        """
        TODO
        
        This function traces back through the back matrix created with the
        align function in order to return the final alignment score and strings.
        
        Parameters:
        	None
        
        Returns:
         	(alignment score, seqA alignment, seqB alignment) : Tuple[float, str, str]
         		the score and corresponding strings for the alignment of seqA and seqB
        """
        #pass
        
        # Length of each sequence
        n = len(self._seqA)
        m = len(self._seqB)

        return (self.alignment_score, self.seqA_align, self.seqB_align)


def read_fasta(fasta_file: str) -> Tuple[str, str]:
    """
    DO NOT MODIFY THIS FUNCTION! IT IS ALREADY COMPLETE!

    This function reads in a FASTA file and returns the associated
    string of characters (residues or nucleotides) and the header.
    This function assumes a single protein or nucleotide sequence
    per fasta file and will only read in the first sequence in the
    file if multiple are provided.

    Parameters:
        fasta_file: str
            name (and associated path if not in current working directory)
            of the Fasta file.

    Returns:
        seq: str
            String of characters from FASTA file
        header: str
            Fasta header
    """
    assert fasta_file.endswith(".fa"), "Fasta file must be a fasta file with the suffix .fa"
    with open(fasta_file) as f:
        seq = ""  # initializing sequence
        first_header = True
        for line in f:
            is_header = line.strip().startswith(">")
            # Reading in the first header
            if is_header and first_header:
                header = line.strip()  # reading in fasta header
                first_header = False
            # Reading in the sequence line by line
            elif not is_header:
                seq += line.strip().upper()  # generating full sequence
            # Breaking if more than one header is provided in the fasta file
            elif is_header and not first_header:
                break
    return seq, header
