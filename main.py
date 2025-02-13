# Import NeedlemanWunsch class and read_fasta function
from align import read_fasta, NeedlemanWunsch

def main():
    """
    This function should
    (1) Align all species to humans and print species in order of most similar to human BRD
    (2) Print all alignment scores between each species BRD2 and human BRD2
    """
    hs_seq, hs_header = read_fasta("./data/Homo_sapiens_BRD2.fa")
    gg_seq, gg_header = read_fasta("./data/Gallus_gallus_BRD2.fa")
    mm_seq, mm_header = read_fasta("./data/Mus_musculus_BRD2.fa")
    br_seq, br_header = read_fasta("./data/Balaeniceps_rex_BRD2.fa")
    tt_seq, tt_header = read_fasta("./data/Tursiops_truncatus_BRD2.fa")

    # TODO Align all species to humans and print species in order of most similar to human BRD
    # using gap opening penalty of -10 and a gap extension penalty of -1 and BLOSUM62 matrix
    
    # Initalise NW algorithm
    # Using the BLOSUM62 matrix and a gap open penalty
    # of -10 and a gap extension penalty of -1.
    NW_alg = NeedlemanWunsch(sub_matrix_file = 'substitution_matrices/BLOSUM62.mat',
                             gap_open = -10,
                             gap_extend = -1
                            )
                            
    # Calculate alignment score similarity to human BRD for each species
    gg_score, hs_gg_align, gg_hs_align = NW_alg.align(hs_seq, gg_seq)
    mm_score, hs_mm_align, mm_hs_align = NW_alg.align(hs_seq, mm_seq)
    br_score, hs_br_align, br_hs_align = NW_alg.align(hs_seq, br_seq)
    tt_score, hs_tt_align, tt_hs_align = NW_alg.align(hs_seq, tt_seq)
    
    # store alignment score results in a list
    alignment_scores = [
    ("Gallus gallus (Chicken)", gg_score),
    ("Mus musculus (Mouse)", mm_score),
    ("Balaeniceps rex (Shoebill)", br_score),
    ("Tursiops truncatus (Dolphin)", tt_score)
                       ]
                       
    # Sort species by alignment score in descending order (most similar first)
    alignment_scores.sort(key=lambda x: x[1], reverse=True)
    
    # TODO Print species ranked by similarity to human BRD2
    # TODO Print all of the alignment score between each species BRD2 and human BRD2
    print("Species ranked by similarity to Human BRD2:")
    for rank, (species, score) in enumerate(alignment_scores, start=1):
         print(f"{rank}. {species}: Alignment Score = {score}")
    
if __name__ == "__main__":
    main()
