import numpy as np
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
'''
# Define the sequences and scoring rules
sequence_S = "AGCGA" # Rows
sequence_T = "ACGAA" # Columns
match_score = 1
mismatch_score = -1
gap_penalty = -2
'''

# TP
sequence_S = "TATGAACTA" # Rows
sequence_T = "TGA" # Columns
#ACG
#ACGGGT
match_score = 1
mismatch_score = 0
gap_penalty = 0


# -------------------------------------- Needleman-Wunsch Algorithm for global alignment --------------------------------------
# Note: The Algorithm Consistes Of 3 Stages :  [Initialization | Matrix Filling | Trace Back]
#
def print_arrow(arrow):
    if arrow == 1:
        return '↑'  # Up arrow
    elif arrow == 2:
        return '←'  # Left arrow
    elif arrow == 3:
        return '↖'  # Diagonal arrow
    else:
        return ' '
    
def needleman_wunsch(seq1, seq2):
    len_seq1 = len(seq1)
    len_seq2 = len(seq2)

    # ----------- [ Initialization ] -----------
    # We Create a matrix to store the alignment scores [ Initialization ]
    score_matrix = np.zeros((len_seq1 + 1, len_seq2 + 1))
    arrow_matrix = np.zeros((len_seq1 + 1, len_seq2 + 1))

    #print("NpZero Score_Matrix:\n",score_matrix)

    # Initialize the first row and column with gap penalties
    for i in range(len_seq1 + 1):
        score_matrix[i][0] = i * gap_penalty
    for j in range(len_seq2 + 1):
        score_matrix[0][j] = j * gap_penalty

    print("Init Matrix Now: \n",score_matrix)

    # ----------- [ Matrix Filling ] -----------
    # Fill in the matrix using dynamic programming
    for i in range(1, len_seq1 + 1):
        for j in range(1, len_seq2 + 1):
            #Match  Top Left      Example D(1,1) = Max (D(0,0),D(0,1),D(1,0))
            match = score_matrix[i - 1][j - 1] + (match_score if seq1[i - 1] == seq2[j - 1] else mismatch_score)
            #Insertion S(a,-)  <-
            insert = score_matrix[i][j - 1] + gap_penalty
            #Deletion S(-,b)  /\
            delete = score_matrix[i - 1][j] + gap_penalty

            # Score_matrix Decision
            # Example D(1,1) = Max (D(0,0),D(0,1),D(1,0))
            score_matrix[i][j] = max(match, insert, delete)

    print("Matrix After DP : \n",score_matrix)
    #print(score_matrix[0][1])
    # Traceback to find the optimal alignment
    alignment_seq1 = ""
    alignment_seq2 = ""

    i, j = len_seq1, len_seq2 #Meaning the Bot Left of the matrix

    while i > 0 and j > 0:
        current_score = score_matrix[i][j]
        # If moving to the left (gap in sequence 2), add a gap in sequence 2.
        if i > 0 and score_matrix[i - 1][j] + gap_penalty == current_score:     #/\
            alignment_seq1 = seq1[i - 1] + alignment_seq1
            alignment_seq2 = "-" + alignment_seq2
            arrow_matrix[i][j] = 1  # Up arrow
            i -= 1
        #If moving upward (gap in sequence 1), add a gap in sequence 1.
        elif j > 0 and score_matrix[i][j - 1] + gap_penalty == current_score: 
            alignment_seq1 = "-" + alignment_seq1
            alignment_seq2 = seq2[j - 1] + alignment_seq2
            arrow_matrix[i][j] = 2  # Left arrow
            j -= 1
        else: # If moving diagonally (match or mismatch), add the corresponding characters to both sequences.
            alignment_seq1 = seq1[i - 1] + alignment_seq1
            alignment_seq2 = seq2[j - 1] + alignment_seq2
            arrow_matrix[i][j] = 3  # Diagonal arrow

            i -= 1
            j -= 1
        
    print("Arrow Matrix:")
    for row in arrow_matrix:
        print(' '.join([print_arrow(arrow) for arrow in row]))
    return alignment_seq1, alignment_seq2
    

# -------------------------------------- Smith-Waterman Algorithm for local alignment --------------------------------------
#Note:  Negative Values Become 0 
def smith_waterman(seq1, seq2):
    len_seq1 = len(seq1)
    len_seq2 = len(seq2)

    # Create a matrix to store the alignment scores
    score_matrix = np.zeros((len_seq1 + 1, len_seq2 + 1))
    arrow_matrix = np.zeros((len_seq1 + 1, len_seq2 + 1))

    print("Init Matrix Now: \n",score_matrix)

    # Unlike Needlemen, The Smith Initializes the matrix with 0's
    max_score = 0
    max_i, max_j = 0, 0

    for i in range(1, len_seq1 + 1):
        for j in range(1, len_seq2 + 1):
            match = score_matrix[i - 1][j - 1] + (match_score if seq1[i - 1] == seq2[j - 1] else mismatch_score)
            insert = score_matrix[i][j - 1] + gap_penalty
            delete = score_matrix[i - 1][j] + gap_penalty
            score_matrix[i][j] = max(0, match, insert, delete) # This insures there are no negative numbers =)

            # T his saves our Max Value for later, we need it in Back-Tracing 
            if score_matrix[i][j] > max_score:
                max_score = score_matrix[i][j]
                max_i, max_j = i, j

    # Traceback to find the optimal alignment
    alignment_seq1 = ""
    alignment_seq2 = ""
    i, j = max_i, max_j
    while i > 0 and j > 0 and score_matrix[i][j] > 0:
        current_score = score_matrix[i][j]
        if current_score == 0:
            break
        if i > 0 and score_matrix[i - 1][j] + gap_penalty == current_score: #/\
            alignment_seq1 = seq1[i - 1] + alignment_seq1
            alignment_seq2 = "-" + alignment_seq2
            arrow_matrix[i][j] = 1  # Up arrow
            i -= 1  
        elif j > 0 and score_matrix[i][j - 1] + gap_penalty == current_score: # <-
            alignment_seq1 = "-" + alignment_seq1
            alignment_seq2 = seq2[j - 1] + alignment_seq2
            
            arrow_matrix[i][j] = 2  # Left arrow
            
            j -= 1
        else:
            alignment_seq1 = seq1[i - 1] + alignment_seq1
            alignment_seq2 = seq2[j - 1] + alignment_seq2

            arrow_matrix[i][j] = 3  # Diagonal arrow
            i -= 1
            j -= 1
    print("Matrix After DP : \n",score_matrix)

    print("Arrow Matrix:")
    for row in arrow_matrix:
        print(' '.join([print_arrow(arrow) for arrow in row]))
    return alignment_seq1, alignment_seq2


#------------------------------------------------------------------------------------------------------------------

print(" --------------- Needleman-Wunsch (Global Alignment) --------------- ")
global_alignment_seq1, global_alignment_seq2 = needleman_wunsch(sequence_S, sequence_T)
print("Result After TraceBack:")
print(global_alignment_seq1)
print(global_alignment_seq2)

print("  --------------- Smith-Waterman (Local Alignment) --------------- ")
local_alignment_seq1, local_alignment_seq2 = smith_waterman(sequence_S, sequence_T)
print("Result After TraceBack:")
print(local_alignment_seq1)
print(local_alignment_seq2)


# ---------------------------- Bio ----------------------------
print("  --------------- Bio (Global Alignment) & (Local Alignment) ---------------")
# ----------  Global ----------
# Using Bio.pairwise2 for sequence alignment
alignments = pairwise2.align.globalxx(sequence_S, sequence_T)
for alignment in alignments:
    print(f"Alignment: {alignment[0]}\nAlignment: {alignment[1]}\n")

# ----------  Local ----------

def smith_waterman_bio(sequence_S, sequence_T, match_score=2, mismatch_score=-1, gap_open_penalty=-0.5, gap_extend_penalty=-0.1):
    alignments = pairwise2.align.localms(sequence_S, sequence_T, match_score, mismatch_score, gap_open_penalty, gap_extend_penalty)
    for alignment in alignments:
        print(format_alignment(*alignment))

# Example usage
sequence_S = "ACGGGT"
sequence_T = "ACG"

print("Smith-Waterman BioPython:")
smith_waterman_bio(sequence_S, sequence_T)