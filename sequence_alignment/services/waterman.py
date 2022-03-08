import pandas as pd
from IPython.core.display import HTML,display



from numpy import full

#build an array of zeroes
seq1 = "CTCCCGCT"
seq2 = "TATCC"

n_rows = len("-"+seq1)
n_columns = len("-"+seq2)

scoring_array = full([n_rows,n_columns],0)
#print("Scoring array:\n",scoring_array)

traceback_array = full([n_rows,n_columns],"-")
#print("Traceback array:\n",traceback_array)

def pretty_table_from_array(data_array, row_labels,col_labels):
    """Show an HTML table from a 2d numpy array"""
    df = pd.DataFrame(data_array,index=row_labels,columns=col_labels)
    return df

row_labels = [label for label in "-"+seq1]
column_labels = [label for label in "-"+seq2]

print("Scoring array:")
display(pretty_table_from_array(scoring_array,row_labels,column_labels))
print("Traceback array:")
display(pretty_table_from_array(traceback_array,row_labels,column_labels))

# scoring 
count = 0
for row_index in range(n_rows):
    for col_index in range(n_columns):    
        scoring_array[row_index,col_index] = count
        count += 0
        
display(pretty_table_from_array(scoring_array,row_labels,column_labels))

#traceback
up_arrow = "\u2191"
right_arrow = "\u2192"
down_arrow = "\u2193"
left_arrow = "\u2190"
down_right_arrow = "\u2198"
up_left_arrow = "\u2196"

print("Up arrow",up_arrow)
print("Left arrow",left_arrow)
print("Up Left arrow",up_left_arrow)

#build an array of zeroes
n_rows = len(seq1) + 1 #need an extra row up top
n_columns = len(seq2) + 1 #need an extra column on the left
row_labels = [label for label in "-"+seq1]
column_labels = [label for label in "-"+seq2]


scoring_array = full([n_rows,n_columns],0)
traceback_array = full([n_rows,n_columns],"-")


#Define Unicode arrows we'll use in the traceback array
up_arrow = "\u2191"
right_arrow = "\u2192"
down_arrow = "\u2193"
left_arrow = "\u2190"
down_right_arrow = "\u2198"
up_left_arrow = "\u2196"

arrow = "-"
gap_penalty = -1
match_bonus = 1
mismatch_penalty = -1
#iterate over columns first because we want to do 
# all the columns for row 1 before row 2
for row in range(n_rows):
    for col in range(n_columns):        
        if row == 0 and col == 0:
            #We're in the upper right corner
            score = 0
            arrow = "-"
        elif row == 0:
            #We're on the first row
            #but NOT in the corner
            
            #Look up the score of the previous cell (to the left) in the score array\
            previous_score = scoring_array[row,col - 1]
            # add the gap penalty to it's score
            score = previous_score + gap_penalty
            arrow = left_arrow
        elif col == 0:
            #We're on the first column but not in the first row
            previous_score = scoring_array[row -1,col]
            score = previous_score + gap_penalty
            arrow = up_arrow
        else: 
            #We're in a 'middle' cell of the alignment
            
            #Calculate the scores for coming from above,
            #from the left, (representing an insertion into seq1)
            cell_to_the_left = scoring_array[row,col-1]
            from_left_score = cell_to_the_left + gap_penalty
             
            #or from above (representing an insertion into seq2)
            above_cell = scoring_array[row-1,col]
            from_above_score = above_cell + gap_penalty
            
            #diagonal cell, representing a substitution (e.g. A --> T)
            diagonal_left_cell = scoring_array[row-1,col-1]
            
            #NOTE: since the table has an extra row and column (the blank ones), 
            #when indexing back to the sequence we want row -1 and col - 1.
            #since row 1 represents character 0 of the sequence.
            if seq1[row-1] == seq2[col-1]:
                diagonal_left_cell_score = diagonal_left_cell + match_bonus
            else:
                diagonal_left_cell_score = diagonal_left_cell + mismatch_penalty
            
            score = max([from_left_score,from_above_score,diagonal_left_cell_score]) 
            #take the max
            
            #make note of which cell was the max in the traceback array 
            #using Unicode arrows
            if score == from_left_score:
                arrow = left_arrow
            elif score == from_above_score:
                arrow = up_arrow
            elif score == diagonal_left_cell_score:
                arrow = up_left_arrow
                
        traceback_array[row,col]=arrow    
        scoring_array[row,col] = score
        
display(pretty_table_from_array(scoring_array,row_labels,column_labels))  
display(pretty_table_from_array(traceback_array,row_labels,column_labels))

def traceback_alignment(traceback_array,seq1,seq2,up_arrow = "\u2191" ,\
                        left_arrow="\u2190",up_left_arrow="\u2196",stop="-"):
    """Align seq1 and seq2 using the traceback matrix and return as two strings
    
    traceback_array -- a numpy array with arrow characters indicating the direction from 
    which the best path to a given alignment position originated
    
    seq1 - a sequence represented as a string
    seq2 - a sequence represented as a string
    up_arrow - the unicode used for the up arrows (there are several arrow symbols in Unicode)
    left_arrow - the unicode used for the left arrows 
    up_left_arrow - the unicode used for the diagonal arrows
    stop - the symbol used in the upper left to indicate the end of the alignment
    """

    n_rows = len(seq1) + 1 #need an extra row up top
    n_columns = len(seq2) + 1 #need an extra row up top
    
    row = len(seq1)
    col = len(seq2)
    arrow = traceback_array[row,col]
    aligned_seq1 = ""
    aligned_seq2 = ""
    alignment_indicator = ""
    while arrow is not "-":
            print("Currently on row:",row)
            print("Currently on col:",col)
            arrow = traceback_array[row,col]
            print("Arrow:",arrow)
            
            if arrow == up_arrow: 
                print("insert indel into top sequence")
                #We want to add the new indel onto the left 
                #side of the growing aligned sequence
                aligned_seq2 = "-"+aligned_seq2 
                aligned_seq1 = seq1[row-1] + aligned_seq1
                alignment_indicator = " "+alignment_indicator
                row -=1
                            
            elif arrow == up_left_arrow:
                print("match or mismatch")
                #Note that we look up the row-1 and col-1 indexes
                #because there is an extra "-" character at the
                #start of each sequence
                seq1_character = seq1[row-1]
                seq2_character = seq2[col-1]
                aligned_seq1 = seq1[row-1] + aligned_seq1
                aligned_seq2 = seq2[col-1] + aligned_seq2
                if seq1_character == seq2_character:
                    alignment_indicator = "|"+alignment_indicator
                else:
                    alignment_indicator = " "+alignment_indicator
                row -=1
                col -=1
                
            elif arrow == left_arrow:
                print("Insert indel into left sequence")
                aligned_seq1 = "-"+aligned_seq1
                aligned_seq2 = seq2[col-1] + aligned_seq2
                alignment_indicator = " "+alignment_indicator
                col -=1
                
            elif arrow == stop:
                break
            else:
                raise ValueError(f"Traceback array entry at {row},{col}: {arrow} is not recognized as an up arrow ({up_arrow}),left_arrow ({left_arrow}), up_left_arrow ({up_left_arrow}), or a stop ({stop}).")
            #print(traceback_array,-row,-col,traceback_array[-row,-col])
            print(aligned_seq1)
            print(alignment_indicator)
            print(aligned_seq2)
            
    return aligned_seq1,aligned_seq2
traceback_alignment(traceback_array,seq1,seq2)


# //////////////////////////////////////

#Build a dict to assign each nucleotide one row or column
#index in the table
nucleotides = "AGCT"

#Step through each nucleotide and give it a row and column index
#using a dictionary with keys = nucleotides and values = indices
nucleotide_indices = {nucleotide:i for i,nucleotide in enumerate(nucleotides)}

#Set up scores 
match_score = 1
#We want separate scores for substitutions that are
#transitions or transversions
transversion_score = -2
transition_score = -1

# Set up a scoring_matrix for each possible substitution
scoring_matrix = full([len(nucleotides),len(nucleotides)],transition_score)

#Fill in the scoring matrix based on whether the new vs. old nucleotide are in the 
#same chemical class (e.g. both purines)
chemical_class = {"A":"Purine","T":"Pyrimidine","C":"Pyrimidine","G":"Purine"}
for nt1 in nucleotides:
    for nt2 in nucleotides:
        #Look up which row/column the 
        #nucleotides are in
        nt1_index = nucleotide_indices[nt1]
        nt2_index = nucleotide_indices[nt2]
        if nt1 == nt2:
            #The nucleotides match
            scoring_matrix[nt1_index][nt2_index] = match_score
            #We can skip further analysis of this pair...
            #We alredy know they match
            continue
        
        nt1_chemical_class = chemical_class[nt1]
        nt2_chemical_class = chemical_class[nt2]
        
        if nt1_chemical_class == nt2_chemical_class:
            #The nucleotides are both pyrimidines or
            #both purines so this is a transition
            scoring_matrix[nt1_index][nt2_index] = transition_score
        else:
            #They are in different chemical classes,
            #so this change is a transversion
            scoring_matrix[nt1_index][nt2_index] = transversion_score
            

#Show the scoring matrix
display(pretty_table_from_array(scoring_matrix,\
        row_labels =[n for n in nucleotides],\
        col_labels = [n for n in nucleotides]))    


# ///////////////////////////////////////////////////////////////////////////////

def score_match(nt1,nt2,scoring_matrix,\
  scoring_matrix_indices={'A': 0, 'G': 1, 'C': 2, 'T': 3}):
    """Return the score for a substitution between nt1 and nt2 based on the scoring matrix
    nt1 -- a string representing the first nucleotide 
    nt2 -- a string representing the second nucleotide
    scoring_matrix -- an N x N numpy array, where N is
      the number of nucleotides (so usually 4x4)
    scoring_matrix_indices -- a dict mapping rows and columns
      of the scoring array to nucleotides 
    
    """
    return scoring_matrix[scoring_matrix_indices[nt1],scoring_matrix_indices[nt2]]


# ////////////////////////////////////////////////////////////////////////////////
AG_score = score_match("A","G",scoring_matrix,nucleotide_indices)
AT_score = score_match("A","T",scoring_matrix,nucleotide_indices)
# print(f"A --> G score:{AG_score}")
# print(f"A --> T score:{AT_score}")

# ////////////////////////////////////////////////////////////////////////////////
def needleman_wunsch(seq1,seq2, scoring_matrix,\
  scoring_matrix_indices={"A":0,"G":0,"G":0,"C":0},\
  scoring_function=score_match, gap_penalty=-1):
    """Perform Needleman Wunsch global alignment on two sequences
    seq1 -- a sequence as a string
    seq2 -- a sequence as a string
    gap_function -- a function that takes no parameters and returns the score for a gap
    scoring_function -- a function that takes two nucleotides and returns a score
    
    """
    #build an array of zeroes 
    n_rows = len(seq1) + 1 #need an extra row up top
    n_columns = len(seq2) + 1 #need an extra column on the left
    scoring_array = full([n_rows,n_columns],0)
    traceback_array = full([n_rows,n_columns],"-")


    #Define Unicode arrows we'll use in the traceback array
    up_arrow = "\u2191"
    right_arrow = "\u2192"
    down_arrow = "\u2193"
    left_arrow = "\u2190"
    down_right_arrow = "\u2198"
    up_left_arrow = "\u2196"

    arrow = "-"
    
    #iterate over columns first because we want to do 
    # all the columns for row 1 before row 2
    for row in range(n_rows):
        for col in range(n_columns):  
            if row == 0 and col == 0:
                #We're in the upper right corner
                score = 0
            elif row == 0:
                #We're on the first row
                #but NOT in the corner

                #Look up the score of the previous cell (to the left) in the score array\
                previous_score = scoring_array[row,col - 1]
                # add the gap penalty to it's score
                score = previous_score + gap_penalty
            elif col == 0:
                #We're on the first column but not in the first row
                previous_score = scoring_array[row -1,col]
                score = previous_score + gap_penalty
            else: 
                #We're in a 'middle' cell of the alignment

                #Calculate the scores for coming from above,
                #from the left, (representing an insertion into seq1)
                cell_to_the_left = scoring_array[row,col-1]
                from_left_score = cell_to_the_left + gap_penalty

                #or from above (representing an insertion into seq2)
                above_cell = scoring_array[row-1,col]
                from_above_score = above_cell + gap_penalty

                #diagonal cell, representing a substitution (e.g. A --> T)
               
                diagonal_left_cell = scoring_array[row-1,col-1]

                #Since the table has an extra row and column (the blank ones), 
                #when indexing back to the sequence we want row -1 and col - 1.
                #since row 1 represents character 0 of the sequence.
                curr_nt_seq1 = seq1[row-1]
                curr_nt_seq2 = seq2[col-1]
                
                #the scoring matrix will tell us the score for matches,
                #transitions and transversions
                diagonal_left_cell_score = diagonal_left_cell + \
                  score_match(curr_nt_seq1,curr_nt_seq2,scoring_matrix)
                score = max([from_left_score,from_above_score,diagonal_left_cell_score]) 
                #take the max
                #make note of which cell was the max in the traceback array 
                #using Unicode arrows
                if score == from_left_score:
                    arrow = left_arrow
                elif score == from_above_score:
                    arrow = up_arrow
                elif score == diagonal_left_cell_score:
                    arrow = up_left_arrow
            
            traceback_array[row,col]=arrow    
            scoring_array[row,col] = score
    return scoring_array,traceback_array
        
        
scoring_array,traceback_array = needleman_wunsch(seq1,seq2,scoring_matrix)
display(pretty_table_from_array(scoring_array,row_labels,column_labels))   
display(pretty_table_from_array(traceback_array,row_labels,column_labels)) 

# /////////////////////////////////////////////////////////////////////////////////////////////////////////////////

