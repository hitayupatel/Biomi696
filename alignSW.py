#! /usr/bin/env python

gap_penalty = -2
match = 1
mismatch = -1

def match_score(alpha, beta):
    '''A function for determining the score between any two bases in alignment'''
    if alpha == beta:
        return match
    elif alpha == '-' or beta == '-':
        return gap_penalty
    else:
        return mismatch


def smith_waternman(sequence1, sequence2):
    '''Main function having Needleman Wunch algorithm'''
    l1 = len(sequence1)  
    l2 = len(sequence2)

    score_matrix = [ [ 0 for i in range(l1+1) ] for j in range(l2+1) ]
        
    max_score = 0
    max_row_number = 0
    max_column_number = 0
    
    #fill the score matrix
    for i in range(1, l2 + 1):
        for j in range(1, l1 + 1):
            match = score_matrix[i - 1][j - 1] + match_score(sequence1[j-1], sequence2[i-1])
            delete = score_matrix[i - 1][j] + gap_penalty
            insert = score_matrix[i][j - 1] + gap_penalty
            if max(match, delete, insert)>=0:
                score_matrix[i][j] = max(match, delete, insert)
                if score_matrix[i][j] > max_score:
                    max_score = score_matrix[i][j]
                    max_row_number = i
                    max_column_number = j
            else:
                score_matrix[i][j] = 0

    #variables to store alignment
    alignment1 = ""
    alignment2 = ""
    
    #traceback
    l2 = max_column_number
    l1 = max_row_number
    print("l2",l2,"l1",l1)
    while l2 > 0 and l1 > 0: 
        current_score = score_matrix[l2][l1]
        if current_score == 0:
		    break
        diagonal_score = score_matrix[l2-1][l1-1]
        up_score = score_matrix[l2][l1-1]
        left_score = score_matrix[l2-1][l1]
        
        #keep track of coming score and calculate new score
        if current_score == diagonal_score + match_score(sequence1[l1-1], sequence2[l2-1]):
            alignment1 += sequence1[l1-1]
            alignment2 += sequence2[l2-1]
            l2 -= 1
            l1 -= 1
        elif current_score == up_score + gap_penalty:
            alignment1 += sequence1[l1-1]
            alignment2 += '-'
            l1 -= 1
        elif current_score == left_score + gap_penalty:
            alignment1 += '-'
            alignment2 += sequence2[l2-1]
            l2 -= 1

    #finalizing alignments
    while l1 > 0:
        alignment1 += sequence1[l1-1]
        alignment2 += '-'
        l1 -= 1
    while l2 > 0:
        alignment1 += '-'
        alignment2 += sequence2[l2-1]
        l2 -= 1

    alignment1 = alignment1[::-1]
    alignment2 = alignment2[::-1]

    return(alignment1, alignment2)


def smith_waternam_with_alignment(alignment1, alignment2):
    pipes = ""
    score = 0
    #calculate score
    for i in range(0,len(alignment1)):
        if alignment1[i] == alignment2[i]:
            pipes += "|"
            score += 1
        elif alignment1[i]=="-" or alignment2[i]=="-":
            pipes += " "
            score -= 2
        else:
            pipes += " "
            score -= 1
    return score, pipes


def main(sequence1, sequence2):
    alignment1, alignment2 = smith_waternman(sequence1, sequence2)
    score, pipes = smith_waternam_with_alignment(alignment1, alignment2)
    print(alignment1 + "\n" + pipes + "\n" + alignment2)
    print("Score=", score)


if __name__ == "__main__":
    main("ggttgacta", "tgttacgg")
