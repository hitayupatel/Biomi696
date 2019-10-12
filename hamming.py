#! /usr/bin/env python

def hamming_distance(str1, str2):
    i, count = 0, 0
    while(i < len(str1)):
        if(str1[i] != str2[i]):
            count += 1
        i += 1
    return count


def main():
    distance = hamming_distance("TTATCGCGCTTTCTTCCAA","ATACCGCGCGTTCGACCAA") 
    print("Sequence 1 : TTATCGCGCTTTCTTCCAA")
    print("Sequence 2 : ATACCGCGCGTTCGACCAA")
    print("The hamming distance is = ", distance)


if __name__ == "__main__":
    main()