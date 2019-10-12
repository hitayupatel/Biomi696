#! /usr/bin/env python

def burrow_wheeler_transform(sequence):
    bwt_matrix = []
    for i in range(len(sequence)):
        prefix, suffix = sequence[i:], sequence[:i]
        bwt_matrix.append(prefix + suffix)
    last_character = [b[-1] for b in sorted(bwt_matrix)]
    transform = ""
    for i in last_character:
        transform += i
    return transform


def main():
    sequence = "AATTGCGCGG$"
    print("Sequence = ", sequence)
    print("Burrow Wheeler Transform =",burrow_wheeler_transform(sequence)) 

if __name__ == '__main__':
    main()