# -*- coding: utf-8 -*-
"""
Written for Software Design 2016 at Olin College of Engineering.  Gene finder takes a sequence of DNA and returns a list of possible amino acids.

@author:  March Saper


"""

import random
from amino_acids import aa, codons, aa_table   # you may find these useful
from load import load_seq


def shuffle_string(s):
    """Shuffles the characters in the input string
        NOTE: this is a helper function, you do not
        have to modify this in any way """
    return ''.join(random.sample(s, len(s)))


def get_complement(nucleotide):
    """ Returns the complementary nucleotide

        nucleotide: a nucleotide (A, C, G, or T) represented as a string
        returns: the complementary nucleotide
    >>> get_complement('A')
    'T'
    >>> get_complement('C')
    'G'
    >>> get_complement('T')
    'A'
    >>> get_complement('G')
    'C'
    """
    two_nucleotide = nucleotide.replace('A', 'Placeholder_one')             # The placeholder waits to change "A" to "T" until after the "T"'s have been counted
    three_nucleotide = two_nucleotide.replace('T', 'A')
    four_nucleotide = three_nucleotide.replace('Placeholder_one', 'T')
    five_nucleotide = four_nucleotide.replace('G', 'Placeholder_two')       # This placeholder waits to change "G" to "C" until after the "C"'s have been counted.
    six_nucleotide = five_nucleotide.replace('C', 'G')
    final_nucleotide = six_nucleotide.replace('Placeholder_two', 'C')  
    return final_nucleotide


def get_reverse_complement(dna):
    """ Computes the reverse complementary sequence of DNA for the specfied DNA
        sequence

        dna: a DNA sequence represented as a string
        returns: the reverse complementary DNA sequence represented as a string
    >>> get_reverse_complement("ATGCCCGCTTT")
    'AAAGCGGGCAT'
    >>> get_reverse_complement("CCGCGTTCA")
    'TGAACGCGG'
    >>> get_reverse_complement("CCTTAAGG")
    'CCTTAAGG'
    """
    complement = get_complement(dna)
    reverse_complement = complement[::-1]  
    return reverse_complement      


def rest_of_ORF(dna):
    """ Takes a DNA sequence that is assumed to begin with a start
        codon and returns the sequence up to but not including the
        first in frame stop codon.  If there is no in frame stop codon,
        returns the whole string.

        dna: a DNA sequence
        returns: the open reading frame represented as a string
    >>> rest_of_ORF("ATGTGAA")
    'ATG'
    >>> rest_of_ORF("ATGAGATAGG")
    'ATGAGA'
    >>> rest_of_ORF("ATATATATAT")
    'ATATATATAT'
    >>> rest_of_ORF("ATGATATG")
    'ATGATATG'
    """
    for i in range(0, len(dna)/3):
        n = i*3                                                                 # You want to go through the function in units of three.  This does that.
        if dna[n:n+3] == 'TAG' or dna[n:n+3] == 'TAA' or dna[n:n+3] == 'TGA':
            return dna[:n]                                                      # Return the dna sequence from the start to the last nucleotide before the stop codon.
        
    return dna    # If dna sequence makes it through the for loop without finding any stop codons, then the whole dna sequence will print.
        

def find_all_ORFs_oneframe(dna):
    """ Finds all non-nested open reading frames in the given DNA
        sequence and returns them as a list.  This function should
        only find ORFs that are in the default frame of the sequence
        (i.e. they start on indices that are multiples of 3).
        By non-nested we mean that if an ORF occurs entirely within
        another ORF, it should not be included in the returned list of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_oneframe("ATGCATGAATGTAGATAGATGTGCCC")
    ['ATGCATGAATGTAGA', 'ATGTGCCC']
    >>> find_all_ORFs_oneframe("ATGCATTAGTTGATGAATTGA")
    ['ATGCAT', 'ATGAAT']
    """
    all_ORFs = []    
    i = 0   
    while i < len(dna): 
        if dna[i:i+3] == 'ATG':    # checks to see if index lands on a start codon.  If it does, trims the string at the corresponding stop codon
            manipulated_string = rest_of_ORF(dna[i:len(dna)])
            i = i + len(manipulated_string)
            all_ORFs.append(manipulated_string)    
        i = i + 3   
    return all_ORFs


def find_all_ORFs(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence in
        all 3 possible frames and returns them as a list.  By non-nested we
        mean that if an ORF occurs entirely within another ORF and they are
        both in the same frame, it should not be included in the returned list
        of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs

        >>> find_all_ORFs("ATGCATGAATGTAG")
        ['ATGCATGAATGTAG', 'ATGAATGTAG', 'ATG']
        >>> find_all_ORFs("GGATGATATAATGTAA")
        ['ATG', 'ATGATA']
        >>> find_all_ORFs("TGACATGTAG")
        ['ATG']
        >>> find_all_ORFs("CTACATGTCA")
        ['ATGTCA']
    """   
    all_ORFs_all_frames = []                                    
    frame_one = find_all_ORFs_oneframe(dna[0:len(dna)])         # a list of ORFs starting in the first reading frame
    frame_two = find_all_ORFs_oneframe(dna[1:len(dna)])         # a list of all ORFs starting in the second reading frame
    frame_three = find_all_ORFs_oneframe(dna[2:len(dna)])       # a list of all ORFs starting in the second reading frame
    all_ORFs_all_frames = frame_one + frame_two + frame_three   # a list of all the above 
    return all_ORFs_all_frames


def find_all_ORFs_both_strands(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence on both
        strands.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
        >>> find_all_ORFs_both_strands("ATGCGAATGTAGCATCAAA")
        ['ATGCGAATG', 'ATGCTACATTCGCAT']
        >>> find_all_ORFs_both_strands("TGACATGTAG")
        ['ATG', 'ATGTCA']
        >>> find_all_ORFs_both_strands("AAAAAAA")
        []
    """
    everything_dna_related = []    
    complement_dna = get_reverse_complement(dna)   
    ORF_on_origional = find_all_ORFs(dna)                      # finds a list of all the ORFs in the dna
    ORF_on_complement = find_all_ORFs(complement_dna)          # finds a list of all the ORFs in the dna complement
    everything_dna_related = ORF_on_origional + ORF_on_complement
    return everything_dna_related


def longest_ORF(dna):
    """ Finds the longest ORF on both strands of the specified DNA and returns it
        as a string
    >>> longest_ORF("ATGCGAATGTAGCATCAAA")
    'ATGCTACATTCGCAT'
    >>> longest_ORF("TGACATGTAG")
    'ATGTCA'
    >>> longest_ORF("AAAAAA")
    ''
    """
    dna_list = find_all_ORFs_both_strands(dna)                 # list of all of the possible ORFs in the input dna
    if len(dna_list) != 0:                                     # if ORFs are found continue on
        longest_ORF_str = max(dna_list, key=len)               # finds the longest ORF in the list 
        return longest_ORF_str    
    else:                                                      # if no ORFs are found return nothing
        return ''


def longest_ORF_noncoding(dna, num_trials):
    """ Computes the maximum length of the longest ORF over num_trials shuffles
        of the specfied DNA sequence

        dna: a DNA sequence
        num_trials: the number of random shuffles
        returns: the maximum length longest ORF 
        >>> longest_ORF_noncoding("AAAAAA", 10)
        0
    """   
    longest_ORF_len = 0    
    for i in range(num_trials):
        shuffled_dna = shuffle_string(dna)                     # produces a new list of dna for each trial
        shuffled_longest_ORF = str(longest_ORF(shuffled_dna))  # finds longest ORF in the newly shuffled dna
        if len(shuffled_longest_ORF) > longest_ORF_len:        # if the length of the longest ORF is greater than the current length of the longest ORF updates the number
            longest_ORF_len = len(shuffled_longest_ORF)
    return longest_ORF_len    


def coding_strand_to_AA(dna):
    """ Computes the Protein encoded by a sequence of DNA.  This function
        does not check for start and stop codons (it assumes that the input
        DNA sequence represents an protein coding region).

        dna: a DNA sequence represented as a string
        returns: a string containing the sequence of amino acids encoded by the
                 the input DNA fragment

        >>> coding_strand_to_AA("ATGCGA")
        'MR'
        >>> coding_strand_to_AA("ATGCCCGCTTT")
        'MPA'
        >>> coding_strand_to_AA("AAAAAAAA")
        'KK'
        >>> coding_strand_to_AA("")
        ''
        >>> coding_strand_to_AA("TAT")
        'Y'
    """
    amino_acid_list = ''

    for i in range(0, len(dna)/3):
        n = i*3    
        current_codon = dna[n:n+3]    
        amino_acid = aa_table[current_codon]                 # The codon is translated into the appropriate amino acid
        amino_acid_list = amino_acid_list + amino_acid    
    return amino_acid_list



def gene_finder(dna):
    """ Returns the amino acid sequences that are likely coded by the specified dna

        dna: a DNA sequence
        returns: a list of all amino acid sequences coded by the sequence dna.
        >>> gene_finder("AAAAAAA")
        []
    """
    threshold = longest_ORF_noncoding(dna, 100)         # computes the threshold length that a sequence of open reading frame must be in order to be considered as a gene
    all_the_dna = find_all_ORFs_both_strands(dna)    
    long_enough_dna = []    
    for item in range(len(all_the_dna)):
        if len(all_the_dna[item]) > threshold:
            long_enough_dna.append(all_the_dna[item])    # if the ORF in the list of ORFs is longer than the threshold it is added to a possible list of long enough ORF

    final_amino_acid_list = []    
    for item in range(len(long_enough_dna)):
        current_dna_strand = long_enough_dna[item]    
        amino_acid = coding_strand_to_AA(current_dna_strand)    
        final_amino_acid_list.append(amino_acid)    

    return final_amino_acid_list


if __name__ == "__main__":
    # import doctest
    # doctest.testmod()

    from load import load_seq
    dna = load_seq("./data/X73525.fa")
    print gene_finder(dna)