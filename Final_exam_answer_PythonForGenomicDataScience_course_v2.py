# -*- coding: utf-8 -*-
"""
Created on Wed Feb 12 14:32:06 2025

@author: linda.lam
Python code for Final exam for Python for Genomic Data Science in Coursera
"""

#import library
import re
#Store each asta file as dictonary; key is the identifier, value is the sequence
def set_sequence_asDict(input_file):
    lines = input_file.split("\n")   
    # Dictionary to store headers as keys and sequences as values
    sequences_dict = {}
    # Temporary variables to store the header and sequence
    current_header = None
    current_sequence = ""
    for line in lines:
        if line.startswith(">"):  # Check if the line is a header
            if current_header:  # Save the previous sequence to the dictionary
                sequences_dict[current_header] = current_sequence
            current_header = line  # Update the current header
            current_sequence = ""  # Reset the sequence for the new header
        else:
            current_sequence += line  # Append the sequence lines
    # After the loop, store the last sequence in the dictionary
    if current_header:
        sequences_dict[current_header] = current_sequence
    return sequences_dict

# gets the lonest and shortest sequence in the fasta file
def getlongest_shortest_seq(sequences_dict):
 
    min_length = float('inf')  # Set to positive infinity initially
    max_length = -float('inf')  # Set to negative infinity initially
    for key, seq in sequences_dict.items():
        seq_length = len(seq)
        # Check if the current sequence is shorter than the minimum length
        if seq_length < min_length:
            min_length = seq_length
        # Check if the current sequence is longer than the maximum length
        if seq_length > max_length:
            max_length = seq_length

    return min_length, max_length    

#returns the longest ORF in the fasta file
def get_longestORF_info_from_fasta(sequences_dict, reading_frame):
    ORF = "ATG"
    stop_codons = ["TAA", "TAG", "TGA"]
    length_ORF = -1
    longest_ORF_start_pos = -1
    for key, seq in sequences_dict.items():
        if reading_frame == 1:
                frame_seq = seq 
        elif reading_frame == 2:
                frame_seq = seq[1:]  
        elif reading_frame == 3:
                frame_seq = seq[2:]         
        # Initialize variables for scanning the frame
        in_ORF = False
        start_index = -1 
        # loop through the sequence in steps of 3 (triplet codons)
        for i in range(0, len(frame_seq) - 2, 3):  # Step by 3 to process codons
            codon = frame_seq[i:i+3]
            
            if not in_ORF:
                if codon == ORF:  # Start codon (ATG) found
                    in_ORF = True
                    start_index = i  # Remember the start index
            elif codon in stop_codons:  # Stop codon found
                length_of_seq = i + 3 - start_index  # Length of ORF
                if length_of_seq > length_ORF:
                    length_ORF = length_of_seq
                    longest_ORF_identifier = key
                    longest_ORF_start_pos = start_index + 3  
                in_ORF = False  # Reset for next ORF scan

    return length_ORF, longest_ORF_identifier, longest_ORF_start_pos

#returns the longest ORF for a perticular identifier
def get_longestORF_info_from_identifier(sequences_dict, sequence_identifier):
    ORF = "ATG"
    stop_codons = ["TAA", "TAG", "TGA"]
    length_ORF = -1
    longest_ORF_start_pos = -1
    longest_ORF_identifier = None
    
    # Add '>' to the sequence_identifier
    sequence_identifier = '>' + sequence_identifier.strip()
# Iterate through the sequences in the dictionary
    for header, seq in sequences_dict.items():     
        # Check if the sequence header starts with the identifier (case-insensitive)
        if re.match(f"^{re.escape(sequence_identifier)}", header.strip(), re.IGNORECASE):
            # Initialize variables for scanning the frame
            in_ORF = False
            start_index = -1

            # Scan through the sequence in steps of 3 (triplet codons)
            for i in range(0, len(seq) - 2, 3):  # Step by 3 to process codons
                codon = seq[i:i+3]

                if not in_ORF:
                    if codon == ORF:  # Start codon (ATG) found
                        in_ORF = True
                        start_index = i  # Remember the start index
                elif codon in stop_codons:  # Stop codon found
                    length_of_seq = i + 3 - start_index  # Length of ORF
                    if length_of_seq > length_ORF:
                        length_ORF = length_of_seq
                        longest_ORF_identifier = header  # Save header of the longest ORF
                        longest_ORF_start_pos = start_index + 3  
                    in_ORF = False  # Reset for next ORF scan

    return length_ORF, longest_ORF_identifier, longest_ORF_start_pos
    

def get_repeats(sequences_dict, length_repeat):
    # store sequence of length n as keys and the number of counts as value
    repeat_counts = {}
    for key, seq in sequences_dict.items():
        if len(seq) >= n:
            for i in range(len(seq) - n + 1):  # Extract substrings of length n
                repeat = seq[i:i+n]
                if repeat in repeat_counts:
                    repeat_counts[repeat] += 1  # Increment count if key exists
                else:
                    repeat_counts[repeat] = 1   # Initialize key if it doesn't exist

    # Find the most frequent repeat
    if repeat_counts:
        most_frequent_repeat = max(repeat_counts, key=repeat_counts.get)
        most_frequent_count = repeat_counts[most_frequent_repeat]
    else:
        most_frequent_repeat = None
        most_frequent_count = 0

    return repeat_counts, most_frequent_repeat, most_frequent_count



"""
run function and print out info
"""
"""
(1) How many records are in the file?
"""
input_file= open("dna2.fasta").read()
number_records = input_file.count(">")
print("The number of records are: ", number_records)   


""""
(2) What are the lengths of the sequences in the file? 
#What is the longest sequence and what is the shortest sequence? 
#Is there more than one longest or shortest sequence? 
#What are their identifiers? 
"""
seq_dict = set_sequence_asDict(input_file)
min_length, max_length = getlongest_shortest_seq(seq_dict)
print("Length of the smallest sequence:" , min_length)
print("Length of the longest sequence:" , max_length )

"""" 
(3 Given an input reading frame on the forward strand (1, 2, or 3) 
your program should be able to identify all ORFs present in each sequence of the FASTA file, 
and answer the following questions: 
what is the length of the longest ORF in the file? 
What is the identifier of the sequence containing the longest ORF? 
What is the starting position of the longest ORF in the sequence that contains it? 
"""
reading_frame = 2
longestORF_info = get_longestORF_info_from_fasta(seq_dict, reading_frame )       
print("\nLength of the longest ORF in the file:" , longestORF_info[0])
print("The identifier of the sequence containing the longest ORF:" + longestORF_info[1])   
print("Start postion of the longest ORF in the file:" ,longestORF_info[2])

"""
(4) For a given sequence identifier, what is the longest ORF contained in the sequence represented by that identifier?
""" 
sequence_identifier = 'gi|142022655|gb|EQ086233.1|16'
length_longestORF = get_longestORF_info_from_identifier(seq_dict, sequence_identifier)[0]
print("\nlongest ORF contained in the sequence represented by identifier", sequence_identifier, ":" ,length_longestORF) 


""""
(5) Given a length n, your program should be able to identify all repeats of length n in all sequences in the FASTA file. 
Your program should also determine how many times each repeat occurs in the file,
 and which is the most frequent repeat of a given length.
"""
n = 7
repeat_counts, most_frequent_repeat, most_frequent_count = get_repeats(seq_dict, n)
# print("Repeat Counts:", repeat_counts)
print("\nMost Frequent Repeat for length of ", n, ": ", most_frequent_repeat, "occurs", most_frequent_count, "times")

n = 12
repeat_counts, most_frequent_repeat, most_frequent_count = get_repeats(seq_dict, n)
# Count how many different sequences occur Max times
max_repeats = [seq for seq, count in repeat_counts.items() if count == most_frequent_count]
print("Number of copies with length", n, "that occurs",most_frequent_count, ":", len(max_repeats))

n = 6
repeat_counts, most_frequent_repeat, most_frequent_count = get_repeats(seq_dict, n)
# print("Repeat Counts:", repeat_counts)
print("Most Frequent Repeat for length of ", n, ": ", most_frequent_repeat, "occurs", most_frequent_count, "times")