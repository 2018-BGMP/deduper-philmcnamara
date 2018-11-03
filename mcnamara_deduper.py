#! /usr/bin/env python3

import argparse

# ARGPARSE


def get_arguments():
    parser = argparse.ArgumentParser(description="Deduper")
    parser.add_argument("-u", "--umis",
                        help="text file input of UMI sequences, one per line",
                        required=True, type=str)
    parser.add_argument("-f", "--file", help="Input SAM file to remove PCR \
                        duplicates", required=True, type=str)
    parser.add_argument("-p", "--paired_end", help="Input SAM file contains \
                        paired-end reads", required=False, action="store_true")
    return parser.parse_args()


args = get_arguments()

if(args.paired_end):
    print("Not functional with paired end reads")
    print("Exiting...")
    exit()

# umidict will hold unique UMIs as keys, with sets holding unique positions
# mapping to those UMIs. After each chromosome is complete the dictionary is
# reset. Reads alignining to the minus strand will be stored as negative values
umidict = {}
with open(args.umis) as f:
    for line in f:
        line = line.strip()
        umidict[line] = set()
# Current chromosome
chromosome = ""

# Functions


def adjust_position(line):
    """Modifies the stored position at index -1 to accurately reflect \
    soft clipping. If insertions/delections/Ns are present, the position will \
    be adjusted to the side of the read with more matches/mismatches"""
    cigar = line[5]

    # Parse the cigar string into a list with elements as the components
    parsed = ""
    for char in range(len(cigar)):
        parsed += cigar[char]
        # Split on every letter except the last
        if((cigar[char] in "SMDIN") and (char != len(cigar) - 1)):
            parsed += "-"
    parsed = parsed.split('-')

    # Forward Strand
    if(line[-1] > 0):
        # Adjust real mapping position if insertions/deletions have occured
        for item in parsed:
            if("D" in item or "I" in item or "N" in item):
                # Position with highest value for M is assumed to be position
                max_matches = 0
                for x in parsed:
                    if("M" in x):
                        matches = int(x[:-1])
                        max_matches = max(matches, max_matches)
                # position of highest M
                pos = parsed.index(str(max_matches) + "M")
                # If not at the beginning of the cigar string we need to adjust
                if(pos != 0 and pos != 1):
                    # Add all left-most mapping to position
                    adjustment = 0
                    for i in range(pos):
                        c = parsed[i]
                        # Exclude soft-clipping adjustment
                        if("M" in c or "D" in c or "I" in c or "N" in c):
                            adjustment += int(parsed[i][:-1])
                    line[-1] += adjustment
        # Only adjust for soft-clipping at 5' end - S is not last char
        if(parsed[0][-1] == "S"):
            adjustment = int(parsed[0][:-1])
            line[-1] += adjustment
    # Reverse Strand
    else:
        for item in parsed:
            if("D" in item or "I" in item or "N" in item):
                # Position with highest value for M is assumed to be position
                max_matches = 0
                for x in parsed:
                    if("M" in x):
                        matches = int(x[:-1])
                        max_matches = max(matches, max_matches)
                # position of highest M, subtract from length to get other end
                pos = len(parsed) - parsed.index(str(max_matches) + "M")
                # If not at the end of the cigar string we need to adjust
                if(pos != 0 and pos != 1):
                    # Add all right-most mapping to position
                    adjustment = 0
                    for i in range(pos):
                        c = parsed[i]
                        # Exclude soft-clipping adjustment
                        if("M" in c or "D" in c or "I" in c or "N" in c):
                            adjustment += int(parsed[i][:-1])
                    line[-1] += adjustment
        for char_index in range(len(cigar)):
            # Only adjust for soft-clipping at 5' end - S is last char
            # Add, since reverse strand position is stored as a negative int
            if(parsed[-1][-1] == "S"):
                adjustment = int(parsed[-1][:-1])
                line[-1] += adjustment


with open(args.file, "r") as input, \
        open(args.file[:-4] + "_uniques.sam", "w") as uniq, \
        open(args.file[:-4] + "_duplicates.sam", "w") as dup, \
        open(args.file[:-4] + "_bad_umi.sam", "w") as badumi:
    for line in input:
        line = line.strip('\n').split('\t')
        if(line[0].startswith("@")):
            for field in line:
                uniq.write(field)
            uniq.write('\n')
        else:
            # Reset the dictionary if necessary
            if(line[2] != chromosome):
                for key in umidict:
                    umidict[key] = set()
            chromosome = line[2]
            # Add the position as an extra field at the end of the list
            # Number of fields is variable, so reference by line[-1]
            line.append(int(line[3]))
            # Extract and save the UMI
            umi = line[0][-8:]
            # Adjust for strand with bitwise flag, reverse strands are negative
            if(int(line[1]) & 16):
                line[-1] *= -1
            # Adjust the read's position
            adjust_position(line)
            # If UMI matches
            if(umi in umidict):
                # Non-duplicate
                if(line[-1] not in umidict[umi]):
                    # Add to set
                    umidict[umi].add(line[-1])
                    # Remove adjusted position before writing
                    line.pop()
                    line = '\t'.join(line)
                    uniq.write(line + '\n')
                # PCR Duplicate
                else:
                    line.pop()
                    line = '\t'.join(line)
                    dup.write(line + '\n')
            # Bad UMI
            else:
                line.pop()
                line = '\t'.join(line)
                badumi.write(line + '\n')
