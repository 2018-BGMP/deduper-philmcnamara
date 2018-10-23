## The Problem:
Given a SAM output file, sorted by leftmost mapping position using SAMtools

We will assume that the library prep was strand-specific.  

Remove all PCR duplicates ("de-dupe") and write the non-duplicate alignments to a new output file
PCR duplicates are sequences where ALL of the below are true: 

        1. Align to the same leftmost mapping position after adjusting for soft-clipping
        2. Align to the same strand ("+" or " -")
        3. Align to the same chromosome
        4. Are barcoded with the same Unique Molecular Identifier (UMI), applied prior to PCR

## Strategy:
SAMtools will sort the reads by chromosome and leftmost mapping position. This allows us to move through the file iteratively, comparing the data on each line. We will read in the data for each sequence, and if all 4 criteria above match the previous read, it will be discarded as a duplicate. The next consideration is how many lines to hold in memory at a time. Because we need to adjust for soft clipping, we cannot be sure a read is not a PCR duplicate until we have read in *n* subsequent reads, where *n* is the maximum possible amount of soft clipping (and alignment position adjustment) that could occur. A sliding window of 30 reads will be our starting point, but this can be increased very easily  (with some performance hits, since we need to iterate through the window for each line of the file).

## Implementation:
The sliding window will be implemented as a 30-node linked list, with each node holding a 24-element array-list. Each node stores a sequence line currently being compared, and the first 21 array-list elements will store the tab-delimited values of the SAM file. We will never need to do index-based lookup, so the convenience of "sliding" the window by removing the head node and adding a new node at the tail should save a good amount of time compared to updating the index of every element.

Certain values used for comparison require function calls and, because we are using a sliding window, 29 of the 30 values will stay the same for every loop. To minimize these calls, we will store the calculated values in extra fields at the end of the array-list. Indexes 0-20 hold the fields from the sequence that will be written to the output file, 21 will hold the strand (+ or -), 22 will hold the adjusted position after considering softclipping, and 23 will hold the UMI. 

A dictionary of UMIs will be created from an input file. We will only consider reads PCR duplicates if their UMIs match AND those UMIs are "real" (they appear in the dictionary). This means alignments with unknown UMIs still appear in the output file. They can be filtered later, but we keep them so the only modifications made by deduper will be the removal of PCR duplicates.

#### Argparse:

-f --file, input file, required

-p --paired, optional (my code errors if you pass this parameter)

-u --umi, umi text file, optional (my code doesn't deal with randomers so it errors out if not provided)

-h --help, print help text

#### Functions:

```python
def get_strand(bit_flag)
'''Input: the bitwise flag of the sequence. Convert to binary and check the 16 bit. If 1, store "-" in list index 21, otherwise store "+"'''

def adjust_softclip(sequence_line)
'''Look at CIGAR string, if left-sided soft-clipping is present, subtract the soft-clipping value from the mapping position. Store the result in array-list index 22. If there is no left-sided softclipping, store the standard mapping position instead'''

def extract_UMI(sequence_line)
'''Slice the last 8 characters off the header field, if the UMI appears in the dictionary store it in array-list index 23.'''

```



### Main:

Open input file for reading, *_deduped output file for writing, _duplicates file for writing

Initialize linked list, read in the first 30 lines of the file that don't start with @, splitting by tabs

Call get_strand(), adjust_softclip(), extract UMI() on each of the 30 initial sequences

**To assess if the read at the head node is a PCR duplicate of any of the following 29:**

target = head node

For node 2 to 30 (Iterate through the linked list):

​	If node's chromosome == target chromosome

​		If node's softclip-adjusted position == target's softclip-adjusted position

​			If node's strand == target's strand

​				If node's UMI == target's UMI

​					**PCR DUPLICATE**

​					**Delete head node, join array-list indexes 0-20 and write to "duplicates"**

​	else

​		**non-duplicate**

​		**Join array-list indexes 0-20 and write to "deduped"**

​		**Delete head node**

The linked list now has 29 elements, with the first sequence removed and the second sequence at **head** 

Call get_strand(), adjust_softclip(), extract UMI() on the next line of the file and store it in a new node at the tail (now position 30)

Continue until we reach the end of the file. 

​					

### Function Unit Tests:

``` python
test_alignment_1 = "NS500451:154:HWKTMBGXX:1:11101:24936:1293:TGAGTGAG	0	2	52159545	36	25M1086N46M	*	0	0	GTCCGATTGCTTCTTTATACTTGACCGAGCTAATGTGGTCCTGCGTTTGCTTCACTCTGAGCATCTCAGGC	6AEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEA	MD:Z:71	NH:i:1	HI:i:1	NM:i:0	SM:i:36	XQ:i:40	X2:i:0	XO:Z:UU	XS:A:-	XG:Z:A".strip('\n').split('\t')

test_alignment_2 = "NS500451:154:HWKTMBGXX:1:11202:20775:15962:CTAGGAAG	31	2	76847388	36	10S61M	*	0	0	ACTTCAGGCACTTTTGCAGGAGGAGCCTCTACTTTCTTAGGAGTAGGCGCTGGTACCTTCTTCTCAGGCAC	EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE	MD:Z:61	NH:i:1	HI:i:1	NM:i:0	SM:i:36	XQ:i:40	X2:i:0	XO:Z:UU	XS:A:-	XA:Z:15513,17091,16907,16715,5761,19621".strip('\n').split('\t')

# Add extra data to the list
get_strand(test_alignment_1)
get_strand(test_alignment_2)
adjust_softclip(test_alignment_1)
adjust_softclip(test_alignment_2)
extract_UMI(test_alignment_1)
extract_UMI(test_alignment_2)

#Both reads are on + strand
assert(test_alignment_1[21] == "+")
assert(test_alignment_2[21] == "-")

#No soft clipping
assert(test_alignment_1[22] == "52159545")
#Adjusted -10 for soft clipping
assert(test_alignment_2[22] == "76847378")

#UMIs
assert(test_alignment_1[23] == "TGAGTGAG")
assert(test_alignment_2[23] == "CTAGGAAG")
```

## File Unit Tests

See "test_input_1"

​			



