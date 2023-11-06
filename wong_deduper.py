#!/usr/bin/env python
# Author: Evelyn Wong
# Date Created: 2023-10-19
# Last Updated: 2023-10-30
''' Description: Given a sorted sam file, this program removes all PCR duplicates and retains a single copy of each read from SAM file.
                Input: sorted sam file <sorted_sam_file.sam> and text file with known list of UMIs <known_umi_file.txt>
                Output: user-defined output sam file, sam stats files which contain counts of chromosomes, unique reads, duplicate reads, and
                        unknown UMIs'''

# importing modules
import argparse
import sys
#import bioinfo if I want to check quality
import re

def get_args():
    parser = argparse.ArgumentParser(description="This function takes in user arguments for demultiplexing. Arguments required: -f, -o, -u")
    parser.add_argument("-f", "--file", help = "designate absolute file path to sorted sam file", required = True)
    parser.add_argument("-o", "--outfile", help = "designates absolute file path to sorted sam output file", required = True)
    parser.add_argument("-u", "--umi", help = "designates file containing list of known UMIs", required = True)
    args = parser.parse_args()
    return args

# Setting up variables 

args = get_args()

file = args.file
out = args.outfile
umi_file = args.umi

line_counter = 0 # checking current line 
uniq_counter = 0 # count number of unique reads
dup_counter = 0 # count number of duplicates
unknown_umi_counter = 0 # count number of unknown UMIs

umi_set = set() # set to hold known umi's from file and used for comparison 
sam_dict: dict = {} # dictionary to output to sam file; key is a tuple (umi, chromosome, strand, position), value is the read 
duplicate_chrom: dict = {} # dictionary to keep track of every new chromosome; reset each time there's new
unique_chrom: dict = {} # dictionary to keep track of every chromosome counts

def known_umi(umi_seq: str, umi_set):
    '''This function returns TRUE for current seq with known UMIs that are a part of the known 96 UMIs given.
    If not part of known list then, returns FALSE'''
    return umi_seq in umi_set

def soft_clipping(cigar_string: str, position: int, strand: str):
    '''This function accounts for any possible soft-clipping and returns the corrected position to be used for comparison in the dictionary.
    Accounts for strandedness: if forward strand, subtract soft-clipped end form 5' end from position given.
    If reverse strand, add soft-clipped end to 5' end position.'''
    position: int
    new_pos: int
    # If there is soft-clipping
    if 'S' in cigar_string:
        # Get line
        cigar_temp = re.findall('\d+[A-Z]{1}', cigar_string)
        # Check orientation
        # if forward stand, check clipped amount and subtract from position
        if strand == '+': # if forward strand
            if 'S' in cigar_temp[0]: # if S in first position (leftmost)
                clipped = int(cigar_temp[0].split('S')[0]) # amount clipped
                
                print(f'position: {position}, clipped: {clipped}')

                new_pos = int(position) - int(clipped) # new position

                print(f'new_pos: {new_pos}')
            else:
                new_pos = int(position) #if not leftmost, return position as-is
        
        # reverse strand, add alignment matches (M), deletions (D), skipped regions (N), and soft-clipping (S) on right-most end to be removed
        else:
            new_pos = int(position) # setting position as current position
            adjust = re.findall('([0-9]+)[MDN]', cigar_string) # look for instances of MDN which is returned as a list
            for i in adjust:
                new_pos += int(i) #adding up the positions 
            soft_rh = re.search('\d+S$', cigar_string) # look for soft-clipping on right-hand side
            if soft_rh: # if True, extract matched text, remove the 'S' at the end, and convert to integer
                new_pos += int(soft_rh.group(0)[:-1]) - 1 # subtract 1 to adjust for position
            else:
                new_pos = new_pos - 1 # no soft-clipping on right-hand side
   
   # No soft-clipping, just check strand orientation and adjust for position 
    else:
        if strand == "+":
            new_pos = int(position)
        else: # reverse strand
            new_pos = int(position) # setting position as current position
            adjust = re.findall('([0-9]+)[MDN]', cigar_string) # look for instances of MDN
            for i in adjust:
                new_pos += int(i) #adding up the positions 
            soft_rh = re.search('\d+S$', cigar_string) # look for soft-clipping on right-hand side
            if soft_rh: # if True, extract matched text, remove the 'S' at the end, and convert to integer
                new_pos += int(soft_rh.group(0)[:-1]) - 1 # subtract 1 to adjust for position
            else:
                new_pos = new_pos - 1 # new position - 1 to adjust for position

    return new_pos

def stranded(flag):
    '''This function checks whether the strand's bitflag is forward or reverse'''
    if ((int(flag) & 16) == 16):
        return '-' #reverse
    else:
        return '+' #forward


# setting up UMI set 
with open(umi_file, "r") as fh:
    for line in fh:
        #remove whitespace
        umi = line.strip()
        #add UMI to set 
        umi_set.add(umi)

#checking umi_set
print(umi_set)

# main 
with open(file, "rt") as in_file, open(f'{out}.sam', "w") as out_file:
    while True:
        line = in_file.readline()
        line_counter +=1

        if not line: 
            break #EOF 

        # Checking progress [C1_SE... file has 18186474]
        #if line_counter % 10 == 0: # for test file
        if line_counter % 10000 == 0:
            print(f'Lines processed so far: {line_counter}')

        if line[0] == "@": # write header to output file 
            out_file.write(line)
            #doppels.write(line + '\n')
        else:
            fields = line.split("\t")
            qname = fields[0] # read name 'e.g. NS500451:154:HWKTMBGXX:1:11101:15364:1139:GAACAGGT"
            flag = fields[1] # bitwise flag
            chrom = fields[2] # chrom/reference name
            pos = fields[3] # leftmost mapping position
            cigar = fields[5] # cigar string 
            umi = qname.split(":")[7] # umi e.g."GAACAGGT"

            #Checking if UMI is in known UMIs
            if not known_umi(umi, umi_set):
                unknown_umi_counter += 1 
                continue
            else:
                # Checking current chromosome in duplicate dictionary
                if chrom in duplicate_chrom.keys():
                    
                    # Checking strandedness:
                    strand = stranded(flag) # either '-' or '+'
                    # Adjusting position
                    new_pos = soft_clipping(cigar, pos, strand)

                    # if current chromosome with specific umi, strand, and corrected position not in dictionary, add to dictionary, and increment unique counter 
                    if (umi, chrom, strand, new_pos) not in sam_dict.keys():
                        sam_dict[umi, chrom, strand, new_pos] = line
                        uniq_counter += 1
                        unique_chrom[chrom] = unique_chrom.get(chrom, 1) + 1
                    # if current chromosome with umi/strand/position is in dictioary, increment duplicate counter for user statistics 
                    else:
                        dup_counter += 1
                        duplicate_chrom[chrom] += 1

                else: # different chromosome + increment
                    duplicate_chrom[chrom] = 1
                    # sam_dict.clear()
                    

    for key, value in sam_dict.items(): # outputting line (dictionary value) to output sam file (unique read)
        out_file.write(f'{value}')

# outputting a user stats files
with open("deduper_statistics.txt", "wt") as stats:
    stats.write(f'Statistics Summary\n')
    stats.write(f'----------------------------------\n')
    stats.write(f'Number of unknown UMIs: {unknown_umi_counter}\n')
    stats.write(f'Number of unique reads: {uniq_counter}\n')
    stats.write(f'Number of duplicates: {dup_counter}\n')
    for chrom, counter_val in duplicate_chrom.items():
        stats.write(f'Duplicate Chrom: {chrom} Number Duplicates: {counter_val}\n')
    stats.write(f'----------------------------------\n')
    for chrom, counter_val in unique_chrom.items():
        stats.write(f'Unique Chrom: {chrom} Number Unique Reads: {counter_val}\n')