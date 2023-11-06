#!/bin/bash
#Author: Evelyn Wong
#Date Created: 2023-11-01
#Description: Bash script to sort given sam file using samtools and outputting with desired name

#SBATCH --account=bgmp                    #REQUIRED: which account to use
#SBATCH --partition=bgmp               #REQUIRED: which partition to use
#SBATCH --mail-user=evew@uoregon.edu     #optional: if you'd like email
#SBATCH --mail-type=ALL                   #optional: must set email first, what type of email you want
#SBATCH --cpus-per-task=8                 #optional: number of cpus, default is 1
#SBATCH --mem=16GB                        #optional: amount of memory, default is 4GB

script=/projects/bgmp/evew/bioinfo/Bi624/Deduper-evelyn-n-wong/part_3/wong_deduper.py
# test_file=/projects/bgmp/evew/bioinfo/Bi624/Deduper-evelyn-n-wong/test_sort.sam
file=/projects/bgmp/evew/bioinfo/Bi624/Deduper-evelyn-n-wong/sam_output/C1_SE_sorted.out.sam
# test_out=/projects/bgmp/evew/bioinfo/Bi624/Deduper-evelyn-n-wong/sam_output/test_out.sam
out=/projects/bgmp/evew/bioinfo/Bi624/Deduper-evelyn-n-wong/part_3/C1_SE_deduped.sam
umi=/projects/bgmp/evew/bioinfo/Bi624/Deduper-evelyn-n-wong/STL96.txt

#/usr/bin/time -v python $script -f $test_file -o $test_out -u $umi
/usr/bin/time -v python $script -f $file -o $out -u $umi