#!/bin/bash
#Author: Evelyn Wong
#Date Created: 2023-10-28
#Description: Bash script to sort given sam file using samtools and outputting with desired name

#SBATCH --account=bgmp                    #REQUIRED: which account to use
#SBATCH --partition=bgmp               #REQUIRED: which partition to use
#SBATCH --mail-user=evew@uoregon.edu     #optional: if you'd like email
#SBATCH --mail-type=ALL                   #optional: must set email first, what type of email you want
#SBATCH --cpus-per-task=8                 #optional: number of cpus, default is 1
#SBATCH --mem=16GB                        #optional: amount of memory, default is 4GB

conda activate bgmp_star #environment with samtools

# dir1=/projects/bgmp/shared/deduper/Dataset1.sam
# out1=/projects/bgmp/evew/bioinfo/Bi624/Deduper-evelyn-n-wong/sam_output/Dataset1.out.sam
# dir2=/projects/bgmp/shared/deduper/Dataset2.sam
# out2=projects/bgmp/evew/bioinfo/Bi624/Deduper-evelyn-n-wong/sam_output/Dataset2.out.sam
# dir3=/projects/bgmp/shared/deduper/Dataset3.sam
# out3=projects/bgmp/evew/bioinfo/Bi624/Deduper-evelyn-n-wong/sam_output/Dataset3.out.sam
dir4=/projects/bgmp/shared/deduper/C1_SE_uniqAlign.sam
out4=/projects/bgmp/evew/bioinfo/Bi624/Deduper-evelyn-n-wong/sam_output/C1_SE_sorted.out.sam
order=/projects/bgmp/evew/bioinfo/Bi624/Deduper-evelyn-n-wong/chrom_order.txt
# test=/projects/bgmp/evew/bioinfo/Bi624/Deduper-evelyn-n-wong/test.sam
# out=/projects/bgmp/evew/bioinfo/Bi624/Deduper-evelyn-n-wong/test_sort.sam
/usr/bin/time -v samtools sort $dir4 -o $out4 -T $order
# /usr/bin/time -v samtools sort $test -o $out 