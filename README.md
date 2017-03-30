# BioProject2
Aligning sequences with Needleman-Wunsch or Smith-Waterman(https://github.com/Droom1995/BioProject2)

# Algorithm description
For both algorithms, standard DP approach with insert/delete/match functions was implemented.  
For further alignment, I remember every step for DP matrix in a different path matrix with next values:   
0 - restart path(for local alignment);  
1 - up(insertion/deletion);  
2 - left(deletion/insertion);  
3 - diagonal (match);  

Algorithms with supporting functions are listed in alignment.py source file. Alignment algorithms accept next parameters: seq1 - first sequence, seq2 - second sequence, letters - letter order for alignment matrix, blosum62 - alignment matrix, gap_p - gap penalty.

# Execution results
./main.py -g -s1 A0PQ23.fasta -s2 Q9CD83.fasta -e blosum62.csv -p 4  
M--T---NR--T--LSREEIRKLDRDLRILVATNGTLTRVLNVVANEEIVVDIINQQLLDVAPKIPELENLKIGRILQRDILLKGQKSGILFVAAESLIVIDLLPTAITTYLTKTHHPIGEIMAASRIETYKEDAQVWIGDLPCWL-ADYGYWDLPKRAVGRRYRIIAGGQPVIITTEYFLRSVFQDTPREELDRCQYSNDIDTRSGDRFVLHGRVFKNL  
MLAVLPEKREMTECHLSDEEIRKLNRDLRILIATNGTLTRILNVLANDEIVVEIVKQQIQDAAPEMDGCDHSSIGRVLRRDIVLKGRRSGIPFVAAESFIAIDLLPPEIVASLLETHRPIGEVMAASCIETFKEEAKVWAGESPAWLELDRRR-NLPPKVVGRQYRVIAEGRPVIIITEYFLRSVFEDNSREEPIRHQRS--VGT-SA-R---SGRSICT-  
545.0  
  
./main.py -l -s1 EU078679.fasta -s2 CH954156.fasta -e blastmatrix.csv -p 4  
TTGACAGTACATAG  
TTGA-AGTTTGTAG  
34.0  
