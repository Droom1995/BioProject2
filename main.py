from Bio import SeqIO
import pandas as pd
import argparse
import alignment as al





argparser = argparse.ArgumentParser(description='Sequence alignment help')

argparser.add_argument('-g', '--global',
                       action="store_true", dest="algo_type",
                       help="global/local", default="False")
argparser.add_argument('-l', '--local',
                       action="store_false", dest="algo_type",
                       help="global/local", default="True")
argparser.add_argument('-s1', '--seq1',
                       action="store", dest="seq1",
                       help="Seq1 file name")
argparser.add_argument('-s2', '--seq2',
                       action="store", dest="seq2",
                       help="Seq1 file name")
argparser.add_argument('-e', '--score_matrix',
                       action="store", dest="blosum",
                       help="Score matrix")
argparser.add_argument('-p', '--p', type=int,
                       action="store", dest="gap_p",
                       help="Gap penalization")
argparser.add_argument('-pe', '--pe',
                       action="store", type=int, dest="gap_pe",
                       help="Gap penalization extended", default=4)
args = argparser.parse_args()

records = list(SeqIO.parse(args.seq1, "fasta"))
inp_seq1 = records[0].seq
records = list(SeqIO.parse(args.seq2, "fasta"))
inp_seq2 = records[0].seq

df1 = pd.read_csv(args.blosum)
gap_penalty = -args.gap_p
gap_penalty_ext = -args.gap_pe
blosum62_init = df1.values
blosum62 = []
letters = []
for arr in blosum62_init:
    arr = arr.tolist()
    letters.append(arr[0])
    arr = arr[1:]
    blosum62.append(arr)
ans = None
if args.algo_type == True:
    ans = al.needlman_wunsch(inp_seq1, inp_seq2, letters, blosum62, gap_penalty)
else:
    ans = al.smith_waterman(inp_seq1, inp_seq2, letters, blosum62, gap_penalty)
print(ans[0])
print(ans[1])
print(ans[2])
