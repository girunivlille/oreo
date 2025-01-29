import hashlib
import sys
import os
import argparse
from Bio import SeqIO

def main():
    parser = argparse.ArgumentParser(description='Integrity checker: checks if files hold the same biological content beside its order. Prints a MD5 for each file and exits 1 if they hold different content.')

    parser.add_argument('--readfiles', type=str, nargs=1, required=False,
                        help='Path of the fasta/fastq files to check, separated by commas (file1,file2..).')
    parser.add_argument('--seq_only', action='count', default=0, required=False,
                        help='Compares only sequences. Will ignore quality scores if one or several files are fastq.')
    parser.add_argument('--rc_sensitivity', action='count', default=0, required=False,
                        help='Sensitivity to reverse complement. If used, a sequence and its reverse complement will not be handled as the same sequence. By default they are handled as the same.')
    args = parser.parse_args()

    md5_size = 128
    output = ""
    all_md5 = []

    for file in args.readfiles[0].split(','):
        md5_sum = "0" * 32
        format = "fastq" if file.endswith("fastq") else "fasta"
        for record in SeqIO.parse(file, format):
            # Handle reverse complement
            if args.rc_sensitivity==0:
                seq = str(record.seq)
                rev_comp = str(record.seq.reverse_complement())
                sequence = min(seq, rev_comp)
                if format=="fastq":
                    if args.seq_only==0 and sequence==rev_comp :
                        qscore = (''.join(str(l) for l in record.letter_annotations["phred_quality"]))[::-1]
                    else:
                        qscore = ''.join(str(l) for l in record.letter_annotations["phred_quality"])
            else:
                sequence = str(record.seq)
                if format=="fastq":
                    qscore = ''.join(str(l) for l in record.letter_annotations["phred_quality"])

            # Hash header+sequence+qscore if format=fastq and not --seq_only else header+sequence 
            if format=="fastq" and args.seq_only==0:
                md5 = hashlib.md5((record.description+sequence+qscore).encode('utf-8')).hexdigest()
                md5_sum = hex((int(md5_sum, 16) + int(md5, 16)) % 2**md5_size)[2:].zfill(32)
            else:
                md5 = hashlib.md5((record.description+sequence).encode('utf-8')).hexdigest()
                md5_sum = hex((int(md5_sum, 16) + int(md5, 16)) % 2**md5_size)[2:].zfill(32)

        output += "\n"+file+" : "+str(md5_sum)
        all_md5.append(md5_sum)

    print(output)
    if len(set(all_md5)) == 1:
        print("\n-------------------------------------\nAll files are identical according to your criteria\n-------------------------------------")
        exit(0)
    else:
        print("\n-------------------------------------\nNot all the files are identical according to your criteria\n-------------------------------------")
        exit(1)

if __name__ == '__main__':
    main()