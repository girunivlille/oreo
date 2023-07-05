import os
import argparse

def main():
    parser = argparse.ArgumentParser(description='Returns a new file containing the sorted reads.')

    parser.add_argument('--reads', type=str, nargs=1, required=True,
                        help='Path of the fasta file containing the reads that have to be sorted.')
    args = parser.parse_args()

    #copying the reads
    print("Copying the reads...")
    command = './copy_readsfile '+args.reads[0]
    os.system(command)

    #minimap2
    print("Mapping the reads...")
    command = 'minimap2/minimap2 -x ava-pb -t16 -k28 -w15 copie-'+args.reads[0]+' copie-'+args.reads[0]+' | gzip -1 > '+args.reads[0][:-3]+'.paf.gz'
    os.system(command)

    #miniasm
    print("Building contigs...")
    command = 'miniasm/miniasm -I1 -F1 -f '+args.reads[0]+' '+args.reads[0][:-3]+'.paf.gz > '+args.reads[0][:-3]+'.gfa'
    os.system(command)

    #reads sorting
    print("Sorting the reads...")
    command = './reads_sorting '+args.reads[0]+' '+args.reads[0][:-3]+'.gfa '+args.reads[0][:-3]+'_sorted.fa '+args.reads[0][:-3]+'.paf.gz'
    os.system(command)
    print("Sorted reads saved in "+args.reads[0][:-3]+"_sorted.fa.")

    #removing useless files (gfa, paf.gz, reads copy)
    print("Removing temporary files...")
    command = 'rm copie-'+args.reads[0]
    os.system(command)
    command = 'rm '+args.reads[0][:-3]+'.gfa'
    os.system(command)
    command = 'rm '+args.reads[0][:-3]+'.paf.gz'
    os.system(command)

if __name__ == '__main__':
    main()