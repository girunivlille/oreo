#!/usr/bin/env python

import os
import sys
import argparse

str_version = 1.0

def main():
    parser = argparse.ArgumentParser(description='Returns a new file containing the sorted reads.')

    parser.add_argument('--reads', type=str, nargs=1, required=False,
                        help='Path of the fasta file containing the reads that have to be sorted.')
    parser.add_argument('-v', '--version', action='count', default=0, required=False,
            help='Prints the used versions of minimap2, miniasm and sort_the_reads.')
    parser.add_argument('--ctg_sort', type=int, nargs=1, required=False, default=[1],
                        help='Algorithm used to sort the contigs : 0:random order, 1:depth-first search, 2:breadth-first search. Default=1')
    parser.add_argument('--opt_minimap', type=str, nargs=1, required=False, default=["-x ava-pb -t16 -k28 -w15"],
                        help='String containing all options to run minimap2. Default=-x ava-pb -t16 -k28 -w15')
    parser.add_argument('--opt_miniasm', type=str, nargs=1, required=False, default=["-I1 -F1"],
                        help='String containing all options to run miniasm. Default=-I1 -F1')
    args = parser.parse_args()

    script_file = os.path.abspath(sys.argv[0])
    script_dir = os.path.dirname(script_file)
    minimap2_file = os.path.join(script_dir,'minimap2','minimap2')
    miniasm_file = os.path.join(script_dir,'miniasm','miniasm')
    copy_readsfile_file = os.path.join(script_dir,'copy_readsfile')
    reads_sorting_file = os.path.join(script_dir,'reads_sorting')

    if len(sys.argv)==1:
        parser.print_help()
    else:
        if args.version>0:
            #afficher les versions de minimap2, miniasm, sort_the_reads
            print("Minimap2 version :")
            os.system(minimap2_file + " --version")
            print("Miniasm version :")
            os.system(miniasm_file + " -V")
            print("Sort_the_reads version:")
            print(str_version)
        else:
            if(args.reads!=None):
                if args.ctg_sort[0] in [0,1,2]:
                    #retrieve the folder from where the file is executed
                    reads_file = os.path.abspath(args.reads[0])
                    readscopy_file = os.path.join(os.path.dirname(reads_file),"copie-"+os.path.basename(reads_file))
                    paf_file = os.path.splitext(reads_file)[0] + '.paf.gz'
                    gfa_file = os.path.splitext(reads_file)[0] + '.gfa'
                    reads_sorted_file = os.path.splitext(reads_file)[0] + '_sorted' + os.path.splitext(reads_file)[1]

                    #copying the reads
                    print("Copying the reads...")
                    command = copy_readsfile_file + ' ' + reads_file + ' ' + readscopy_file

                    os.system(command)

                    #minimap2
                    print("Mapping the reads...")
                    command = minimap2_file+' '+args.opt_minimap[0]+' '+readscopy_file+' '+readscopy_file+' | gzip -1 > '+paf_file
                    os.system(command)

                    #miniasm
                    print("Building contigs...")
                    #command = 'miniasm/miniasm '+args.opt_miniasm[0]+' -f '+args.reads[0]+' '+args.reads[0][:-3]+'.paf.gz > '+args.reads[0][:-3]+'.gfa'
                    command = miniasm_file+' '+args.opt_miniasm[0]+' -f '+readscopy_file+' '+paf_file+' > '+gfa_file
                    os.system(command)

                    #reads sorting
                    print("Sorting the reads...")
                    command = reads_sorting_file + ' ' + reads_file + ' ' + gfa_file + ' ' + reads_sorted_file + ' ' + paf_file + ' ' + str(args.ctg_sort[0])
                    os.system(command)
                    print("Sorted reads saved in "+reads_sorted_file)

                    #removing useless files (gfa, paf.gz, reads copy)
                    print("Removing temporary files...")
                    command = 'rm '+readscopy_file
                    os.system(command)
                    command = 'rm '+gfa_file
                    os.system(command)
                    command = 'rm '+paf_file
                    os.system(command)
                else:
                    print("Wrong number for --ctg_sort. Please write 0, 1 or 2.")
            else:
                parser.print_help()
                print("Please fill the --reads argument.")

if __name__ == '__main__':
    main()
