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

    if len(sys.argv)==1:
        os.system("python3 sort_the_reads.py -h")
    else:
        if args.version>0:
            #afficher les versions de minimap2, miniasm, sort_the_reads
            print("Minimap2 version :")
            os.system("minimap2/minimap2 --version")
            print("Miniasm version :")
            os.system("miniasm/miniasm -V")
            print("Sort_the_reads version:")
            print(str_version)
        else:
            if(args.reads!=None):
                if args.ctg_sort[0]==0 or args.ctg_sort[0]==1 or args.ctg_sort[0]==2:
                    #retrieve the folder from where the file is executed
                    path = args.reads[0].split("/")
                    position = sys.argv[0].split("/")
                    if len(position)==1:
                        folder ="./"
                    else:
                        folder=position[0]
                        for i in range(1, len(position)-1):
                            folder=os.path.join(folder, position[i])
                    copyfile ="copie-"+path[-1]

                    #copying the reads
                    print("Copying the reads...")
                    command = os.path.join(folder,'copy_readsfile')+' '+args.reads[0]
                    os.system(command)

                    #minimap2
                    print("Mapping the reads...")
                    command = os.path.join(folder, 'minimap2', 'minimap2')+' '+args.opt_minimap[0]+' '+copyfile+' '+copyfile+' | gzip -1 > '+args.reads[0][:-3]+'.paf.gz'
                    os.system(command)

                    #miniasm
                    print("Building contigs...")
                    #command = 'miniasm/miniasm '+args.opt_miniasm[0]+' -f '+args.reads[0]+' '+args.reads[0][:-3]+'.paf.gz > '+args.reads[0][:-3]+'.gfa'
                    command = os.path.join(folder, 'miniasm', 'miniasm')+' '+args.opt_miniasm[0]+' -f '+copyfile+' '+args.reads[0][:-3]+'.paf.gz > '+args.reads[0][:-3]+'.gfa'
                    os.system(command)

                    #reads sorting
                    print("Sorting the reads...")
                    command = os.path.join(folder, 'reads_sorting')+' '+args.reads[0]+' '+args.reads[0][:-3]+'.gfa '+args.reads[0][:-3]+'_sorted.fa '+args.reads[0][:-3]+'.paf.gz '+str(args.ctg_sort[0])
                    os.system(command)
                    print("Sorted reads saved in "+args.reads[0][:-3]+"_sorted.fa.")

                    #removing useless files (gfa, paf.gz, reads copy)
                    print("Removing temporary files...")
                    command = 'rm '+copyfile
                    os.system(command)
                    command = 'rm '+args.reads[0][:-3]+'.gfa'
                    os.system(command)
                    command = 'rm '+args.reads[0][:-3]+'.paf.gz'
                    os.system(command)
                else:
                    print("Wrong number for --ctg_sort. Please write 0, 1 or 2.")
            else:
                os.system("python3 sort_the_reads.py -h")
                print("Please fill the --reads argument.")

if __name__ == '__main__':
    main()
