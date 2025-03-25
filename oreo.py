#!/usr/bin/env python

import os
import sys
import argparse
import subprocess
import pandas as pd

oreo_version = 1.0

def main():
    parser = argparse.ArgumentParser(description='OReO : Returns a new file containing the sorted reads.')

    parser.add_argument('--reads', type=str, nargs=1, required=True,
                        help='[MANDATORY] Path of the fasta/fastq file containing the reads that have to be sorted.')
    parser.add_argument('--format', type=str, nargs=1, required=True, default=["fastq"],
                        help='[MANDATORY] Format of the reads file: fasta or fastq. Default=fastq')
    parser.add_argument('--techno', type=str, nargs=1, required=True, default=["ont"],
                        help='[MANDATORY] Sequencing technology : ont (ont), pacbio (pb). Default=ont')
    parser.add_argument('--ctgs_reads', type=str, nargs=1, required=False, default=["None"],
                        help='Path of the fasta/fastq file containing the reads that will be used to construct contigs. Default= reads file in --reads.')
    parser.add_argument('--output_dir', type=str, nargs=1, required=False, default=["None"],
                        help='Name of the output directory that will be created (relative path). Default=PathToReadsFile/ReadsFileName_sorted')
    parser.add_argument('-k', '--keep', action='count', default=0, required=False, 
                        help='Keep temporary files.')
    parser.add_argument('--rev_comp', type=int, nargs=1, required=False, default=[1],
                        help='Reads are the same in the output and the input, separated by strands (0) or reads are all put in the same strand in the output (1). Default=1')
    parser.add_argument('--ctg_sort', type=int, nargs=1, required=False, default=[1],
                        help='Algorithm used to sort the contigs : random order (0), depth-first search (1), breadth-first search (2). Default=1')
    parser.add_argument('--opt_minimap_ava', type=str, required=False, default="None",
                        help='String containing all options to run all-vs-all minimap2 (miniasm input). Please write it this way: --opt_minimap_ava="TheOptionsYouWant". Default= -t32 -k21 -w15 for ONT and -t32 -k28 -w100 for HiFi')
    parser.add_argument('--opt_minimap_reads_vs_ctgs', type=str, required=False, default="None",
                        help='String containing all options to run reads vs contigs minimap2 (miniasm input). Please write it this way: --opt_minimap_reads_vs_ctgs="TheOptionsYouWant". Default= -k21 -w15 for ONT and -k28 -w100 for HiFi')
    parser.add_argument('--opt_miniasm', type=str, required=False, default="-I1 -F1",
                        help='String containing all options to run miniasm. Please write it this way: --opt_miniasm="TheOptionsYouWant". Default=-I1 -F1')
    parser.add_argument('-t', '--memtime', action='count', default=0, required=False,
                        help='Create a csv file with a time and memory summary for each step (suffix _memtime).')
    parser.add_argument('-v', '--version', action='count', default=0, required=False,
            help='Prints the used versions of minimap2, miniasm and OReO.')
    args = parser.parse_args()

    # Retrieve executables paths for oreo, minimap and miniasm
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
            # Display minimap2, miniasm, OReO versions
            print("Minimap2 version :")
            os.system(minimap2_file + " --version")
            print("Miniasm version :")
            os.system(miniasm_file + " -V")
            print("OReO version:")
            print(oreo_version)
        else:
            if(args.reads!=None):
                if args.ctg_sort[0] in [0,1,2]:
                    # Retrieve the folder from where the file is executed
                    reads_file = os.path.abspath(args.reads[0])

                    # Check if the output directory can be created
                    out_dir_absolute_path = os.path.splitext(reads_file)[0] + '_sorted' if args.output_dir[0]=="None" else os.path.join(os.getcwd(), args.output_dir[0])
                    if os.path.exists(out_dir_absolute_path) and os.path.isdir(out_dir_absolute_path):
                        print("Error: the output directory already exists. Please consider deleting it, renaming it, or specifying an other name for the output directory with the argument --output_dir.")
                        sys.exit(1)
                    else:
                        os.makedirs(out_dir_absolute_path)
                    
                    # Name temporary output files put in the output directory
                    base_output_path = os.path.join(out_dir_absolute_path, os.path.splitext(os.path.basename(reads_file))[0])
                    readscopy_file = os.path.join(out_dir_absolute_path, "copie-"+os.path.basename(reads_file))
                    paf_file = base_output_path + '.paf.gz'
                    gfa_file = base_output_path + '.gfa'
                    gfa_links_file = base_output_path + '_links.gfa'
                    ctgs_file = base_output_path + '_contigs.fasta'
                    reads_vs_ctgs_paf_file = base_output_path + '_reads_vs_ctgs.paf'
                    sorted_reads_vs_ctgs_paf_file = base_output_path + '_reads_vs_ctgs_sorted.paf'

                    # Retrieve the path of the reads that are used to construct the contigs
                    ctg_reads_file = reads_file if args.ctgs_reads[0]=="None" else os.path.abspath(args.ctgs_reads[0])

                    # Name output files put in the output directory (sorted reads and positional informations)
                    reads_sorted_file = base_output_path + '_sorted' + os.path.splitext(reads_file)[1]
                    reverse_order_file = base_output_path + '_reverse_order.txt'

                    # Name and initialize the file containing RAM and time informations for each step
                    oreo_ram_time = base_output_path + '_oreo_memtime.csv'
                    if args.memtime>0:
                        with open(oreo_ram_time, "w") as csv_memtime:
                            csv_memtime.write("step,time,memory\n")

                    # Set minimap and miniasm options
                    opt_minimap_ava = "-t10 -k21 -w15" if (args.opt_minimap_ava=="None" and args.techno[0]=="ont") else ("-t10 -k28 -w100" if (args.opt_minimap_ava=="None" and args.techno[0]=="pb") else args.opt_minimap_ava)
                    opt_minimap_reads_vs_ctgs = "-k21 -w15" if (args.opt_minimap_reads_vs_ctgs=="None" and args.techno[0]=="ont") else ("-k28 -w100" if (args.opt_minimap_reads_vs_ctgs=="None" and args.techno[0]=="pb") else args.opt_minimap_reads_vs_ctgs)

                    memtime_separator = [',', ',', ''] if args.memtime>0 else [':time=', 's,max_RAM=', 'k']
                    memtime_output = " -a -o "+oreo_ram_time if args.memtime>0 else ''

                    # minimap2 : all-vs-all of the reads
                    print("\n\n***************Mapping the reads (step 1/7)...")
                    command = '/usr/bin/time -f minimap_ava'+ memtime_separator[0] +'%U'+memtime_separator[1]+'%M'+memtime_separator[2]+ memtime_output +' '+minimap2_file+' -x ava-'+args.techno[0]+' '+opt_minimap_ava+' '+ctg_reads_file+' '+ctg_reads_file + ' | gzip -1 > '+paf_file
                    print(command+"\n\n")
                    os.system(command)

                    # miniasm to construct contigs using minimap2 output
                    print("\n\n***************Building contigs (step 2/7)...")
                    command = '/usr/bin/time -f miniasm'+ memtime_separator[0] +'%U'+memtime_separator[1]+'%M'+memtime_separator[2]+ memtime_output +' ' +miniasm_file+' '+args.opt_miniasm+' -f '+ctg_reads_file+' '+paf_file+' > '+gfa_file
                    print(command+"\n\n")
                    os.system(command)

                    # Retrieve miniasm contigs and links between contigs
                    print("\n\n***************Retrieving contigs and links (step 3/7)...")
                    command = "grep -c '^S' "+gfa_file
                    print(command)
                    nb_ctgs_proc = subprocess.run(command, shell=True, check=True, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                    nb_ctgs = int(nb_ctgs_proc.stdout.strip())
                    command = '/usr/bin/time -f ctgs_retrieve'+ memtime_separator[0] +'%U'+memtime_separator[1]+'%M'+memtime_separator[2]+ memtime_output +' awk \'/^S/ {print ">"$2"\\n"$3}\' '+gfa_file+' | fold > '+ ctgs_file
                    print(command)
                    os.system(command)
                    command = '/usr/bin/time -f links_retrieve'+ memtime_separator[0] +'%U'+memtime_separator[1]+'%M'+memtime_separator[2]+ memtime_output + " grep '^L' "+gfa_file+" > "+gfa_links_file
                    print(command+"\n\n")
                    os.system(command)

                    # Copy reads for the mapping step
                    print("\n\n***************Copying the reads (step 4/7)...")
                    command = '/usr/bin/time -f reads_copy'+ memtime_separator[0] +'%U'+memtime_separator[1]+'%M'+memtime_separator[2]+ memtime_output+' '+copy_readsfile_file + ' ' + reads_file + ' ' + readscopy_file
                    print(command+"\n\n")
                    os.system(command)
                    
                    # minimap2 : map the reads on the contigs
                    print("\n\n***************Mapping the reads on the contigs (step 5/7)...")
                    command = '/usr/bin/time -f minimap_reads_vs_ctgs'+ memtime_separator[0] +'%U'+memtime_separator[1]+'%M'+memtime_separator[2]+ memtime_output+ ' ' +minimap2_file + ' ' + opt_minimap_reads_vs_ctgs + ' -x map-' + args.techno[0] + ' ' + ctgs_file + ' ' + readscopy_file + ' > '+reads_vs_ctgs_paf_file
                    print(command+"\n\n")
                    os.system(command)

                    # Sort minimap output by increasing position of the reads in each contig
                    print("\n\n***************Sorting mapping output (step 6/7)...")
                    sort = " -k5,5 " if args.rev_comp[0]==0 else " "
                    command = "/usr/bin/time -f paf_sort"+ memtime_separator[0] +'%U'+memtime_separator[1]+'%M'+memtime_separator[2]+ memtime_output+" sort"+sort+"-k6,6 -k8,8n "+reads_vs_ctgs_paf_file+' > '+ sorted_reads_vs_ctgs_paf_file 
                    print(command+"\n\n")
                    os.system(command)

                    # Final reads sorting using minimap output
                    print("\n\n***************Sorting the reads (step 7/7)...")
                    command = '/usr/bin/time -f reads_sort'+ memtime_separator[0] +'%U'+memtime_separator[1]+'%M'+memtime_separator[2]+ memtime_output+' '+reads_sorting_file + ' ' + reads_file + ' ' + args.format[0] + ' ' + gfa_links_file + ' '+ str(nb_ctgs) + ' ' + reads_sorted_file + ' ' + sorted_reads_vs_ctgs_paf_file + ' '+ str(args.rev_comp[0]) + ' '+ reverse_order_file + ' ' + str(args.ctg_sort[0])
                    print(command+"\n\n")
                    os.system(command)
                    print("Sorted reads saved in "+reads_sorted_file)

                    if args.memtime>0:
                        df = pd.read_csv(oreo_ram_time)
                        total_memory = df['memory'].max()
                        total_time = df['time'].sum()

                        with open(oreo_ram_time, "a") as csv_memtime:
                            csv_memtime.write(f"total (sum, max),{total_time},{total_memory}\n")


                    # Remove useless files (gfa, paf.gz, reads copy)
                    if args.keep==0:
                        print("\n\nRemoving temporary files...")
                        for file in [readscopy_file, gfa_file, paf_file, gfa_links_file, sorted_reads_vs_ctgs_paf_file, reads_vs_ctgs_paf_file, ctgs_file]:
                            os.remove(file)

                else:
                    print("Wrong number for --ctg_sort. Please write 0, 1 or 2.")
            else:
                parser.print_help()
                print("Please fill the --reads argument.")

if __name__ == '__main__':
    main()