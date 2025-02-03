#!/usr/bin/env python

import os
import sys
import argparse
import subprocess

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
    parser.add_argument('--output', type=str, nargs=1, required=False, default=["None"],
                        help='Output path. Default=readspath_sorted.fasta')
    parser.add_argument('-k', '--keep', action='count', default=0, required=False, 
                        help='Keep temporary files.')
    parser.add_argument('--rev_comp', type=int, nargs=1, required=False, default=[1],
                        help='Reads are the same in the output and the input, separated by strands (0) or reads are all put in the same strand in the output (1). Default=1')
    parser.add_argument('--ctg_sort', type=int, nargs=1, required=False, default=[1],
                        help='Algorithm used to sort the contigs : random order (0), depth-first search (1), breadth-first search (2). Default=1')
    parser.add_argument('--opt_minimap_ava', type=str, required=False, default="-t32 -k21 -w15",
                        help='String containing all options to run all-vs-all minimap2 (miniasm input). Please write it this way: --opt_minimap_ava="TheOptionsYouWant". Default= -t32 -k21 -w15')
    parser.add_argument('--opt_minimap_reads_vs_ctgs', type=str, required=False, default="-k21 -w15",
                        help='String containing all options to run reads vs contigs minimap2 (miniasm input). Please write it this way: --opt_minimap_reads_vs_ctgs="TheOptionsYouWant". Default= -k21 -w15')
    parser.add_argument('--opt_miniasm', type=str, required=False, default="-I1 -F1",
                        help='String containing all options to run miniasm. Please write it this way: --opt_miniasm="TheOptionsYouWant". Default=-I1 -F1')
    parser.add_argument('-t', '--memtime', action='count', default=0, required=False,
                        help='Creates files with time and memory summary (suffix _memtime).')
    parser.add_argument('-v', '--version', action='count', default=0, required=False,
            help='Prints the used versions of minimap2, miniasm and OReO.')
    args = parser.parse_args()

    script_file = os.path.abspath(sys.argv[0])
    script_dir = os.path.dirname(script_file)
    minimap2_file = os.path.join(script_dir,'minimap2','minimap2')
    miniasm_file = os.path.join(script_dir,'miniasm','miniasm')
    copy_readsfile_file = os.path.join(script_dir,'copy_readsfile')
    reads_sorting_file = os.path.join(script_dir,'reads_sorting')
    log_to_plot_file = os.path.join(script_dir,'log_to_plot.py')

    if len(sys.argv)==1:
        parser.print_help()
    else:
        if args.version>0:
            #afficher les versions de minimap2, miniasm, OReO
            print("Minimap2 version :")
            os.system(minimap2_file + " --version")
            print("Miniasm version :")
            os.system(miniasm_file + " -V")
            print("OReO version:")
            print(oreo_version)
        else:
            if(args.reads!=None):
                if args.ctg_sort[0] in [0,1,2]:
                    #retrieve the folder from where the file is executed
                    reads_file = os.path.abspath(args.reads[0])
                    readscopy_file = os.path.join(os.path.dirname(reads_file),"copie-"+os.path.basename(reads_file))
                    ctg_reads_file = reads_file if args.ctgs_reads[0]=="None" else os.path.abspath(args.ctgs_reads[0])
                    minimap_benchmark = os.path.splitext(reads_file)[0]+'_minimap_all_vs_all_memtime.txt'
                    miniasm_benchmark = os.path.splitext(reads_file)[0]+'_miniasm_memtime.txt'
                    minimap_map_ctgs_benchmark = os.path.splitext(reads_file)[0]+'_minimap_reads_vs_contigs_memtime.txt'
                    ctgs_retrieve_benchmark = os.path.splitext(reads_file)[0]+'_ctgs_retrieve_memtime.txt'
                    sorting_paf_benchmark = os.path.splitext(reads_file)[0]+'_paf_sort_memtime.txt'
                    sorting_benchmark = os.path.splitext(reads_file)[0]+'_sort_memtime.txt'
                    paf_file = os.path.splitext(reads_file)[0] + '.paf.gz'
                    gfa_file = os.path.splitext(reads_file)[0] + '.gfa'
                    gfa_links_file = os.path.splitext(reads_file)[0] + '_links.gfa'
                    ctgs_file = os.path.splitext(reads_file)[0] + '_contigs.fasta'
                    reads_vs_ctgs_paf_file = os.path.splitext(reads_file)[0] + '_reads_vs_ctgs.paf'
                    sorted_reads_vs_ctgs_paf_file = os.path.splitext(reads_file)[0] + '_reads_vs_ctgs_sorted.paf'
                    log_file = os.path.splitext(reads_file)[0] + '_oreo_log.txt'
                    reverse_order_file = os.path.splitext(reads_file)[0] + '_reverse_order.txt'
                    pie_file = os.path.splitext(reads_file)[0] + '_map_reads_pie.png'
                    if args.output[0]=="None":
                        reads_sorted_file = os.path.splitext(reads_file)[0] + '_sorted' + os.path.splitext(reads_file)[1]
                    else:
                        reads_sorted_file = args.output[0]

                    # minimap2 : all-vs-all of the reads
                    print("Mapping the reads...")
                    if args.memtime>0:
                        command = '/usr/bin/time -v '+minimap2_file+' -x ava-'+args.techno[0]+' '+args.opt_minimap_ava+' '+ctg_reads_file+' '+ctg_reads_file+' 2> '+minimap_benchmark+' | gzip -1 > '+paf_file
                    else:
                        command = minimap2_file+' -x ava-'+args.techno[0]+' '+args.opt_minimap_ava+' '+ctg_reads_file+' '+ctg_reads_file+' | gzip -1 > '+paf_file
                    os.system(command)

                    # miniasm to construct contigs using minimap2 output
                    print("Building contigs...")
                    if args.memtime>0:
                        command = '/usr/bin/time -v '+miniasm_file+' '+args.opt_miniasm+' -f '+ctg_reads_file+' '+paf_file+' > '+gfa_file+' 2> '+miniasm_benchmark
                    else:
                        command = miniasm_file+' '+args.opt_miniasm+' -f '+ctg_reads_file+' '+paf_file+' > '+gfa_file
                    os.system(command)

                    # Retrieve miniasm contigs and links
                    print("Retrieving contigs and links...")
                    command = "grep -c '^S' "+gfa_file
                    nb_ctgs_proc = subprocess.run(command, shell=True, check=True, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                    nb_ctgs = int(nb_ctgs_proc.stdout.strip())
                    if args.memtime>0:
                        command = '/usr/bin/time -v awk \'/^S/ {print ">"$2"\\n"$3}\' '+gfa_file+' | fold > '+ ctgs_file +' 2> '+ctgs_retrieve_benchmark
                    else:
                        command = 'awk \'/^S/ {print ">"$2"\\n"$3}\' ' + gfa_file + ' | fold > ' + ctgs_file
                    os.system(command)
                    command = "grep '^L' "+gfa_file+" > "+gfa_links_file
                    os.system(command)

                    # Copie des reads ?
                    print("Copying the reads...")
                    command = copy_readsfile_file + ' ' + reads_file + ' ' + readscopy_file
                    os.system(command)
                    
                    # minimap2 : map the reads on the contigs
                    print("Mapping the reads on the contigs...")
                    if args.memtime>0:
                        command = '/usr/bin/time -v '+minimap2_file+' '+args.opt_minimap_reads_vs_ctgs+' -x map-'+args.techno[0]+' '+ctgs_file+' '+readscopy_file+' 2> '+minimap_map_ctgs_benchmark+' > '+reads_vs_ctgs_paf_file
                    else:
                        command = minimap2_file+' '+args.opt_minimap_reads_vs_ctgs+' -x map-'+args.techno[0]+' '+ctgs_file+' '+readscopy_file+' > '+reads_vs_ctgs_paf_file
                    os.system(command)

                    # Sort minimap output by increasing position of the reads in each contig
                    print("Sorting mapping output...")
                    sort = " -k5,5 " if args.rev_comp[0]==0 else " "
                    if args.memtime>0:
                        command = "/usr/bin/time -v sort"+sort+"-k6,6 -k8,8n "+reads_vs_ctgs_paf_file+' > '+ sorted_reads_vs_ctgs_paf_file +' 2> '+sorting_paf_benchmark
                    else:
                        command = "sort"+sort+"-k6,6 -k8,8n "+reads_vs_ctgs_paf_file+' > '+ sorted_reads_vs_ctgs_paf_file
                    os.system(command)

                    # Final reads sorting using minimap output
                    print("Sorting the reads...")
                    if args.memtime>0:
                        command = '/usr/bin/time -v '+reads_sorting_file + ' ' + reads_file + ' ' + args.format[0] + ' ' + gfa_links_file + ' '+ str(nb_ctgs) + ' ' + reads_sorted_file + ' ' + sorted_reads_vs_ctgs_paf_file + ' ' + log_file + ' '+ str(args.rev_comp[0]) + ' '+ reverse_order_file + ' ' + str(args.ctg_sort[0])+' 2> '+sorting_benchmark
                    else:
                        command = reads_sorting_file + ' ' + reads_file + ' ' + args.format[0] + ' ' + gfa_links_file + ' '+ str(nb_ctgs) + ' ' + reads_sorted_file + ' ' + sorted_reads_vs_ctgs_paf_file + ' ' + log_file + ' '+ str(args.rev_comp[0]) + ' '+ reverse_order_file + ' ' + str(args.ctg_sort[0])
                    os.system(command)
                    print("Sorted reads saved in "+reads_sorted_file)

                    #removing useless files (gfa, paf.gz, reads copy)
                    if args.keep==0:
                        print("Removing temporary files...")
                        command = 'rm '+readscopy_file
                        os.system(command)
                        command = 'rm '+gfa_file
                        os.system(command)
                        command = 'rm '+paf_file
                        os.system(command)
                        command = 'rm '+gfa_links_file
                        os.system(command)
                        command = 'rm '+sorted_reads_vs_ctgs_paf_file
                        os.system(command)
                        command = 'rm '+reads_vs_ctgs_paf_file
                        os.system(command)
                        command = 'rm '+ctgs_file
                        os.system(command)
                        command = 'rm '+log_file
                        os.system(command)
                        

                    # Plot infos in the log file
                    #command = "python3 " + log_to_plot_file + " " + log_file + " " + pie_file
                    #os.system(command)
                else:
                    print("Wrong number for --ctg_sort. Please write 0, 1 or 2.")
            else:
                parser.print_help()
                print("Please fill the --reads argument.")

if __name__ == '__main__':
    main()