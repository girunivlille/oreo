# OReO: Optimised Read Order

OReO takes a FASTA/FASTQ file as an input and reorders it according to sequences order in the genome (i.e. sequence similarity).

## Installation 

```
git clone --recurse-submodules https://github.com/girunivlille/oreo
cd oreo
cd minimap2 && make
cd ..
cd miniasm && make
cd ..
make
```

## Usage

`python3 oreo.py --reads readsfile.fa --format fasta|fastq --techno ont|pb [OPTIONS]`

Options:
* `--output` - name of the sorted read file
* `--rev_comp` - reads are in the same orientation in the output and in the input, separated by strands (0) or reads are reverse complemented when necessary in the output (1). Default=1
* `--ctgs_reads` - path of the fasta/fastq file containing the reads that will be used to construct contigs. Default= reads file in --reads.
* `--ctg_sort` - algorithm used to sort the contigs : random order (0), depth-first search (1), breadth-first search (2). Default=1
* `--opt_minimap_ava` - string containing all options to run all-vs-all minimap2 (miniasm input). Please write it this way: --opt_minimap_ava="TheOptionsYouWant". Default= -t32 -k21 -w15
* `--opt_minimap_reads_vs_ctgs` - string containing all options to run reads vs contigs minimap2 (miniasm input). Please write it this way: --opt_minimap_reads_vs_ctgs="TheOptionsYouWant". Default= -k21 -w15
* `--opt_miniasm` - string containing all options to run miniasm. Please write it this way: --opt_miniasm="TheOptionsYouWant". Default=-I1 -F1
* `-k, --keep` - keep temporary files (minimap, miniasm)
* `-t, --memtime` - store '_memtime' files containing each step's memory and time usage.
* `-h, --help` - print help

## Integrity cheking

One can verify that the reads in the sorted file are the same that in the initial file using the following command:

`python3 integrity_checker.py --readfiles readfile.fasta,readfile_sorted.fasta [OPTIONS]`

Options:
*
