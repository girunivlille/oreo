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

## **RAM usage**

**If you have limited RAM (or in order to save time), please use the --ctgs_reads option described in the following lines and specify a subset of the reads to construct the contigs.**
If the sample only or mostly contains one specie, 30X of the longest reads can be sufficient. 
If you do not use this option, the tool will automatically use the same file that needs to be sorted to construct the contigs and it might take time and RAM. 


## Usage

`python3 oreo.py --reads readsfile.fasta|readfile.fastq --format fasta|fastq --techno ont|pb [OPTIONS]`

Options:
* `--output_dir` - name of the output directory that will be created (relative path). Default=PathToReadsFile/ReadsFileName_sorted
* `--rev_comp` - reads are in the same orientation in the output and in the input, separated by strands (0) or reads are reverse complemented when necessary in the output (1). Default=1
* `--ctgs_reads` - path of the fasta/fastq file containing the reads that will be used to construct contigs. It is recommanded to use at least 30X of the longest reads. Default= reads file in --reads.
* `--ctg_sort` - algorithm used to sort the contigs : random order (0), depth-first search (1), breadth-first search (2). Default=1
* `--opt_minimap_ava` - string containing all options to run all-vs-all minimap2 (miniasm input). Please write it this way: --opt_minimap_ava="TheOptionsYouWant". Default= -t32 -k21 -w15 for ONT and -t32 -k28 -w100 for HiFi
* `--opt_minimap_reads_vs_ctgs` - string containing all options to run reads vs contigs minimap2 (miniasm input). Please write it this way: --opt_minimap_reads_vs_ctgs="TheOptionsYouWant". Default= -k21 -w15 for ONT and -k28 -w100 for HiFi
* `--opt_miniasm` - string containing all options to run miniasm. Please write it this way: --opt_miniasm="TheOptionsYouWant". Default=-I1 -F1
* `-k, --keep` - keep temporary files (minimap, miniasm)
* `-t, --memtime` - Create a csv file with a time and memory summary for each step (suffix _memtime)
* `-h, --help` - print help

OReO outputs the sorted read file and a file named `readfile`_reverse_order.txt. The second file can be used to retrieve the initial file from the sorted file.
It can be deleted if the original order or original strands will not be needed afterwards.


## Integrity cheking

One can verify that the reads in the sorted file are the same that in the initial file using the following command:

`python3 integrity_checker.py --readfiles readfile.fasta,readfile_sorted.fasta [OPTIONS]`

Options:
* `--seq_only` - compares only sequences. Will ignore quality scores if one or several files are fastq.
* `--rc_sensitivity` - sensitivity to reverse complement. If used, a sequence and its reverse complement will not be handled as the same sequence. By default they are handled as the same.

This can be extended to more than two files. By default, this checker considers that a sequence holds the same information as its reverse complement. This can be disabled by the --rc_sensitivity option.
The checker exits with 1 if at least 2 files hold different content, else it exits with 0.

## Retrieve the original order

The original file can be retrieved from the sorted file and the `readfile`_reverse_order.txt file using the command:

`./reverse_order readfile_sorted.fasta|readfile_sorted.fastq fasta|fastq readfile_reverse_order.txt outputfile.fasta|outputfile.fastq`

It retrieves both the original order and the original strands.

## Example

An example is available in the directory `example`. The reads were simulated with Badread from a phage genome.
To try OReO, one can use the following lines:

```
cd oreo
python3 oreo.py --reads example/badread_lambdavirus_ont_50X.fasta --format fasta --techno ont --opt_minimap_ava='-t10 -k15 -w10' --opt_miniasm='-1 -2' --opt_minimap_reads_vs_ctgs='-k15 -w10'
python3 integrity_checker.py --readfiles example/badread_lambdavirus_ont_50X.fasta,example/badread_lambdavirus_ont_50X_sorted.fasta
./reverse_order example/badread_lambdavirus_ont_50X_sorted.fasta fasta example/badread_lambdavirus_ont_50X_reverse_order.txt example/badread_lambdavirus_ont_50X_desorted.fasta
```
