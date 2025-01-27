# OReO - Optimized Read Order

OReO take a FASTA/FASTQ file as an input and reorders it according to seuqneces order in the genome (i.e. sequence similarity).

## Compilation

git clone --recurse-submodules https://github.com/girunivlille/oreo
cd minimap2 && make
cd ..
cd miniasm && make
cd ..
make

## Usage
#### Arguments 
--reads : Fasta file containing the reads to sort
<code>
python3 oreo.py --reads readsfile.fa [OPTIONS]
</code>
