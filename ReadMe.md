# OReO - Optimized Read Order
Sorting reads to bring k-mers closer
## Compilation

git clone --recurse-submodules https://github.com/girunivlille/sort_the_reads

cd minimap2 && make

cd ..

cd miniasm && make

cd ..

make

## Usage
#### Arguments 
--reads : Fasta file containing the reads to sort

<code>
python3 sort_the_reads.py --reads readsfile.fa
</code>
