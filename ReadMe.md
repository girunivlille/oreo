# OReO - Optimized Read Order
Sorting reads to bring k-mers closer
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
