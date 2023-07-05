git clone --recurse-submodules https://github.com/girunivlille/sort_the_reads

cd minimap2 && make  

cd ..

cd miniasm && make 

cd ..

make
