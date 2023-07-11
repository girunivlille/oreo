import os
import argparse
import time

def main():
    parser = argparse.ArgumentParser(description='Returns a new file containing the sorted reads.')

    parser.add_argument('--reads', type=str, nargs=1, required=True,
                        help='Path of the fasta file containing the reads that have to be sorted.')
    args = parser.parse_args()
    read_filename = args.reads[0]


    read_filename_copy = os.path.splitext(read_filename)[0]+"_copy"+os.path.splitext(read_filename)[1]
    read_filename_paf = os.path.splitext(read_filename)[0] + '.paf.gz'
    read_filename_gfa = os.path.splitext(read_filename)[0] + '.gfa'
    read_filename_sorted = os.path.splitext(read_filename)[0] + '_sorted' + os.path.splitext(read_filename)[1]

    t1 = time.time()


    # copying the reads
    print("Copying the reads...")
    fread = open(read_filename,"r")
    fread_copy = open(read_filename_copy,"w")
    i = 0
    line = fread.readline()
    while line:
        if line[0] == ">":
            fread_copy.write(">."+str(i)+"\n")
            i += 1
        else:
            fread_copy.write(line)
        line = fread.readline()
    fread.close()
    fread_copy.close()



    t2 = time.time()
    print("Done in "+str(round(t2-t1,2))+" s.")

    #minimap2



    print("Mapping the reads...")
    command = 'minimap2 -x ava-pb -t16 -k28 -w15 '+read_filename_copy+' '+read_filename_copy+' | gzip -1 > '+read_filename_paf
    os.system(command)

    t3 = time.time()
    print("Done in "+str(round(t3-t2,2))+" s.")


    #miniasm
    print("Building contigs...")
    command = 'miniasm -I1 -F1 -f '+read_filename+' '+read_filename_paf+' > '+read_filename_gfa
    os.system(command)

    t4 = time.time()
    print("Done in "+str(round(t4-t3,2))+" s.")


    #reads sorting
    print("Sorting the reads...")
    command = './reads_sorting '+read_filename+' '+read_filename_gfa+' '+read_filename_sorted+' '+read_filename_paf
    os.system(command)
    print("Sorted reads saved in "+read_filename_sorted)

    t5 = time.time()
    print("Done in "+str(round(t5-t4,2))+" s.")


    #removing useless files (gfa, paf.gz, reads copy)
    print("Removing temporary files...")
    command = 'rm '+read_filename_copy
    os.system(command)
    command = 'rm '+read_filename_gfa
    os.system(command)
    command = 'rm '+read_filename_paf
    os.system(command)

    t6 = time.time()
    print("Done in "+str(round(t6-t5,2))+" s.")


if __name__ == '__main__':
    main()
