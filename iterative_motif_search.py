#takes as input a fasta file with >WEIGHTS as the first line
#generates a random set of seq weights drawn from u(0.1,1), to adjust for uncertainty in p-value calls
#calls meme iteratively to search for motifs in those sequences
#each iteration, counts the number of overlaps at each position of the promoter in the specified sequence ('myseq')

import subprocess
import os
from Bio.Motif.Parsers import MEME

#import modules
import random 
import sys
import time

#changeable parameters
set_size = 14
promoter_length_of_myseq = 105

#unchangeable parameters
weight_min = 1
weight_max = 10

#initialize promoter array as a list
promoter_array = list([0] * promoter_length_of_myseq)

#define sequence of interest
myseq = "GBAA_pXO1_0164"

#number of iterations of algorithm
iterations = 1000

#number of times meme finds a motif in the sequence specified as "myseq"
found = 0

#number of times meme fails to find a motif in the sequence specified as "myseq"
not_found = 0

#get current directory to switch back and forth with output
current_directory = os.getcwd()

#make string of current directory and meme to be able to switch
directory = [current_directory,"/meme_out"]

#define meme directory
output_directory = ''.join(directory)

weights_array = list()

#main for loop for the algorithm, over the number of iterations
#steps: 1) create (quasi-)random seq weights, 2)calls meme from python, 3) parse meme results, 4) update array
for k in range(iterations): 
    
    weights_array[:] = []    
    
    for i in range(set_size - 1): 

        rand = random.randint(weight_min, weight_max)
        rand_weight = rand * 0.1
        weights_array.append(rand_weight)

    #create weight string to add to top of new file
    weights = ' '.join("%s" % i for i in weights_array)
    weights_string = ">WEIGHTS "
    weights_string += weights

    print weights_string
    
    fasta_sequences = "co2atxaseqs.fasta"
    fasta_with_weights = "seqweights.fasta"

    fin = open(fasta_sequences, "r") 
    fout = open(fasta_with_weights, "w") 

    for line in fin:
        fout.write( line.replace (">WEIGHTS", weights_string) )

    fin.close()
    fout.close()
    
    #call command line meme
    subprocess.call(["meme", "-dna", "-revcomp", "seqweights.fasta", "-minw", "8", "-maxw", "14"])

    os.chdir(output_directory)
    f = open("meme.txt")

    record = MEME.read(f)
    for motif in record.motifs:
        for instance in motif.instances:
            if instance.sequence_name == myseq:
                print instance.sequence_name, instance.start, instance.length 
                found += 1
                for i in range(instance.length):
                    position = instance.start+i - 1 #subtract one to adjust for meme indexing from 1 rather than 0
                    promoter_array[position] += 1
            else: 
                print "not found"
                not_found += 1
                
    os.chdir(current_directory)
    
percent_found = found / (found+not_found)        
#can pipe this to your preferred location        
print promoter_array
print percent_found