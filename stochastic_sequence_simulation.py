random_nucleotide_simulation.py

#generates a text file in fasta format of a given number of sequences of random nucleotides with a given length, 
#a given motif length, a given number of mismatches in each realization of the motif, 
#and a given number of sequences containing the motif
#also allows for a number of motifs at the end (typically one) to have no mismatches in it 

#import modules
import random 
import sys

set_size = 30
length = 100
motif_set_size = 0
motif_no_mismatch_size = 0
motif = list("AAACCGGGTT")
motiflength = len(motif)
number_mismatch = 2

def simulate_sequence_with_motif(length):

    print "inside loop"
    motif = list("AAACCGGGTT")
    dna = ["A", "C", "G", "T"]
    sequence_list = list() #define variable
    for i in range(length):
        random_base = random.choice(dna)
        sequence_list.append(random_base)
    motif_start = random.randint(0, length - motiflength)
    print motif_start
    mod_motif = simulate_mismatch(motif) 
    sequence_list[motif_start:(motiflength+motif_start)] = mod_motif[0:motiflength]
    sequence = "".join(sequence_list)
    return sequence
        
def simulate_sequence_with_motif_no_mismatch(length):

    print "inside loop"
    motif = list("AAACCGGGTT")
    dna = ["A", "C", "G", "T"]
    sequence_list = list() #define variable
    for i in range(length):
        random_base = random.choice(dna)
        sequence_list.append(random_base)
    motif_start = random.randint(0, length - motiflength)
    print motif_start
    sequence_list[motif_start:(motiflength+motif_start)] = motif[0:motiflength]
    sequence = "".join(sequence_list)
    return sequence
    
#generates a specified number of mismatches within a motif    
def simulate_mismatch(motif):    
    
    mismatches = 0
    position = -1
    #mismatchflag = repeat_to_length('0', motiflength) 
    mismatchflag = [0] * motiflength
    print mismatchflag
    while mismatches < number_mismatch: #loop until 2 mismatches created
        print motif
        if position == -1: #flag indicating that a new position needs to be sampled
            position = random.randint(0, motiflength - 1) #sample a position along the motif
            print position
        else: 
            if mismatchflag[position] == 0: #if the base at this position hasn't been changed yet, change it
                dna = ['A', 'C', 'G', 'T'] 
                random_base = random.choice(dna)
                print random_base
                if random_base != motif[position]: #if random base doesn't change the base, then choose base again
                    motif[position] = random_base 
                    mismatchflag[position] = 1 
                    mismatches += 1 
                    print mismatches 
                    position = -1 
            else: 
                position = -1 #if the position has already been changed, choose a new one
    return motif
    
def simulate_sequence(length): 
    dna = ["A", "C", "G", "T"]
    sequence_list = list() #define variable
    for i in range(length):
        random_base = random.choice(dna)
        sequence_list.append(random_base)
    sequence = "".join(sequence_list)
    return sequence    

sequenceset = [] 

#add sequences with no motif 
for i in range(set_size-motif_set_size):

    sequenceset.append(simulate_sequence(length))
   
#add sequences with motif with mismatches  
for i in range(motif_set_size-motif_no_mismatch_size): 

    sequenceset.append(simulate_sequence_with_motif(length))
    
#add sequences (typically just the last one) with motif without any mismatches    
for i in range(motif_no_mismatch_size): 

    sequenceset.append(simulate_sequence_with_motif_no_mismatch(length))
        
filename = "randomseqs.fasta"

f = open(filename,"w") #open only for writing; existing file with same name replaced

f.write(">WEIGHTS")
f.write(" ")
f.write("\n")    

j=0

#numbers and writes sequences to a fasta file
for sequence in sequenceset:
    f.write(">sequence")
    a = str(j)
    f.write(a)
    f.write("\n")
    f.write(sequence)
    f.write("\n")
    j += 1