# python
# find the longest open reading frame of mRNA

import sys

def translateORF(sequence, start, end):
	
	codonTable = {
		"UUU":"F", "UUC":"F", "UUA":"L", "UUG":"L",
		"UCU":"S", "UCC":"S", "UCA":"S", "UCG":"S",
		"UAU":"Y", "UAC":"Y", "UAA":"STOP", "UAG":"STOP",
		"UGU":"C", "UGC":"C", "UGA":"STOP", "UGG":"W",
		"CUU":"L", "CUC":"L", "CUA":"L", "CUG":"L",
		"CCU":"P", "CCC":"P", "CCA":"P", "CCG":"P",
		"CAU":"H", "CAC":"H", "CAA":"Q", "CAG":"Q",
		"CGU":"R", "CGC":"R", "CGA":"R", "CGG":"R",
		"AUU":"I", "AUC":"I", "AUA":"I", "AUG":"M",
		"ACU":"T", "ACC":"T", "ACA":"T", "ACG":"T",
		"AAU":"N", "AAC":"N", "AAA":"K", "AAG":"K",	
		"AGU":"S", "AGC":"S", "AGA":"R", "AGG":"R",	
		"GUU":"V", "GUC":"V", "GUA":"V", "GUG":"V",
		"GCU":"A", "GCC":"A", "GCA":"A", "GCG":"A",
		"GAU":"D", "GAC":"D", "GAA":"E", "GAG":"E",
		"GGU":"G", "GGC":"G", "GGA":"G", "GGG":"G",} # dictionary of RNA codons 
    
	protein = []
#	print start, end
	for i in range(start,end,3):
        	codon = sequence[i:i+3]
#		print codon
        	if codon in codonTable:
            		aminoacid = codonTable[codon]
            		protein.append(aminoacid)
 #       	else:
 #           		protein.append("N")
	return "".join(protein)


def getORF(sequence):
	start_codon = 'AUG'
	stop_codon = ['UAA', 'UGA', 'UAG']
	start_codon_index = 0
	end_codon_index = 0
	start_codon_found = False
	
#	print("length of sequence is" + str(len(sequence)))
 
    	orfs = []
	orflength = 0
 
    	for j in range(0, 3):
#		print ("frame" + str(j+1))
        	for indx in range(j, len(sequence), 3):
            		current_codon = sequence[indx:indx+3]
			#print current_codon
            		if current_codon in start_codon and not start_codon_found:
                		start_codon_found = True
                		start_codon_index = indx
            		if current_codon in stop_codon and start_codon_found:
                		end_codon_index = indx
                		length = end_codon_index - start_codon_index + 1
			
#				print "orflength "+str(orflength)
#				print length, start_codon_index, end_codon_index
	               		if length > orflength:
					frame = j+1
					orf_start = start_codon_index
					orf_end = end_codon_index
					orflength = length
#                		start_codon_found = False
 
        	start_codon_index = 0
        	end_codon_index = 0
        	start_codon_found = False
 
#    	print len(orfs), orfs
#	print orflength, orf_start, orf_end
#	print ("frame" + str(frame))
	return orf_start, orf_end


infile = sys.argv[1]

f = open(infile, "r")
file = f.readlines()

sequences = []
seq = ""
seqId = []

for f in file:
	if not f.startswith('>'):
		f = f.replace(" ", "")      # remove all spaces and newline from the text 
        	f = f.replace("\n", "")
        	seq = seq + f               # ... then form a long sequence

	else:	
		sequences.append(seq)
        	name = f.replace("\n", "")
		seqId.append(name)
		seq = ""

sequences.append(seq)
sequences = sequences[1:]

#	print(name)
#	print(seq)
#	getORF(seq)

# print(seqId)

for i in range(0,len(sequences)):
	print(seqId[i])
#	print(sequences[i])
	start_pos, end_pos = getORF(sequences[i])
	aaseq = translateORF(sequences[i], start_pos, end_pos)
	
	print aaseq


