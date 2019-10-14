# Screens a list of fasta files for hits
# cd "C:/Users/Angela/Documents/KHR Project/Seq Files"
# transpose column of gb names in excel, then save as .csv and open in notepad


#1. import modules
from Bio import SeqIO
from Bio.Seq import Seq
import re

#3.5 assigning variables
# primers:
#primer1 = 'GATGTCCACGAGGTCTCT'
#primer2 = 'CGTACGCTGCAGGTCGAC'
#ref_seq = 'GAAAGTTACGATGTCCACGAGGTCTCTGGTTAATGAGA'
myfile = str('PR_604115-3_A3__C10.seq')
input_handle = open(myfile, 'r')
#2. make fxn for finding primer1 location
#Finds location of pirmer 1
def find_p1r_loc (ref_seq):
    p1r_loc = 0
    primer1r = 'AGAGACCTCGTGGACATC'
    p1r_loc = ref_seq.find(primer1r)
    return p1r_loc 
#Finds location of pirmer 2
def find_p2r_loc (ref_seq):
    p1r_loc = 0
    primer2r = 'GTCGACCTGCAGCGTACG'
    p2r_loc = ref_seq.find(primer2r)
    return p2r_loc 
#Finds location of pirmer 1F
def find_p1f_loc (ref_seq):
    p1f_loc = 0
    primer1f = 'GATGTCCACGAGGTCTCT'
    p1f_loc = ref_seq.find(primer1f)
    return p1f_loc 
#Finds location of pirmer 2F
def find_p2f_loc (ref_seq):
    p1f_loc = 0
    primer2f = 'CGTACGCTGCAGGTCGAC'
    p2f_loc = ref_seq.find(primer2f)
    return p2f_loc

#Defines the file as a fasta so the seq can be read
for seq_record in SeqIO.parse(input_handle, "fasta") :
    dnaseq = seq_record.seq
    print(dnaseq)
#Defines the functions above as values that can be printed
    yo_ma = find_p1r_loc(dnaseq)
    print(yo_ma)
    yo_yo_ma = find_p2r_loc(dnaseq)
    print(yo_yo_ma)
    yo_yo_yo_ma = find_p1f_loc(dnaseq)
    print(yo_yo_yo_ma)
    yo_yo_yo_yo_ma = find_p2f_loc(dnaseq)
    print(yo_yo_yo_yo_ma)
input_handle.close()


#3. Make fxn for finding primer orientation
#def find_primer_orientation ()


##
###4. Find sequence between primers
###5. Get user input for list of file names for searching 
###6. perform search on files 
###7. Output CSV with file name, primer 1 orientation, primer 2 orientation, sequence (btw primers) length, sequence (btw primers) output
##
##myfilenames = input("Specify saved fasta-formatted genbank ids, separated only by commas:\n").split(',')
##newfile = open("findKP4Genes.csv", "w")
##newfile.write("GenBankID,IRR1_loc,IRR1_orient,SSM4_loc,SSM4_orient,KHR_loc,KHR_orient\n")
##khrlist = []
##
##for myfilename in myfilenames:
##    fa_filename = str(myfilename) + ".fasta"
##    input_handle  = open(fa_filename, "r")
##    for seq_record in SeqIO.parse(input_handle, "fasta") :
##        dnaseq = seq_record.seq
##        full_list = []
##        output_list = []
##
##        # scan
##        full_list.extend(find_hits (IRR1f))
##        full_list.extend(find_hits (SSM4f))
##        full_list.extend(find_hits (KHRf))
##        print(myfilename + ' full list: ' + str(full_list))
##            
##        #write info to line in file
##        newfile.write(myfilename + ',' + str(full_list[0]) + ',' + str(full_list[1]) + ',' + str(full_list[2]) + ',' + str(full_list[3]) + ','  + str(full_list[4]) + ','  + str(full_list[5]) + '\n')
##        if full_list[5] != 'NA':
##            khrlist.append(myfilename)
##    input_handle.close()
##
##print ('Sequences positive for KHR: ')
##print (khrlist)
##
##newfile.close()
##
##
