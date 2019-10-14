# Screens a list of .seq files for barcodes and other relevant information

# IMPORT MODULES
from Bio import SeqIO
from Bio.Seq import Seq
import csv

# FUNCTION DEFINITIONS
def oligo_info (oligo):
    'Takes an string variable "oligo" as an argument and returns list of 2 values (oligo location, oligo orientation)'
    # declare variables
    oligo_pos = 0		# an integer storing oligo position
    oligo_ori = 0		# an integer storing oligo orientation
    pos_ori = [] 		# initialize a list of 2 values (oligo position, oligo orientation)
    # start function code
    oligo_pos = dnaseq.find(oligo)	# finds location of "oligo" within "dnaseq" and returns position if found, -1 if not found
    if oligo_pos != -1:				# a conditional for when oligo is found in dnaseq
        oligo_ori = 1				# stores oligo orientation
        pos_ori = [oligo_pos, oligo_ori]	# stores oligo position and orientation
    else:
        rc_oligo = Seq(oligo).reverse_complement()	# a string storing the reverse-complement sequence of oligo
        oligo_pos = dnaseq.find(rc_oligo)	# finds location of rev-com oligo within dnaseq and returns position if found, -1 if not found
        if oligo_pos != -1:					# a conditional for when reverse-complement oligo is found in dnaseq
            oligo_ori = -1					# stores oligo orientation
            pos_ori = [oligo_pos, oligo_ori]	# stores oligo position and orientation
        else:
            pos_ori = ['NA', 'NA']			# the output if neither forward or reverse orientation of oligo is found in dnaseq
    print(pos_ori)
    return pos_ori
def barcode_seq (primer_info_list):
    'Takes a list of information regarding poition and orientation of primers ("primer_info_list") and returns a string storing barcode sequence'
    # declare variables
    bar_start = 0 #Defines the barcode starting position as numerical value
    barcode = '' #Defines the barcode itself as a string 
	# start function code
    if primer_info_list[1] == 1:
        bar_start = primer_info_list[0]+len(primer1)
        barcode = dnaseq[bar_start:(bar_start+20)]
    elif primer_info_list[1] == -1:
        bar_start = primer_info_list[0]
        barcode = (dnaseq[bar_start-20:(bar_start)]).reverse_complement()
    elif primer_info_list[3] == 1:
        bar_start = primer_info_list[2]-20
        barcode = dnaseq[bar_start:(bar_start+20)]
    elif primer_info_list[3] == -1:
        bar_start = primer_info_list[2]+len(primer2)+20
        barcode = (dnaseq[bar_start-20:(bar_start)]).reverse_complement()
    else:
        barcode = 'NA'
    print(barcode)
    return barcode
    
# START OF MAIN FUNCTION
primer1 = 'GATGTCCACGAGGTCTCT'		# a string storing nucleotide sequence of one primer
primer2 = 'CGTACGCTGCAGGTCGAC'	# a string storing nucleotide sequence of another primer
myfilenames = input("Specify filenames WITHOUT extension (.seq is assumed), separated only by commas: \n").split(',') # get user input
# start csv file
with open('find_barcode_results.csv', mode='w') as barcode_file:
    fieldnames = ['filename', 'primer1_pos', 'primer1_ori', 'primer2_pos', 'primer2_ori', 'barcode_loc', 'barcode_ori', 'barcode_seq', 'barcode_len']
    barcode_writer = csv.writer(barcode_file)
    barcode_writer.writerow(fieldnames)
    # write info to line in csv file
    for myfilename in myfilenames:			# perform loop recursively on all files listed by the user
        seq_filename = str(myfilename) + ".seq"
        input_handle  = open(seq_filename, "r")
        for seq_record in SeqIO.parse(input_handle, "fasta") :
            dnaseq = seq_record.seq			# a string variable to store the dna sequence from individual file
            # build output list
            output_list = []		# initialize a list of values to be added to csv file (one line of csv)
            output_list.append(myfilename)			# first value in list is filename
            output_list.extend( oligo_info(primer1) )	# add primer 1 info to list
            output_list.extend( oligo_info(primer2) )	# add primer 2 info to list
            primer_info_list1 = [output_list[1],output_list[2],output_list[3],output_list[4]] # stores a list of neccesary input for barcode_seq() function
            print(primer_info_list1)
            barcode = str(barcode_seq (primer_info_list1))	# assign barcode sequence to string variable "barcode"
            output_list.extend( oligo_info(barcode) )	# add barcode info to list
            output_list.append(barcode)			# add barcode sequence to list
            output_list.append(len(barcode) )		# add barcode length to list
            #write info to csv file
            barcode_writer.writerow(output_list)  	
    input_handle.close()	# close single file being read
barcode_file.close()	# close newly-written csv file

# print closing message
print ('CSV file created called \"find_barcode_results.csv\"')


