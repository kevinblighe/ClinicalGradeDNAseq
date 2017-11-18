
#Name:			Kevin Blighe
#Email:			kevinblighe@outlook.com / kevin@clinicalbioinformatics.co.uk
#Date:			20th July 2017
#Function(s):		Splits multi-allelic sites; sets unique ID; modifies indels; marks uninterpretable variants as LowQual
#Optional function(s):	Changes ID field in VCF to CHR:POS:REF:VAR:ZYGOSITY

import sys
import os

#Check the number of command line arguments
if not len(sys.argv)==4:
	print "\nError:\tincorrect number of command-line arguments"
	print "Syntax:\tFixInDels.py [Input VCF] [Output VCF] [Overwrite ID field with CHR:POS:REF:VAR:ZYGOSITY? TRUE/FALSE]\n"
	sys.exit()

if sys.argv[1]==sys.argv[2]:
	#boolOverwrite = raw_input("\nInput file is the same as the output file - overwrite input file (sim/nao)?\n")
	print "Error:\tInput file is the same as the output file - choose a different output file\n"
	sys.exit()

if not (sys.argv[3].upper()=="TRUE") and not (sys.argv[3].upper()=="FALSE"):
	print("\nError:\tArgument 4 must be TRUE/FALSE\n")
	sys.exit()

#File input
fileInput = open(sys.argv[1], "r")

#File output for split multi-alleles
fileOutputTemp = open("temp.vcf.tmp", "w")

#Loop through each line in the input file, and split multiallelic sites
print ("Splitting multi-allelic sites...")
for strLine in fileInput:
	#Strip the endline character from each input line
	strLine = strLine.rstrip("\n")

	#The '#' character in VCF format indicates that the line is a header. Ignore these and just output to the new file
	if strLine.startswith("#"):
		fileOutputTemp.write(strLine + "\n")
	else:
		#Split the tab-delimited line into an array
		strArray = [splits for splits in strLine.split("\t") if splits is not ""]

		#Check first if it's multiallelic
		#Multi-allelic variants will have 2 calls in the VAR field, separated by a comma (',')
		if "," in strArray[4]:
			strVars = [splits for splits in strArray[4].split(",") if splits is not ""]
			iNumMultialleles = len(strVars)

			for i in range(0, (iNumMultialleles)):
				strArray[4] = strVars[i]
				fileOutputTemp.write("\t".join(strArray) + "\n")
		else:
			fileOutputTemp.write("\t".join(strArray) + "\n")
print ("Done.")

#Close the files
fileInput.close()
fileOutputTemp.close()

###

#Use the temporary file with split multi-allelic sites as new input
fileInputTemp = open("temp.vcf.tmp", "r")

#File output
fileOutput = open(sys.argv[2], "w")

#Loop through each line in the input file
print "Fixing insertion and deletion (indel) annotation..."
for strLine in fileInputTemp:

	#A boolan value to say whether the record is a CN event from PCR or not
	boolCN = "FALSE"

	strZygosity = "NA"

	#Strip the endline character from each input line
	strLine = strLine.rstrip("\n")

	#The '#' character in VCF format indicates that the line is a header. Ignore these and just output to the new file
	if strLine.startswith("#"):
		fileOutput.write(strLine + "\n")
	else:
		#Split the tab-delimited line into an array
		strArray = [splits for splits in strLine.split("\t") if splits is not ""]

		#Determine the number of characters in the reference base call
		iRefBaseLength = len(strArray[3])

		#Determine the number of characters in the variant base call
		iVarBaseLength = len(strArray[4])

		#If the length of the referance base is greater than 1, it's a deletion in VCF format
		if iRefBaseLength>1:
			#Remove the first base from the reference base call
			strNewRefBase = strArray[3][(1):]

			#Change the old reference base to the new base call
			strArray[4] = strNewRefBase

			#Change the variant to 'del' + reference bases that were deleted
			strArray[4] = strArray[4] = "-"

			#Change the base-position by 1 base (subtraction)
			strArray[1] = str(int(strArray[1]) - 1)

			#If the user requested to change the ID field with CHR:POS:REF:VAR:ZYGOSITY
			if sys.argv[3].upper()=="TRUE":
				#Determine zygosity
				strArrayGenotypeInfo = [splits for splits in strArray[9].split(":") if splits is not ""]

				if strArrayGenotypeInfo[0]=="0/1" or strArrayGenotypeInfo[0]=="1/0":
					strZygosity = "Het"
				elif strArrayGenotypeInfo[0]=="1/1":
					strZygosity = "Hom"
				#If the genotype is missing, look further at the downsampled reads
				#75% reads
				elif strArrayGenotypeInfo[0]=="0/0":
					#Determine zygosity
	                                strArrayGenotypeInfo_75pcReads = [splits for splits in strArray[10].split(":") if splits is not ""]

					if strArrayGenotypeInfo_75pcReads[0]=="0/1" or strArrayGenotypeInfo_75pcReads[0]=="1/0":
	                                        strZygosity = "Het"
        	                        elif strArrayGenotypeInfo_75pcReads[0]=="1/1":
	                                        strZygosity = "Hom"
					#50% reads
					elif strArrayGenotypeInfo[0]=="0/0":
						#Determine zygosity
        	                                strArrayGenotypeInfo_50pcReads = [splits for splits in strArray[11].split(":") if splits is not ""]

                	                        if strArrayGenotypeInfo_50pcReads[0]=="0/1" or strArrayGenotypeInfo_50pcReads[0]=="1/0":
                        	                        strZygosity = "Het"
                                	        elif strArrayGenotypeInfo_50pcReads[0]=="1/1":
                                        	        strZygosity = "Hom"
						#25% reads
						elif strArrayGenotypeInfo[0]=="0/0":
	                                                #Determine zygosity
							strArrayGenotypeInfo_25pcReads = [splits for splits in strArray[12].split(":") if splits is not ""]
							if strArrayGenotypeInfo_25pcReads[0]=="0/1" or strArrayGenotypeInfo_25pcReads[0]=="1/0":
                                                        	strZygosity = "Het"
							elif strArrayGenotypeInfo_25pcReads[0]=="1/1":
								strZygosity = "Hom"
				#Else if a copy number event
				elif strArrayGenotypeInfo[0]=="CNV":
					boolCN = "TRUE"

				if boolCN=="FALSE":
					#If the CHR column is in format 'chr[1-22XY]' or '[1-22XY]'
					if "chr" in strArray[0]:
						strArray[2] = strArray[0] + ":" + strArray[1] + ":" + strArray[3] + ":" + strArray[4] + ":" + strZygosity
					else:
						strArray[2] = "chr" + strArray[0] + ":" + strArray[1] + ":" + strArray[3] + ":" + strArray[4] + ":" + strZygosity

					#If we were unable to determine the zygosity, it's a low quality call; therefore, set the QUAL column as such
					if strZygosity=="NA":
						strArray[6] = "LowQual"
					else:
						strArray[6] = "."

			#Output the modified line
			fileOutput.write("\t".join(strArray) + "\n")

		#If the length of the variant base is greater than 1, it's an insertion in VCF format
		elif iVarBaseLength>1:
			#Remove the first base from the variant base call
			strNewVarBase = strArray[4][(1):]

			#Change the old reference base to '-'
			strArray[3] = "-"

			#Change the variant to 'ins' + reference bases that were deleted
			strArray[4] = strNewVarBase

			#Change the base-position by 1 base (adition)
			strArray[1] = str(int(strArray[1]) + 1)

			#If the user requested to change the ID field with CHR:POS:REF:VAR:ZYGOSITY
			if sys.argv[3].upper()=="TRUE":
				#Determine zygosity
				strArrayGenotypeInfo = [splits for splits in strArray[9].split(":") if splits is not ""]

				if strArrayGenotypeInfo[0]=="0/1" or strArrayGenotypeInfo[0]=="1/0":
					strZygosity = "Het"
				elif strArrayGenotypeInfo[0]=="1/1":
					strZygosity = "Hom"
				#If the genotype is missing, look further at the downsampled reads
                                #75% reads
                                elif strArrayGenotypeInfo[0]=="0/0":
                                        #Determine zygosity
                                        strArrayGenotypeInfo_75pcReads = [splits for splits in strArray[10].split(":") if splits is not ""]

                                        if strArrayGenotypeInfo_75pcReads[0]=="0/1" or strArrayGenotypeInfo_75pcReads[0]=="1/0":
                                                strZygosity = "Het"
                                        elif strArrayGenotypeInfo_75pcReads[0]=="1/1":
                                                strZygosity = "Hom"
                                        #50% reads
                                        elif strArrayGenotypeInfo[0]=="0/0":
                                                #Determine zygosity
                                                strArrayGenotypeInfo_50pcReads = [splits for splits in strArray[11].split(":") if splits is not ""]

                                                if strArrayGenotypeInfo_50pcReads[0]=="0/1" or strArrayGenotypeInfo_50pcReads[0]=="1/0":
                                                        strZygosity = "Het"
                                                elif strArrayGenotypeInfo_50pcReads[0]=="1/1":
                                                        strZygosity = "Hom"
                                                #25% reads
                                                elif strArrayGenotypeInfo[0]=="0/0":
                                                        #Determine zygosity
                                                        strArrayGenotypeInfo_25pcReads = [splits for splits in strArray[12].split(":") if splits is not ""]
                                                        if strArrayGenotypeInfo_25pcReads[0]=="0/1" or strArrayGenotypeInfo_25pcReads[0]=="1/0":
                                                                strZygosity = "Het"
                                                        elif strArrayGenotypeInfo_25pcReads[0]=="1/1":
                                                                strZygosity = "Hom"
				#Else if a copy number event
                                elif strArrayGenotypeInfo[0]=="CNV":
                                	boolCN = "TRUE"

				if boolCN=="FALSE":
					#If the CHR column is in format 'chr[1-22XY]' or '[1-22XY]'
					if "chr" in strArray[0]:
						strArray[2] = strArray[0] + ":" + strArray[1] + ":" + strArray[3] + ":" + strArray[4] + ":" + strZygosity
					else:
						strArray[2] = "chr" + strArray[0] + ":" + strArray[1] + ":" + strArray[3] + ":" + strArray[4] + ":" + strZygosity

					#If we were unable to determine the zygosity, it's a low quality call; therefore, set the QUAL column as such
                                        if strZygosity=="NA":
                                                strArray[6] = "LowQual"
                                        else:
                                                strArray[6] = "."

			#Output the modified line
			fileOutput.write("\t".join(strArray) + "\n")

		#Else if it's a single nucleotide variant
		else:
			#If the user requested to change the ID field with CHR:POS:REF:VAR:ZYGOSITY
			if sys.argv[3].upper()=="TRUE":
				#Determine zygosity
				strArrayGenotypeInfo = [splits for splits in strArray[9].split(":") if splits is not ""]

				if strArrayGenotypeInfo[0]=="0/1" or strArrayGenotypeInfo[0]=="1/0":
					strZygosity = "Het"
				elif strArrayGenotypeInfo[0]=="1/1":
					strZygosity = "Hom"
				#If the genotype is missing, look further at the downsampled reads
                                #75% reads
                                elif strArrayGenotypeInfo[0]=="0/0":
                                        #Determine zygosity
                                        strArrayGenotypeInfo_75pcReads = [splits for splits in strArray[10].split(":") if splits is not ""]

                                        if strArrayGenotypeInfo_75pcReads[0]=="0/1" or strArrayGenotypeInfo_75pcReads[0]=="1/0":
                                                strZygosity = "Het"
                                        elif strArrayGenotypeInfo_75pcReads[0]=="1/1":
                                                strZygosity = "Hom"
                                        #50% reads
                                        elif strArrayGenotypeInfo[0]=="0/0":
                                                #Determine zygosity
                                                strArrayGenotypeInfo_50pcReads = [splits for splits in strArray[11].split(":") if splits is not ""]

                                                if strArrayGenotypeInfo_50pcReads[0]=="0/1" or strArrayGenotypeInfo_50pcReads[0]=="1/0":
                                                        strZygosity = "Het"
                                                elif strArrayGenotypeInfo_50pcReads[0]=="1/1":
                                                        strZygosity = "Hom"
                                                #25% reads
                                                elif strArrayGenotypeInfo[0]=="0/0":
                                                        #Determine zygosity
                                                        strArrayGenotypeInfo_25pcReads = [splits for splits in strArray[12].split(":") if splits is not ""]
                                                        if strArrayGenotypeInfo_25pcReads[0]=="0/1" or strArrayGenotypeInfo_25pcReads[0]=="1/0":
                                                                strZygosity = "Het"
                                                        elif strArrayGenotypeInfo_25pcReads[0]=="1/1":
                                                                strZygosity = "Hom"
	
                                #Else if a copy number event
                                elif strArrayGenotypeInfo[0]=="CNV":
					boolCN = "TRUE"

				if boolCN=="FALSE":
					#If the CHR column is in format 'chr[1-22XY]' or '[1-22XY]'
					if "chr" in strArray[0]:
						strArray[2] = strArray[0] + ":" + strArray[1] + ":" + strArray[3] + ":" + strArray[4] + ":" + strZygosity
					else:
						strArray[2] = "chr" + strArray[0] + ":" + strArray[1] + ":" + strArray[3] + ":" + strArray[4] + ":" + strZygosity

					#If we were unable to determine the zygosity, it's a low quality call; therefore, set the QUAL column as such
                                        if strZygosity=="NA":
                                                strArray[6] = "LowQual"
                                        else:
                                                strArray[6] = "."

				fileOutput.write("\t".join(strArray) + "\n")
print ("Done.")

#Close the input and output file
fileInput.close()
fileOutput.close()

#Delete the temporary file
print("Removing temporary file...")
os.remove("temp.vcf.tmp")
print("End.")

