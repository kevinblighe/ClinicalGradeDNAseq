
#!/bin/bash ;

######################################################### ;
#Pre-liminary step:  Safeguarding the script's operation# ;
######################################################### ;

echo "Checking command line parameters..." ;

Read1="${1}" ;
Read2="${2}" ;
Ref_FASTA="${3}" ;
RunNumber="${4}" ;
PatientID="${5}" ;
BEDfile="${6}" ;
TrimmingQualityReadEnds="${7}" ;
TrimmingReadLengthMin="${8}" ;
TrimmingAdaptor="${9}" ;
FilterReadDepthCutoff="${10}" ;
FilterMappingQualityCutOff="${11}" ;
CallingStringency="${12}" ;
Results_root="${13}" ;
User="${14}" ;

#Check that 14 parameters have been passed
if [ $# -ne 14 ] ;
then
        echo "Error - incorrect number of parameters. Consult the relevant standard operating procedure for correct usage." ;
	exit 1 ;
fi

#Check that the file $Read1 exists
test -e "${Read1}" ;
if [ $? -ne 0 ] ;
then
        echo "Error - ${Read1} does not exist. Please check the complete file path and re-run." ;
	exit 1 ;
fi

#Check that the file $Read2 exists
test -e "${Read2}" ;
if [ $? -ne 0 ] ;
then
        echo "Error - ${Read2} does not exist. Please check the complete file path and re-run." ;
        exit 1 ;
fi

#Check that the file $Ref_FASTA exists
test -e "${Ref_FASTA}" ;
if [ $? -ne 0 ] ;
then
        echo "Error - ${Ref_FASTA} does not exist. Please check the complete file path and re-run." ;
        exit 1 ;
fi

#Check that the file $BEDfile exists
test -e "${BEDfile}" ;
if [ $? -ne 0 ] ;
then
        echo "Error - ${BEDfile} does not exist. Please check the complete file path and re-run." ;
        exit 1 ;
fi

#Check that various parameters are integers
if [[ ! "${TrimmingQualityReadEnds}" =~ ^[0-9]*$ ]] ;
then
	echo "Error - "${TrimmingQualityReadEnds}" must be an integer." ;
	exit 1 ;
fi
if [[ ! $TrimmingReadLengthMin =~ ^[0-9]*$ ]] ;
then
        echo "Error - "${TrimmingReadLengthMin}" must be an integer." ;
        exit 1 ;
fi
if [[ ! $FilterReadDepthCutoff =~ ^[0-9]*$ ]] ;
then
        echo "Error - "${FilterReadDepthCutoff}" must be an integer." ;
        exit 1 ;
fi
if [[ ! $FilterMappingQualityCutOff =~ ^[0-9]*$ ]] ;
then
        echo "Error - "${FilterMappingQualityCutOff}" must be an integer." ;
        exit 1 ;
fi

#Check that the specified trimming adaptor is correct
if [ $TrimmingAdaptor != "illumina" ] ;
then
	if [ $TrimmingAdaptor != "nextera" ] ;
	then
		if [ $TrimmingAdaptor != "small_rna" ] ;
		then
		        echo "Error - "${TrimmingAdaptor}" must be one of illumina|nextera|small_rna." ;
		        exit 1 ;
		fi
	fi
fi

#Check that the specified trimming adaptor is correct
if [ $CallingStringency != "relaxed" ] ;
then
        if [ $CallingStringency != "normal" ] ;
        then
		echo "Error - "${CallingStringency}" must be one of normal|relaxed." ;
		exit 1 ;
        fi
fi

echo "Done." ;

echo "Running master analysis script on `date`..." ;

mkdir -p "${Results_root}" ;
mkdir -p -m 777 "${Results_root}"/"${RunNumber}" ;
mkdir -p -m 777 "${Results_root}"/"${RunNumber}"/"${PatientID}" ;
mkdir -p -m 777 "${Results_root}"/"${RunNumber}"/"${PatientID}"/tmp ;

echo "Log file is "${Results_root}"/"${RunNumber}"/"${PatientID}"/Master.log" ;

/home/ubuntu/pipeline/AnalysisMasterVersion1.sh "${Read1}" "${Read2}" "${Ref_FASTA}" "${RunNumber}" "${PatientID}" "${BEDfile}" "${TrimmingQualityReadEnds}" "${TrimmingReadLengthMin}" "${TrimmingAdaptor}" "${FilterReadDepthCutoff}" "${FilterMappingQualityCutOff}" "${CallingStringency}" "${Results_root}" "${User}" &> "${Results_root}"/"${RunNumber}"/"${PatientID}"/Master.log ;

echo `unix2dos "${Results_root}"/"${RunNumber}"/"${PatientID}"/Master.log` ;

echo "Transferring results via sFTP to local storage..." ;

checkCredentials=1 ;
remoteDir="/remote/SAMBA/share/" ;

#The first sshpass command will return 0 if successful, >0 if not
#The loop will continue until correct credentials are supplied and access gained
while [[ checkCredentials -eq 1 ]] ;
do
	echo "Enter authentication username:" ;
	read username ;
	echo "Enter password:" ;
	read -s SSHPASS ;
	eval "export SSHPASS='""$SSHPASS""'" ;

	sshpass -e sftp $username@XXX.XXX.XXX.XXX << !
!

	#Check the exit code of the sshpass command
	if [ $? -eq 0 ]
	then
		checkCredentials=0 ;
	fi
done

#Make the actual transfer
sshpass -e sftp $username@XXX.XXX.XXX.XXX << !
        cd "${remoteDir}" ;

        mkdir "${remoteDir}"/"${RunNumber}"/ ;
	mkdir "${remoteDir}"/"${RunNumber}"/"${PatientID}"/ ;

        put -r "${Results_root}"/"${RunNumber}"/"${PatientID}"/* "${remoteDir}"/"${RunNumber}"/"${PatientID}"/ ;
!

echo "Finished on `date`." ;

exit 0 ;
