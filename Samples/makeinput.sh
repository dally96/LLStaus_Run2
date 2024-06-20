#! /bin/bash

#For use run 'sh makeinput.sh /TTto.../Run3Summer22EE.../MINIAODSIM'

#Use this to split the DAS entry for the sample name
IFS="/"

#Stores the parts of the entry into an array
read -ra arr <<< "$1"

#file_name used for the name of the output file to use for the tagger
#sample is to rejoin the parts of the argument together
file_name=""
sample=""

#Appends the file_name and sample with the parts of the sample in the array
#If statements are to make sure there are no underscores at the beginning or end of the file name
for val in "${arr[@]}"
do 
  if [[ $val = "${arr[0]}" ]]
  then
    continue
  fi
  if [[ $val = "${arr[-1]}" ]]
  then
    file_name+=$val
  else
    file_name+=$val
    file_name+="_"
  fi
  sample+="/"
  sample+=$val 
done

#Reset the IFS split character so sample doesn't get split when printing it into the file
IFS=""

#Make the output file a txt file
file_name+=".txt"

#Make the output file
touch ${file_name}

echo isMC=True disTauTagOutputOpt=1 era=2022 >> ${file_name}
echo "" >> ${file_name}
echo ${arr[1]} $sample >> ${file_name}

printf "File for sample $sample has been made\n"



