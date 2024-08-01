# !/bin/bash

# Check if the correct number of arguments is provided
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <directory> <file>"
    exit 1
fi

DIR=$1

if [[ "${DIR}" != */ ]]; then
    DIR="${DIR}/"
fi

FILE=$2

NAME=$(echo $FILE | cut -d'.' -f1)

# Run fasterq-dump for the specified files
fasterq-dump ${DIR}${FILE} -O ${DIR}${NAME} -t /scratch/jjacobs/tmp/

# Check if fasterq-dump was successful
if [ $? -ne 0 ]; then
    echo "fasterq-dump failed"
    exit 1
fi

# Run make_fasta.py with the generated FASTQ files
python make_fasta.py ${DIR}${NAME} 

# Check if make_fasta.py was successful
if [ $? -ne 0 ]; then
    echo "make_fasta.py failed"
    exit 1
fi