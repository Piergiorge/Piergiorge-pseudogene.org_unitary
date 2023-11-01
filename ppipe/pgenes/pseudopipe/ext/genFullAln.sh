#!/bin/sh
if [ ! -d "$1" ]; then echo 'Please specify pgenes directory'; exit 0; fi
if [ -z "$2" ]; then echo 'Please specify output file'; exit 0; fi
if [ -f "$2" ]; then rm -f $2; fi

for i in `cd $1; pwd`/pgenes/*/pgenes/*.all.fa
do 
	gzip -c $i >> $2;
done

echo Finished generating pgene full alignment
