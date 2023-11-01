#/bin/sh
if [ $# -lt 5 ]
then
        echo "Usage: ppipe [output dir] [masked dna dir] [input dna dir] [input pep dir] [exon mask dir] [PBS?]"
        exit 1
fi

. `dirname $0`/../bin/env.sh

outDir=`if [ ! -d $1 ]; then mkdir -p $1; fi;cd $1;pwd`
rmkDir=`if [ -d $2 ]; then echo \`cd $2;pwd\`; else echo $2; fi`
pepDir=`if [ -d $4 ]; then echo \`cd $4;pwd\`; else echo $4; fi`
dnaTmp=`if [ -d $3 ]; then echo \`cd $3;pwd\`/%s.fa; else echo $3; fi`
emkTmp=`if [ -d $5 ]; then echo \`cd $5;pwd\`/%s_exLocs; else echo $5; fi`

inputDNA=$outDir/dna/dna_all.fa
inputPEP=$outDir/pep/pep_all.fa

jobsExec=$sqDummy; if [ "$6" = "1" ]; then jobsExec=$sqDedicated; fi


echo Making directories
cd $outDir
if [ ! -d dna ]; then mkdir dna; fi
if [ ! -d pep ]; then mkdir pep; fi
if [ ! -d blast ]
then 
	mkdir blast
	cd blast
	mkdir stamps
	mkdir output
	mkdir status
	mkdir processed
	cd ..
fi
if [ ! -d pgenes ]
then
        mkdir pgenes
        cd pgenes
        mkdir minus plus
	cd minus
	mkdir log stamp
	cd ../plus
	mkdir log stamp
        cd ../../
fi

echo Copying sequences
if [ ! -f $inputDNA ]
then
	if [ -f $rmkDir ]
	then 
		inputDNA=$rmkDir
	else
		cat $rmkDir/*.fa > $inputDNA
	fi
fi
if [ ! -f $inputPEP ]
then
	if [ -f $pepDir ]
        then
		inputPEP=$pepDir
        else
		cat $pepDir/*.fa > $inputPEP
	fi
fi

echo Fomatting the DNAs
cd dna
if [ ! -f formatdb.log ]
then
	if [ -f $rmkDir -a -f $rmkDir.nin ]
	then
		for i in $rmkDir.*
		do
			ext=`echo $i | sed 's/.*\.\(.*\)/\1/'`
			if [ ! -e $inputDNA.$ext ]
			then
				ln -s $i $inputDNA.$ext
			fi
		done
		touch formatdb.log
	else
		$formatDB -i $inputDNA -o T -p F
	fi
fi
cd ..

echo Preparing the blast jobs
cd blast

if [ ! -f jobs ]
then
	pepNum=`grep '>' $inputPEP | wc | sed 's/\s*\([0-9]*\).*/\1/'`
	$pythonExec $fastaSplitter $inputPEP $((($pepNum+359)/360)) 'split%04d'
	for f in `ls split*`
	do
		echo "\"( cd $(pwd); touch stamps/${f}.Start ; ( $blastExec -p tblastn -m 8 -z 3.1e9 -e .1 -d $inputDNA -i $f -o output/${f}.Out ; touch stamps/${f}.Stamp ) >status/${f}.Status 2>&1 )\"" 
	done > jobs

	$pythonExec $jobsExec jobs
	
	echo Finished blast
else
	echo Skipping blast
fi

echo Processing blast output
cd processed
if [ ! -f ./processed.stamp ]
then
        echo "\"(cd $(pwd); $pythonExec $blastHandler $inputPEP  'split\d{4}.Out\Z' ../output; touch processed.stamp)\"" > jobs

	$pythonExec $jobsExec jobs
	
	echo Finished processing blast output
else
        echo Skipping the processing of blast output
fi

echo Running Pseudopipe on both strands
cd ../../pgenes

for t in 'M' 'P'
do
	echo 'Working on '$t' strand'

	if [ $t = 'M' ]
	then
		cd minus
	else
		cd plus
	fi

	echo "export BlastoutSortedTemplate=${outDir}/blast/processed/%s_${t}_blastHits.sorted;export ChromosomeFastaTemplate=${dnaTmp};export ExonMaskTemplate=${emkTmp};export ExonMaskFields='2 3';export FastaProgram=${fastaExec};export ProteinQueryFile=${inputPEP}" > setenvPipelineVars

	ms=$(cd ../../blast/processed ; for f in *_${t}_*sorted; do echo ${f/_${t}_blastHits.sorted/}; done)

	for c in $ms
	do
        	echo "\"(cd $(pwd); source setenvPipelineVars; touch stamp/$c.Start ; $pythonExec $pseudopipe $c > log/$c.log 2>&1; touch stamp/$c.Stop)\""
	done > jobs

	$pythonExec $jobsExec jobs

	echo Finished Pseudopipe on strand $t

	cd ..
done

echo Generating final results
outFilePrefix=$outDir/pgenes/`basename $outDir`_pgenes
$genPgeneResult $outDir $outFilePrefix.txt
$genFullAln     $outDir $outFilePrefix.align.gz

echo Finished running Pseudopipe
