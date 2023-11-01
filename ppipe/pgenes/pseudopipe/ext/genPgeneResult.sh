#/bin.sh
if [ $# -lt 2 ]
then
	echo usage: ppResult [output dir] [output file]
	exit 0
fi

echo -e "#chr\\tstart\\tend\\tstrand\\tquery\\tfrac\\tins\\tdel\\tshift\\tstop\\texpect\\tident\\tpolya\\ttype" > $2

(for t in 'minus' 'plus'
do
	for i in `find $1/pgenes/${t}/pgenes -name '*.all.gff'`
	do
		sed '1d' $i | cut -f 1,5-8,14-21,29 | grep -v 'GENE' | sed 's/\(.*PSSD\)[0-9]/\1/' | sed 's/\tFP$/\tFRAG/'
	done
done) >> $2
