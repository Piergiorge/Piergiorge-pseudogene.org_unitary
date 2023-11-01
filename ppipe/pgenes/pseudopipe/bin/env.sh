#!/bin/sha

if [ ! -z "$PSEUDOPIPE_ENV" ]; then source $PSEUDOPIPE_ENV; return; fi

# Pseudopipe configuration
export PSEUDOPIPE_HOME=`cd \`dirname $0\`/../; pwd`
export pseudopipe=$PSEUDOPIPE_HOME/core/runScripts.py
export genPgeneResult=$PSEUDOPIPE_HOME/ext/genPgeneResult.sh
export genFullAln=$PSEUDOPIPE_HOME/ext/genFullAln.sh
export fastaSplitter=$PSEUDOPIPE_HOME/ext/splitFasta.py
export sqDummy=$PSEUDOPIPE_HOME/ext/sqDummy.py
export blastHandler=$PSEUDOPIPE_HOME/core/processBlastOutput.py
export extractExLoc=$PSEUDOPIPE_HOME/core/extractKPExonLocations.py

# Python configuration
export pythonExec=/home/bp272/bin/Python-2.6.6/python

# Alignment tools configuration
export formatDB=/home/bp272/bin/blast-2.2.13/bin/formatdb
export blastExec=/home/bp272/bin/blast-2.2.13/bin/blastall
export fastaExec=/home/bp272/bin/fasta-35.1.5/tfasty35
