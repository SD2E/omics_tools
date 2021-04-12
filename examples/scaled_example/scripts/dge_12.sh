#!/usr/bin/env bash
source /broad/software/scripts/useuse
use R-3.5
WD="/data/."
RESULTS=results
RUNDIR="$WD/$RESULTS"
mkdir $RUNDIR
cd $RUNDIR

Rscript $WD/scripts/dge_12.r > $WD/$RESULTS/dge_12_log.txt
