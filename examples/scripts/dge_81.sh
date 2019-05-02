#!/usr/bin/env bash
source /broad/software/scripts/useuse
use R-3.5

WD="/btl/foundry/users/alex/20190228_novel_chassis/run_DGE/"
RESULTS=$JOB_ID
RUNDIR="$WD$RESULTS"
mkdir $RUNDIR
cd $RUNDIR

Rscript $WD/scripts/dge_81.r > $WD/$RESULTS/dge_81_log.txt
