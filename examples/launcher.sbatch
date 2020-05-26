#!/bin/bash
#SBATCH -J dge
#SBATCH -o dge_%j.out
#SBATCH -e dge_%j.err
#SBATCH -t 24:00:00
#SBATCH -p all
#SBATCH -N 2 -n 7
#SBATCH -A SD2E-Community

ml tacc-singularity

## you'll need the container the first time you run this,
## but can comment out this block after the first run
singularity pull docker://jurrutia/omics_tools:0.1.3
##

ml launcher
export LAUNCHER_PLUGIN_DIR=$LAUNCHER_DIR/plugins
export LAUNCHER_RMI=SLURM
export LAUNCHER_JOB_FILE=launcher_file.sh
export LAUNCHER_WORKDIR=`pwd`

$LAUNCHER_DIR/paramrun
