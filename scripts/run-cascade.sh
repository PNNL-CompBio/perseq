#!/bin/sh
#SBATCH --account=emsls50064
#SBATCH --partition=long
#SBATCH --time=10000
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --job-name="perseq"
#SBATCH --output="logs/%A.out"
#SBATCH --error="logs/%A.err"

CONFIG="config/config.yml"
CLUSTERCONFIG="config/cluster.yml"
LOGS=logs

if [ ! -d "$LOGS" ]; then
    mkdir -p "$LOGS"
fi

snakemake \
    --configfile config/config.yml \
    --jobs 5000 \
    --use-conda \
    --printshellcmds \
    --cluster "'sbatch --account {cluster.account} --partition {cluster.partition} --exclusive --nodes 1 --time {cluster.time} --job-name {cluster.name} --output {cluster.output} --error {cluster.error}'" \
    --cluster-config cluster.yml
