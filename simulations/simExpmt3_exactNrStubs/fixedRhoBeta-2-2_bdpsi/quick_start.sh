sbatch --wrap="bash prepare.sh --nsims 100;bash submit_all_runs.sh --nsims 100 --seqlength 1000"
# bash prepare.sh --nsims 10