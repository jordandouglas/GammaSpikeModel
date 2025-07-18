sbatch --wrap="bash prepare.sh --nsims 200 --seqlength 500;bash submit_all_runs.sh --nsims 200"
