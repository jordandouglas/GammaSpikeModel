#!/bin/bash
#SBATCH --job-name=rep1
#SBATCH --output=rep1-%j.out
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1G
#SBATCH --time=05:00:00

cd "$SLURM_SUBMIT_DIR/templates/rep1"
java -jar ~/beastMCMC/GammaSpikeModelRhoSampling_origin_jar/GammaSpikeModelRhoSampling.jar -overwrite -D seqLength=1000 -df var.json -DFout seqsim.out.xml ../../seqsim.xml
Rscript ../../scripts/putSequencesIntoJson.R
java -jar ~/beastMCMC/GammaSpikeModelRhoSampling_origin_jar/GammaSpikeModelRhoSampling.jar -overwrite -df var.seq.json -DFout run.out.xml ../../run.xml
