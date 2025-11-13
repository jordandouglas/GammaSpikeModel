#!/bin/bash

# Initialize variable
NSIMS=""
SEQLEN=1000 # Default value

# Argument parsing
while [[ "$#" -gt 0 ]]; do
  case $1 in
    --nsims)
      NSIMS="$2"
      shift
      ;;
    --seqlength)
      SEQLEN="$2"
      shift 2
      ;;
    -h|--help)
      echo "Usage: $0 --nsims N [--seqlength L]"
      exit 0
      ;;
    *)
      echo "Unknown parameter passed: $1"
      exit 1
      ;;
  esac
  shift
done

# Check that --nsims was provided
if [[ -z "$NSIMS" ]]; then
  echo "Error: --nsims is required."
  echo "Usage: $0 --nsims N"
  exit 1
fi

# Validate it's a positive integer
if ! [[ "$NSIMS" =~ ^[1-9][0-9]*$ ]]; then
  echo "Error: --nsims must be a positive integer."
  exit 1
fi

echo "Number of simulations: $NSIMS"
echo "Sequence length: $SEQLEN"

# Move into templates directory
cd templates || { echo "Directory 'templates' not found."; exit 1; }

# Loop through replicates
for ((i=1; i<=NSIMS; i++)); do
  rep_dir="rep$i"

  if [ ! -d "$rep_dir" ]; then
    echo "Directory $rep_dir not found, skipping..."
    continue
  fi

  cd "$rep_dir"

  # Create SLURM job script
  job_script="job_rep${i}.sh"

  cat > "$job_script" <<EOF
#!/bin/bash
#SBATCH --job-name=rep${i}
#SBATCH --output=rep${i}-%j.out
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1G
#SBATCH --time=12:00:00

cd "\$SLURM_SUBMIT_DIR/templates/$rep_dir"
java -jar ~/beastMCMC/GammaSpikeModelRhoSampling_origin_jar/GammaSpikeModelRhoSampling.jar -overwrite -D seqLength=$SEQLEN -df var.json -DFout seqsim.out.xml ../../seqsim.xml
Rscript ../../scripts/putSequencesIntoJson.R
java -jar ~/beastMCMC/GammaSpikeModelRhoSampling_origin_jar/GammaSpikeModelRhoSampling.jar -overwrite -df var.seq.json -DFout run.out.xml ../../run.xml
EOF

  chmod +x "$job_script"
  echo "Submitting $job_script in templates/$rep_dir ..."
  sbatch "$job_script"

  cd ..
done