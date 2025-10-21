#!/bin/bash

# Default values
NSIMS=100
BURNIN=0

# Parse command-line arguments
while [[ $# -gt 0 ]]; do
  key="$1"
  case $key in
    --nsims)
      NSIMS="$2"
      shift 2
      ;;
    --burnin)
      BURNIN="$2"
      shift 2
      ;;
    *)
      echo "Unknown option: $1"
      echo "Usage: $0 [--nsims N] [--burnin B]"
      exit 1
      ;;
  esac
done

echo "Number of simulations: $NSIMS"
echo "Burnin proportion: $BURNIN"

# ---- Prepare templates ----
rm -r templates
mkdir -p templates

# ---- Simulate parameters ----
rm -r truth
mkdir -p truth
cd truth

# ---- Simulate parameters ----
java -jar ~/beastMCMC/GammaSpikeModelRhoSampling_feastFunction_jar/GammaSpikeModelRhoSampling.jar -overwrite -D nSims=$NSIMS ../truth.xml

# ---- Clean up .trees files ----
for file in ./*.trees; do
    awk '
    BEGIN {
        skip = 0
        intrees = 0
    }

    /^[[:space:]]*Begin taxa;/ { skip = 1; next }
    /^[[:space:]]*End;/ { if (skip == 1) { skip = 0; next } }
    skip == 1 { next }

    /^[[:space:]]*Begin trees;/ { intrees = 1; print; next }
    intrees == 1 && skip == 0 && /^[[:space:]]*Translate/ { skip = 2; next }
    skip == 2 { if ($0 ~ /;/) { skip = 0 }; next }

    { print }
    ' "$file" > "${file}.cleaned"

    mv "${file}.cleaned" "$file"
done

# ---- Add metadata ----
Rscript ../scripts/addTreeMetadata.R

cd ..

# ---- Populate JSONs with R ----
Rscript scripts/populate.R --burnin $BURNIN --nsims $NSIMS