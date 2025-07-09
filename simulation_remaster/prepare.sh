#!/bin/bash

# Default values
NSIMS=50
BURNIN=0
SEQLEN=100
#NTAXA=20

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
    --seqlength)
      SEQLEN="$2"
      shift 2
      ;;
    # --ntaxa)
    #   NTAXA="$2"
    #   shift 2
    #   ;;
    *)
      echo "Unknown option: $1"
      # echo "Usage: $0 [--nsims N] [--burnin B] [--seqlength L] [--ntaxa T]"
      echo "Usage: $0 [--nsims N] [--burnin B] [--seqlength L]"
      exit 1
      ;;
  esac
done

echo "Number of simulations: $NSIMS"
echo "Burnin proportion: $BURNIN"
echo "Sequence length: $SEQLEN"
# echo "Number of taxa: $NTAXA"

# ---- Prepare templates ----
rm -r templates
mkdir -p templates

# ---- Simulate parameters ----
rm -r truth
mkdir -p truth
cd truth

# ---- Simulate parameters ----
beast -overwrite -D nSims=$NSIMS ../truth.xml

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

# ---- Simulate datasets ----
cd templates

for rep in rep*/;
do
	cd $rep
	beast -overwrite -D seqLength=$SEQLEN -df var.json ../../seqsim.xml
	echo "---------------"
	echo $rep
	echo "---------------"
	Rscript ../../scripts/putSequencesIntoJson.R
	cd ../
done

cd ../
