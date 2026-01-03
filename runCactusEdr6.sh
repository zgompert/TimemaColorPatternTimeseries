#!/bin/bash 
#SBATCH --time=96:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=96
#SBATCH --mem=724000
#SBATCH --account=gompert
#SBATCH --qos=gompert-grn
#SBATCH --partition=gompert-grn
#SBATCH --job-name=cactus
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=zach.gompert@usu.edu

module load cactus

cd /scratch/general/nfs1/u6000989/cactus

# Max number of concurrent jobs
MAX_JOBS=4

# Input files and associated genome IDs
comps=(
"cactusStripe_TcrGSH1_TcrE240028H1.txt"
"cactusStripe_TcrGSH1_TcrE240029H1.txt"
"cactusStripe_TcrGSH1_TcrE240030H1.txt"
"cactusStripe_TcrE240028H1_TcrE240028H2.txt"
"cactusStripe_TcrE240029H1_TcrE240029H2.txt"
"cactusStripe_TcrE240030H1_TcrE240030H2.txt"
"cactusStripe_TcrE240016H1_TcrE240028H1.txt"
"cactusStripe_TcrE240016H1_TcrE240029H1.txt"
"cactusStripe_TcrE240016H1_TcrE240030H1.txt"
)

ids1=(
  "t_cris_h_gs1"
  "t_cris_h_gs1"
  "t_cris_h_gs1"
  "t_cris_e_240028h1"
  "t_cris_e_240029h1"
  "t_cris_e_240030h1"
  "t_cris_e_240016h1"
  "t_cris_e_240016h1"
  "t_cris_e_240016h1"
)

ids2=(
  "t_cris_e_240028h1"
  "t_cris_e_240029h1"
  "t_cris_e_240030h1"
  "t_cris_e_240028h2"
  "t_cris_e_240029h2"
  "t_cris_e_240030h2"
  "t_cris_e_240028h1"
  "t_cris_e_240029h1"
  "t_cris_e_240030h1"
)

# Loop over jobs
for i in "${!comps[@]}"; do
  file="${comps[$i]}"
  id1="${ids1[$i]}"
  id2="${ids2[$i]}"

  if [[ "$file" =~ cactusStripe_([a-zA-Z0-9]+)_([a-zA-Z0-9]+)\.txt ]]; then
    base="${BASH_REMATCH[1]}_${BASH_REMATCH[2]}"
  else
    echo "Failed to parse: $file"
    continue
  fi

  (
    cactus "jobStore_$base" \
      "/uufs/chpc.utah.edu/common/home/gompert-group4/data/timema/hic_genomes/comp_aligns/$file" \
      "cactus${base}.hal" --maxCores 24

    ~/source/hal/bin/halSynteny --queryGenome "$id2" --targetGenome "$id1" \
      "cactus${base}.hal" "out_synteny_${base}.psl"
  ) &

  # Limit the number of background jobs
  while (( $(jobs -rp | wc -l) >= MAX_JOBS )); do
    wait -n
  done
done

# Wait for all remaining background jobs to finish
wait


