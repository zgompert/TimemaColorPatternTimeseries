#!/bin/bash 
#SBATCH --time=168:00:00
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
MAX_JOBS=6

# Input files and associated genome IDs
comps=(
  "cactusStripe_TcrGSH1_TcrE240087H1.txt"
  "cactusStripe_TcrGSH1_TcrE240087H2.txt"
  "cactusStripe_TcrGUSH2_TcrE240087H1.txt"
  "cactusStripe_TcrGUSH2_TcrE240087H2.txt"
  "cactusStripe_TcrE240087H1_TcrE240087H2.txt"
  "cactusStripe_TcrE240016H1_TcrE240087H1.txt"
  "cactusStripe_TcrE240016H2_TcrE240087H1.txt"
  "cactusStripe_TcrE240016H1_TcrE240087H2.txt"
  "cactusStripe_TcrE240016H2_TcrE240087H2.txt"
  "cactusStripe_TcrE240038H1_TcrE240087H1.txt"
  "cactusStripe_TcrE240038H1_TcrE240087H2.txt"
  "cactusStripe_TcrGSH2_TcrE240087H1.txt"	
  "cactusStripe_TcrGSH2_TcrE240087H2.txt"	
)

ids1=(
  "t_cris_h_gs1"
  "t_cris_h_gs1"
  "t_cris_h_gus2"
  "t_cris_h_gus2"
  "t_cris_e_240087h1"
  "t_cris_e_240016h1"
  "t_cris_e_240016h2"
  "t_cris_e_240016h1"
  "t_cris_e_240016h2"
  "t_cris_e_240038h1"
  "t_cris_e_240038h1"
  "t_cris_h_gs2"
  "t_cris_h_gs2"
)

ids2=(
  "t_cris_e_240087h1"
  "t_cris_e_240087h2"
  "t_cris_e_240087h1"
  "t_cris_e_240087h2"
  "t_cris_e_240087h2"
  "t_cris_e_240087h1"
  "t_cris_e_240087h1"
  "t_cris_e_240087h2"
  "t_cris_e_240087h2"
  "t_cris_e_240087h1"
  "t_cris_e_240087h2"
  "t_cris_e_240087h1"
  "t_cris_e_240087h2"
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


