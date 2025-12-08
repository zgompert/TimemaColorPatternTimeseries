# TimemaColorPatternTimeseries
Combined analyses of color and pattern in Timema cristinae, genetic mapping, comparative alignments and time series

# Timema cristinae comparative alignments

I am conducting comparative alignments of all of our phased (haplotype resolved) *Timema cristinae* genomes. This is an ongoing enterprise. I am starting with a bunch of pairwise alignments, but also am trying various approaches to align many genomes together.

## Genomes

We have four phased genomes that were part of (Gompert et al. 2025)[https://www.science.org/doi/full/10.1126/science.adp3745]. We also have newer phased genomes from Edinburgh genomics, which are in `/uufs/chpc.utah.edu/common/home/gompert-group4/data/timema/hic_genomes/edinburgh`. The cen* genomes are from Dovetail, the rest are Edinburh. Here is a summary of where things stand:

| ID | Location | Phenotype | Cactus aligns | SibeliaZ aligns |
|---------|-----|---------|:-:|:-:-|
| cen4119 | VP  | Stripe  | Y | Y |
| cen4280 | VP  | Green   | Y | Y |
| cen4120 | R12 | Green   | Y | Y |
| cen4122 | R23 | Stripe  | Y | Y |
| 24_0016 | VP  | Green   | Y | Y |
| 24_0028 | VP  | Green   | N | N |
| 24_0029 | VP  | Green   | N | N |
| 24_0030 | VP  | Green   | N | N |
| 24_0038 | VP  | Melanic | Y | Y |
| 24_0039 | VP  | Melanic | Y | Y |
| 24_0072 | R12 | Stripe  | Y | Y |
| 24_0073 | R12 | Green   | Y | Y |
| 24_0087 | VP  | Stripe  | Y | Y |
| 24_0175 | FH  | Stripe  | Y | Y |
| 24_0176 | FH  | Stripe  | Y | Y |

My first step with each genome is to split the fasta into files per haplotype and then to run repeat masking. This is done with `repeatmasker` (version 4.0.7); here is an example with six phased genomes:

```bash
#!/bin/sh 
#SBATCH --time=96:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --mem=384000
#SBATCH --account=gompert
#SBATCH --qos=gompert-grn
#SBATCH --partition=gompert-grn
#SBATCH --job-name=repeat
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=zach.gompert@usu.edu

module load repeatmasker

#version 4.0.7
cd /uufs/chpc.utah.edu/common/home/gompert-group4/data/timema/hic_genomes/repeat_mask

## run repeat masker on each genome sequencei
## uses library from the 2020 Science paper developed by Victor

MAX_JOBS=4

genomes=(
	"/uufs/chpc.utah.edu/common/home/gompert-group4/data/timema/hic_genomes/edingburgh/24_0028/Hap1Chr.fasta"
	"/uufs/chpc.utah.edu/common/home/gompert-group4/data/timema/hic_genomes/edingburgh/24_0028/Hap2Chr.fasta"
	"/uufs/chpc.utah.edu/common/home/gompert-group4/data/timema/hic_genomes/edingburgh/24_0029/Hap1Chr.fasta"
	"/uufs/chpc.utah.edu/common/home/gompert-group4/data/timema/hic_genomes/edingburgh/24_0029/Hap2Chr.fasta"
	"/uufs/chpc.utah.edu/common/home/gompert-group4/data/timema/hic_genomes/edingburgh/24_0030/Hap1Chr.fasta"
	"/uufs/chpc.utah.edu/common/home/gompert-group4/data/timema/hic_genomes/edingburgh/24_0030/Hap2Chr.fasta"
)

for file in "${genomes[@]}"; do

	RepeatMasker -s -e ncbi -xsmall -pa 24 -lib RepeatLibMergeCentroidsRM.lib $file &
	# Limit the number of background jobs
	while (( $(jobs -rp | wc -l) >= MAX_JOBS )); do
		wait -n
  	done
done

# Wait for all remaining background jobs to finish
wait
```
  
## Comparative alignments with Cactus



I ran repeat masking on them and did a big series of preliminary alignments, see [SynPlotTcrEd.R](SynPlotTcrEd.R)

- I then extracted chromosome 8 from all of the phased, *T. cristinae* genomes, see [extractCh8.pl](extractCh8.pl) and `/uufs/chpc.utah.edu/common/home/gompert-group4/data/timema/hic_genomes/comp_aligns/chr8haplotypes`

- I used SibeliaZ to do an initial alignment of this (it is fast but not as good at deeper divergence)

```bash
#!/bin/bash 
#SBATCH --time=96:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=96
#SBATCH --mem=724000
#SBATCH --account=gompert
#SBATCH --qos=gompert-grn
#SBATCH --partition=gompert-grn
#SBATCH --job-name=sibelia
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=zach.gompert@usu.edu

## sibelia comparative alignment for ch8

cd /uufs/chpc.utah.edu/common/home/gompert-group4/data/timema/hic_genomes/comp_aligns/chr8haplotypes

## see https://github.com/medvedevgroup/SibeliaZ
## for options
## -a = 24 genomes * 4 * 2
~/source/SibeliaZ/build/bin/bin/sibeliaz -k 25 -a 192 -b 200 -m 50 -t 180 ch8_t_cris_e_240016h1.fasta ch8_t_cris_e_240072h1.fasta ch8_t_cris_e_240175h1.fasta ch8_t_cris_h_gus1.fasta ch8_t_cris_e_240016h2.fasta ch8_t_cris_e_240072h2.fasta ch8_t_cris_e_240175h2.fasta ch8_t_cris_h_gus2.fasta ch8_t_cris_e_240038h1.fasta ch8_t_cris_e_240073h1.fasta ch8_t_cris_e_240176h1.fasta ch8_t_cris_r_gs1.fasta ch8_t_cris_e_240038h2.fasta ch8_t_cris_e_240073h2.fasta ch8_t_cris_e_240176h2.fasta ch8_t_cris_r_gs2.fasta ch8_t_cris_e_240039h1.fasta ch8_t_cris_e_240087h1.fasta ch8_t_cris_h_gs1.fasta ch8_t_cris_r_gus1.fasta ch8_t_cris_e_240039h2.fasta ch8_t_cris_e_240087h2.fasta ch8_t_cris_h_gs2.fasta ch8_t_cris_r_gus2.fasta
```
- The results are in `/uufs/chpc.utah.edu/common/home/gompert-group4/data/timema/hic_genomes/comp_aligns/chr8haplotypes/sibeliaz_out`

- I converted the aligment file to mfa format using [maftools](https://github.com/dentearl/mafTools/tree/master?tab=readme-ov-file)

```bash
 ~/source/mafTools/bin/mafToFastaStitcher --maf sibeliaz_out/alignment.maf --seqs ch8_combined.fasta --breakpointPenalty 5 --interstitialSequence 20  --outMfa output.mfa
```

- Next I used [trimAl](https://trimal.readthedocs.io/en/latest/usage.html) to remove any positions from the alignment with gaps.

```bash
## auto, not used
## trimal -in output.mfa -out clean.mfa -gappyout
## remove positions with any gaps
trimal -in output.mfa -out clean.mfa -gt 1.0
```

- I used castertree topology with caster

```bash
~/../gompert-group5/projects/LycAdmix/Caster/ASTER-Linux/bin/caster-site -i clean.mfa -o cout_tree --thread 24
```

- Finally, I added ML branch lengths to the tree in R, see [AddBlen.R](AddBlen.R).

```R
## adds ml branch lengths to the tree
library(phangorn)
library(ape)

packageVersion("phangorn")
#[1] ‘2.12.1’
packageVersion("ape")
#[1] ‘5.8’

## read tree, just topology, from caster
tre<-read.tree("topo_cout_tree")

## read the alignment and convert to phyDat object
aln<-read.dna("clean.mfa",format="fasta")
alnPh<-phyDat(aln,type="DNA")

## add arbitrary branch lengths as a starting point
treb<-compute.brlen(tre,1)

## compute the initial likelihood and then optimize branch lenghts under JC
## given the topoloty
fit<-pml(tree=treb,data=alnPh)
fit_optimized <- optim.pml(fit, model = "JC", optEdge = TRUE, optGamma = FALSE, optInv = FALSE)

## I might want to increase these slightly so they are overestimates
fit_optimized$tree$edge.length<-fit_optimized$tree$edge.length * 1.5

## write the tree
write.tree(phy=fit_optimized$tree,file="topoB_cout_tree")

## root on refugio and write again
refugio<-fit_optimized$tree$tip.label[9:16]
trr<-root(fit_optimized$tree,outgroup=refugio)
write.tree(phy=trr,file="topoRoot_cout_tree")
```

- Lastly, I used this to make the Cactus infile, I had to manually add the ancestors and fix the Refugio root

```
(((((((t_cris_e_240087h2:0.4426538755,t_cris_h_gs2:0.4582198529):0.07364615261,t_cris_e_240087h1:0.4635711029)AVPS:0.1559547772,(((t_cris_e_240016h2:0.3290234023,t_cris_h_gs1:0.2991242597):0.04056500028,t_cris_e_240039h1:0.4349549599):0.03885041141,(t_cris_e_240038h1:0.292011903,t_cris_e_240038h2:0.3390838267):0.06219885817)AVPM:0.03646449823):0.06323898609,t_cris_e_240039h2:0.4271548927):0.02428566702,(((t_cris_e_240176h1:0.304560471,t_cris_e_240175h1:0.3542393032):0.06988151601,t_cris_e_240175h2:0.2887467482):0.09896247756,t_cris_e_240176h2:0.3692333657)AFHS:0.05223943831):0.032424558,((t_cris_h_gus1:0.2903114095,t_cris_h_gus2:0.4758871493):0.1156325651,t_cris_e_240016h1:0.4679979238)AVPG:0.05056619986)AH154:0.05505660008,(((((t_cris_r_gus2:0.3435745311,t_cris_r_gus1:0.3564212078):0.1041618822,t_cris_e_240073h1:0.3693590572):0.04537837434,t_cris_e_240073h2:0.415594649)ARG:0.07449178421,((t_cris_r_gs2:0.346746946,t_cris_e_240072h1:0.3992755297):0.06593359072,t_cris_r_gs1:0.2781758046)ARS:0.1200039446):0.1166655065,t_cris_e_240072h2)AREF:0.4529215631);

t_cris_e_240016h1 /uufs/chpc.utah.edu/common/home/gompert-group4/data/timema/hic_genomes/edingburgh/24_0016/Hap1Chr.fasta.masked
t_cris_e_240016h2 /uufs/chpc.utah.edu/common/home/gompert-group4/data/timema/hic_genomes/edingburgh/24_0016/Hap2Chr.fasta.masked
t_cris_e_240038h1 /uufs/chpc.utah.edu/common/home/gompert-group4/data/timema/hic_genomes/edingburgh/24_0038/Hap1Chr.fasta.masked
t_cris_e_240038h2 /uufs/chpc.utah.edu/common/home/gompert-group4/data/timema/hic_genomes/edingburgh/24_0038/Hap2Chr.fasta.masked
t_cris_e_240039h1 /uufs/chpc.utah.edu/common/home/gompert-group4/data/timema/hic_genomes/edingburgh/24_0039/Hap1Chr.fasta.masked
t_cris_e_240039h2 /uufs/chpc.utah.edu/common/home/gompert-group4/data/timema/hic_genomes/edingburgh/24_0039/Hap2Chr.fasta.masked
t_cris_e_240072h1 /uufs/chpc.utah.edu/common/home/gompert-group4/data/timema/hic_genomes/edingburgh/24_0072/Hap1Chr.fasta.masked
t_cris_e_240072h2 /uufs/chpc.utah.edu/common/home/gompert-group4/data/timema/hic_genomes/edingburgh/24_0072/Hap2Chr.fasta.masked
t_cris_e_240073h1 /uufs/chpc.utah.edu/common/home/gompert-group4/data/timema/hic_genomes/edingburgh/24_0073/Hap1Chr.fasta.masked
t_cris_e_240073h2 /uufs/chpc.utah.edu/common/home/gompert-group4/data/timema/hic_genomes/edingburgh/24_0073/Hap2Chr.fasta.masked
t_cris_e_240087h1 /uufs/chpc.utah.edu/common/home/gompert-group4/data/timema/hic_genomes/edingburgh/24_0087/Hap1Chr.fasta.masked
t_cris_e_240087h2 /uufs/chpc.utah.edu/common/home/gompert-group4/data/timema/hic_genomes/edingburgh/24_0087/Hap2Chr.fasta.masked
t_cris_e_240175h1 /uufs/chpc.utah.edu/common/home/gompert-group4/data/timema/hic_genomes/edingburgh/24_0175/Hap1Chr.fasta.masked
t_cris_e_240175h2 /uufs/chpc.utah.edu/common/home/gompert-group4/data/timema/hic_genomes/edingburgh/24_0175/Hap2Chr.fasta.masked
t_cris_e_240176h1 /uufs/chpc.utah.edu/common/home/gompert-group4/data/timema/hic_genomes/edingburgh/24_0176/Hap1Chr.fasta.masked
t_cris_e_240176h2 /uufs/chpc.utah.edu/common/home/gompert-group4/data/timema/hic_genomes/edingburgh/24_0176/Hap2Chr.fasta.masked
t_cris_h_gs1 /uufs/chpc.utah.edu/common/home/gompert-group4/data/timema/hic_genomes/t_crist_gs_hap_cen4119/HiRise/Hap1/chroms_final_assembly.fasta.masked
t_cris_h_gs2 /uufs/chpc.utah.edu/common/home/gompert-group4/data/timema/hic_genomes/t_crist_gs_hap_cen4119/HiRise/Hap2/chroms_final_assembly.fasta.masked
t_cris_h_gus1 /uufs/chpc.utah.edu/common/home/gompert-group4/data/timema/hic_genomes/t_crist_gus_hap_cen4280/HiRise/Hap1/chroms_final_assembly.fasta.masked
t_cris_h_gus2 /uufs/chpc.utah.edu/common/home/gompert-group4/data/timema/hic_genomes/t_crist_gus_hap_cen4280/HiRise/Hap2/chroms_final_assembly.fasta.masked
t_cris_r_gs1 /uufs/chpc.utah.edu/common/home/gompert-group4/data/timema/hic_genomes/t_crist_refug_stripe/HiRise/hap1/chroms_final_assembly.fasta.masked
t_cris_r_gs2 /uufs/chpc.utah.edu/common/home/gompert-group4/data/timema/hic_genomes/t_crist_refug_stripe/HiRise/hap2/chroms_final_assembly.fasta.masked
t_cris_r_gus1 /uufs/chpc.utah.edu/common/home/gompert-group4/data/timema/hic_genomes/t_crist_refug_green/HiRise/hap1/chroms_final_assembly.fasta.masked
t_cris_r_gus2 /uufs/chpc.utah.edu/common/home/gompert-group4/data/timema/hic_genomes/t_crist_refug_green/HiRise/hap2/chroms_final_assembly.fasta.masked
```
```bash
#!/bin/bash 
#SBATCH --time=168:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=96
#SBATCH --mem=524000
#SBATCH --account=gompert
#SBATCH --qos=gompert-grn
#SBATCH --partition=gompert-grn
#SBATCH --job-name=cactus
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=zach.gompert@usu.edu

module load cactus

cd /scratch/general/nfs1/u6000989/cactus

cactus jobStore_Prog /uufs/chpc.utah.edu/common/home/gompert-group4/data/timema/hic_genomes/comp_aligns/chr8haplotypes/ch8_cactus.txt cactusTcrAll8.hal --maxCores 72
```
