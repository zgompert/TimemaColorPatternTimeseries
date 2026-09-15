#!/usr/bin/perl
#
# bwa aln and samse
#


use Parallel::ForkManager;
my $max = 28;
my $pm = Parallel::ForkManager->new($max);

my $genome  = "/uufs/chpc.utah.edu/common/home/gompert-group4/data/timema/hic_genomes/t_crist_gus_hap_cen4280/HiRise/Hap2/chroms_final_assembly.fasta.masked";

FILES:
foreach $fq (@ARGV){
	$pm->start and next FILES; ## fork
	if ($fq =~ m/(\S+)\.fastq/){
		$ind = "fha2011_$1";
	}
	else {
                        die "Failed to match $file\n";
                }
                system "bwa aln -n 5 -l 20 -k 2 -t 1 -q 10 -f a_$ind".".sai $genome $fq\n";
                system "bwa samse -n 1 -r \'\@RG\\tID:$ind\\tPL:ILLUMINA\\tLB:$ind\\tSM:$ind"."\' -f a_$ind".".sam $genome a_$ind".".sai $fq\n";
           $pm->finish;
        
}

$pm->wait_all_children;
