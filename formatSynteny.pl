#!/usr/bin/perl
#
# reformats block_coords.txt to facilitate plotting

open(IN, "blocks_coords.txt") or die "failed to read the input\n";

## parse file
$sid = 0;
$first = 1;
$blk = 0;
while(<IN>){
	chomp;
	if(m/^Seq/ and $first==1){
		$sid = 1;
	} 
	elsif(m/^-/){
		$sid = 0;
		$blk = 0;
		$first = 0;
	}
	elsif($sid==1){ ## in seq ID section
		@line = split(/\s+/,$_);
		$ss{$line[0]} = $line[2];
	}
       	elsif(m/^Block\s+\#(\d+)/){
		$bid = $1;
		$blk = 1;
	} 
	elsif(m/^\d+/){
		@line = split(/\s+/,$_);
		$aln{$bid}{$line[0]}{'strand'} = $line[1];
		$aln{$bid}{$line[0]}{'start'} = $line[2];
		$aln{$bid}{$line[0]}{'end'} = $line[3];
	}
}
close(IN);
open(OUT, "> syn_blocks.txt") or die "failed to write\n";

## print header
print OUT "block";
foreach $sp (sort keys %ss){
	$sid = $ss{$sp};
	print OUT " $sid.strand $sid.start $sid.end";
}
print OUT "\n";

## print synteny blocks
foreach $bid (sort keys %aln){
	print OUT "$bid";
	foreach $sp (sort keys %ss){
		if(defined $aln{$bid}{$sp}){
			print OUT " $aln{$bid}{$sp}{'strand'}";
			print OUT " $aln{$bid}{$sp}{'start'}";
			print OUT " $aln{$bid}{$sp}{'end'}";
		}
		else{ 
			print OUT " NA NA NA";
		}
	}
	print OUT "\n";
}
close(OUT);
