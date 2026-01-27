use strict;
use File::Basename;
die("Argument: InFile RefFile OutFile\n") if ((@ARGV < 3) or (@ARGV >3));
	
my $InFQ=$ARGV[0];
my $InQC=$ARGV[1];
my $dir=$ARGV[2];
########## QC ARG[1]
my @files=glob($InQC);
my %ID;
my %reads;
my %length;
my %GC;
my %Q20;
my %Q30;
for my $file(@files)
{
#	print $file."\n";
	open In,"<",$file;
	my $id;
	my $r1;
	my $l1;
	my $gc1;
	my $flag=0;
	my %hash;
	while(<In>){
		chomp;
		if(/^Filename\s*(\w.*)\.fastq/){
				$id=$1;
		}
		if(/^Total Sequences\s*(\d+)/){
				$r1=$1;
		}
		if(/^Sequence length\s*(\d.*)/){
				$l1=$1;
		}
		if(/^\%GC\s*(\d.*)/){
				$gc1=$1;
		}
		if(/^\%GC\s*(\d.*)/){
				$gc1=$1;
		}
		if(/^#Quality/){
			$flag=1;	
			next;
		}
		if($flag==1){
 			my @F=split;
			$hash{$F[0]}=$F[1];	
		}
		if(/^>>END_MODULE/){
			$flag=0;
			next;	
		}
		if(/^#Length/){
			$flag=2;
			open Out,">",$dir.'/'.$id.".LengthCount.txt";	
			printf Out "Length\tCount\n";
			next;
		}
		if($flag==2){
			my @FF=split;
			printf Out "%s\t%s\n",$FF[0],$FF[1];
		}

	}
	close Out;
	$ID{$id}=$id;
	$reads{$id}=$r1;
	$length{$id}=$l1;
	$GC{$id}=$gc1;

	my $all=0;
	my $q20=0;
	my $q30=0;
	$all+=$hash{$_} foreach keys %hash;
	$q20+=$hash{$_} foreach 0..20;
	$q30+=$hash{$_} foreach 0..30;
	$q20=1-$q20/$all;
	$q30=1-$q30/$all;


	$Q20{$id}=$q20;
	$Q30{$id}=$q30;
}




########## FQ Avglen # arg[1]
my @files1=glob($InFQ);
my %Avglen;
for my $file(@files1)
{
#	print $file."\n";
	open In,"<",$file;
	$file=~/.*\/(.*)\.fastq/;
	my $ff=$1;
	my $seq=$reads{$ff};
	my $bp=0;
	my $flag=0;
	my $count_A;
	my $count_T;
	my $count_C;
	my $count_G;
	while(my $line=<In>){
		chomp;

		if($line=~/^@/){
			$flag=1;
			next;
		} 
		if($flag==1)
		{
			$count_A=($line =~ s/A/A/g);
			$count_T=($line =~ s/T/T/g);
			$count_C=($line =~ s/C/C/g);
			$count_G=($line =~ s/G/G/g);
			$bp=$bp+$count_A+$count_T+$count_C+$count_G;
		}
		if($line=~/^+/){
			$flag=0;
			next;
		}


	}
    if(exists($reads{$ff})){
		$Avglen{$ff}=$bp/$seq;
	}
}
open Out1,">",$dir."/Result_Table_Info_Merge.txt";	
printf Out1 "SampleID\tMerge Reads\tLength(bp)\tAvglen(bp)\tGC\tQ20\tQ30\n";
for my $id(sort keys %ID)
	{
		printf Out1 "%s\t%s\t%s\t%s\t%s\t%s\t%s\n",$ID{$id},$reads{$id},$length{$id},$Avglen{$id},$GC{$id},$Q20{$id},$Q30{$id};
	}

close Out1;


