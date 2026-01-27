use strict;
use File::Basename;
die("Argument: InFile RefFile OutFile\n") if ((@ARGV < 2) or (@ARGV >2));
	

my $In3=$ARGV[0];
my $dir=$ARGV[1];


########## In3 Avglen Clean Tags
my @files3=glob($In3);
my %Avglen;
my %count3;
my %length;
for my $file(@files3)
{
#	print $file."\n";
#	if($file=~/3_filtered_fa\/(.*gut.*fa)/){
#		next
#	}
	open In,"<",$file;
	my $ff;
	if($file=~/3_filter_fq\/(.*)\.fa/){
		$ff=$1;
	}
	print $ff."\n";
#	if($file=~/3_filtered_fa\/(.*)\.sub/){
#		$ff=$1;
#	}
#	if($file=~/3_filtered_fa\/(.*)\.join/){
#		$ff=$1;
#	}
	my $cc=0;
	my $bp=0;
	my $count_A;
	my $count_T;
	my $count_C;
	my $count_G;
	my $temp;
	while(my $line=<In>){
		chomp;
		if($line=~/^>/){
			$cc=$cc+1;
		} else {
			$count_A=($line =~ s/A/A/g);
			$count_T=($line =~ s/T/T/g);
			$count_C=($line =~ s/C/C/g);
			$count_G=($line =~ s/G/G/g);
			$temp=$count_A+$count_T+$count_C+$count_G;
			$bp=$bp+$temp;
			if(exists($length{$temp})){
				$length{$temp}=$length{$temp}+1;
			} else {
				$length{$temp}=1;
			}
		}
	}
	$Avglen{$ff}=$bp/$cc;
	$count3{$ff}=$cc;
}
open Out1,">",$dir."/Result_Table_Info_End.txt";	
printf Out1 "SampleID\tSubsample Tags\tAvglen_Sub(bp)\n";
for my $id(sort keys %count3)
	{
		printf Out1 "%s\t%s\t%s\n",$id,$count3{$id},$Avglen{$id};
	}

close Out1;
################### length 分布统计
open Out2,">",$dir."/Result_Length_End.txt";		
printf Out2 "Length\tNum\n";
for my $id(sort keys %length)
	{
		printf Out2 "%s\t%s\n",$id,$length{$id};
	}

close Out2;

