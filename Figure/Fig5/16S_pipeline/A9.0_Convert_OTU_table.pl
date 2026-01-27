use strict;
use File::Basename;
die("Argument: otu_tab\n") if (!(@ARGV == 1));

my $rawFile=shift;
my ($fn, $dir, undef) = fileparse($rawFile, qw/.tab/);
open In,"<",$rawFile;
open Out,">",$dir."/".$fn.".abundance.txt";

my @Indi;
my %Data;
while(<In>)
{
	chomp;
	my @Line= split /\t/;
        if(/^\#OTU\s+ID/)
        {
		@Indi=@Line;
		$Indi[0]="ID";
		pop @Indi;
        }
	elsif(/^\d+/)
        {
		my $ID=$Line[$#Line];
		$ID=~s/\s+//g;
		if($ID=~/g_\w+;s__\w+/)
		{
			$ID=~s/;s__/\./g;
		}
		$ID=~s/\w__//g;
		$ID=~s/;/\|/g;
		$ID=~s/\|+$//g;
		for my $i(1..($#Line-1))
                {
                        $Data{$ID}{$Indi[$i]}=$Data{$ID}{$Indi[$i]}+$Line[$i];
                }
	}
}

for my $i(0..$#Indi)
{
	printf Out "%s\t",$Indi[$i];
}
printf Out "\n";

for my $id(sort keys %Data)
{
	printf Out "%s\t",$id;
	for my $i(1..$#Indi)
	{
		printf Out "%s\t",$Data{$id}{$Indi[$i]};
	}
	printf Out "\n";
}

