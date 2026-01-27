use strict;
use File::Basename;
die("Argument: InFile OutFile Grouping_File Class_col SubClass_col\n") if ((@ARGV < 3) or (@ARGV >5));
if(@ARGV == 4) 
{
	my $InFile=$ARGV[0];
	my $OutFile=$ARGV[1];
	my $GroupFile=$ARGV[2];
	my $grpCol=$ARGV[3]-1;

	my ($fn, $dir, undef) = fileparse($InFile, qw/.txt/);

	my %Groups;
	my $level;
	open Ref,"<",$GroupFile;
	my $count=0;
	while(<Ref>)
	{
		chomp;
        	if($count ==0)
        	{
                	my @temp=split /\t/;
	                $level=$temp[$grpCol];
			$count++
	        }
 	       else
		{
                	my @temp=split /\t/;
                	$Groups{$temp[0]}=$temp[$grpCol];
			$count++;
        	}
	}

	open In,"<",$InFile;
	open Out,">",$OutFile;
	my $FirstLine=<In>;
	my @line1=split /\t/,$FirstLine;
	my %gr;
	my $gr_c=0;
	for my $id(@line1){
		$gr{$gr_c}=$id;
		$gr_c++;
	}

	printf Out "Class\t";
	for my $i(1..$#line1)
	{
		if((exists $Groups{$line1[$i]}) ){
			printf Out "%s\t",$Groups{$line1[$i]};
		}
	}
	printf Out "\n";
	for my $i(0..$#line1)
	{
		if($i==0){
			printf Out "%s\t",$gr{$i};
		} else {
			
			if((exists $Groups{$gr{$i}}) ){
				printf Out "%s\t",$gr{$i};
			}
		}
	}
	printf Out "\n";
	while(<In>)
	{
		my @line2=split /\t/;
		for my $i(0..$#line2)
		{
			if($i==0){ 
				printf Out "%s\t",$line2[$i];
			} else {
			
				if((exists $Groups{$gr{$i}}) ){
					printf Out "%s\t",$line2[$i];
				}
			}
		}
		printf Out "\n";	
	}
}

if(@ARGV == 5)
{
        my $InFile=$ARGV[0];
        my $OutFile=$ARGV[1];
        my $GroupFile=$ARGV[2];
        my $class=$ARGV[3]-1;
	my $subclass=$ARGV[4]-1;

        my %Groups;
        my $L1;
	my $L2;
        open Ref,"<",$GroupFile;
        my $count=0;
        while(<Ref>)
        {
                if($count ==0)
                {
                        my @temp=split /\s+/;
                        $L1=$temp[$class];
			$L2=$temp[$subclass];
                        $count++
                }
               else
                {
                        my @temp=split /\s+/;
                        $Groups{$temp[0]}{"class"}=$temp[$class];
			$Groups{$temp[0]}{"subclass"}=$temp[$class].".".$temp[$subclass];
                        $count++;
                }
        }
	open In,"<",$InFile;
        open Out,">",$OutFile;
        my $FirstLine=<In>;

        my @temp =split /\s+/,$FirstLine;
        printf Out "Class\t";
        for my $i(1..$#temp)
        {
                printf Out "%s\t",$Groups{$temp[$i]}{"class"};
        }
        printf Out "\n";
	printf Out "SubClass\t";
        for my $i(1..$#temp)
        {
                printf Out "%s\t",$Groups{$temp[$i]}{"subclass"};
        }
        printf Out "\n";
        printf Out "%s",$FirstLine;



        while(<In>)
        {
                printf Out "%s",$_;
        }
}
