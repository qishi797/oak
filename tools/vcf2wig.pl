## transfer TCGA VCF 4.1 to file in wiggle format(TCGA level II data): coordinate + methylation_value(%)
##forward and reverse strand are combined
## author: Yaping Liu  lyping1986@gmail.com


my $input_file_name = $ARGV[0];
my $type = $ARGV[1];

my $use_age = "USAGE: perl vcf2wig.pl input_file_name [CG]";
if($ARGV[0] eq ""){
	print "$use_age\n";
	exit(1);
}
my $cpg_name_output = $input_file_name;
$cpg_name_output =~ s/\.vcf//;
$cpg_name_output = $cpg_name_output.".$type.wig";
open(OUT,">$cpg_name_output") or die;
#my $descript;
if($type eq "CG"){
	$type = "CG"; 
	$descript = "CG"; ##only output homozygous CpG
}
else{
	if($type eq ""){
		$descript = "Cytosine";
		$type = "\\w+"; ##by default, output all sites in VCF file
	}
	else{
		$descript = $type;
	}
}

#my $head_line = "track name=$cpg_name_output type=bedDetail description=\"$descript methylation level\" visibility=3";
my $head_line = "track type=wiggle_0 name=".$input_file_name." description=\"$descript methylation level\"";
print OUT "$head_line\n";


my $chr_flag = "";
my $pre = "";
open(FH,"<$input_file_name") or die;
while(<FH>){
	$line=$_;
	chomp($line);
	next if $line =~ /^#/;
	my @splitin = split "\t", $line;
	next unless ($splitin[6] eq "PASS");
	if($splitin[0] ne $chr_flag){
			$chr_flag = $splitin[0];
			my $chr_head_line = "variableStep chrom=$chr_flag";
			print OUT "$chr_head_line\n";
		}

	if($splitin[9] =~ /:(\d+):($type):(\d+):/){
		my $num_c = $1;
		my $context = $2;
		my $num_t = $3;


		if($num_c + $num_t != 0){
			$methy = $num_c/($num_c + $num_t);
			$methy = sprintf("%.2f",100*$methy);
			my $cor = $splitin[1];
			if($cor eq $pre){
				next;
			}
			my $out_line = "$cor\t$methy";
			print OUT "$out_line\n";
			$pre = $cor;
		}
		
	}
	
}

close(FH);
close(OUT);
print "finished!\n";
