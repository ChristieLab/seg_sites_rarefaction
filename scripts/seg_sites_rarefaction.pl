use List::Util qw/sum/;
use List::Util qw/shuffle/;
use Data::Dumper;

if ($#ARGV == 4) {
  $file = $ARGV[0];
  $file2 = $ARGV[1];
  $file3 = $ARGV[2];
  $mcmcs = $ARGV[3];
  $outfile = $ARGV[4];
} 
else {
  die;
}

open(FILE, "<$file") or die $!;
open(FILE2, "<$file2") or die $!;
open(FILE3, "<$file3") or die $!;
open(OUTFILE, '>', $outfile) or die $!;

my @output;

$header = <FILE>;
chomp($header);
@header = split(/\t/, $header);

my @stat = stat $file;



# read in the subfacet details
my @subfacets;
my $tracker = 0;
while(<FILE3>){
  chomp($_);
  $pop = $_;
  $pop =~ s/\R//g;
  $subfacets[$tracker] = $pop;
  $tracker++;
}

close FILE3;

# init results
my %results;
@pops = uniq(@subfacets);
foreach $pop (@pops){
  @{$results{$pop}} = (0) x $mcmcs;
}

my $locus_tracker = 0;

print("Starting...\n");
while (<FILE>) {
  if($locus_tracker % 10000 == 0){
      print("Loci: ", $locus_tracker, "\n");
  }
  
  #print("=============================================\n");
  
  $genotypes = $_;
  chomp($genotypes);
  
  $g = <FILE2>;
  $g =~ s/\R//g;
  next if($g < 1);
  #print("g = ", $g, "\n");

  

  @genotypes = split(/\t/, $genotypes);
  
  $n = sum(@genotypes);
  
  # figure out our population
  $tpop = $subfacets[$locus_tracker];
  
  my @genoarray = ();
  
  # build sampling array
  $i = 0;
  for($i; $i <= $#genotypes; $i++){
    if($genotypes[$i] != 0){
      #print($i, " filling with ", $genotypes[$i], " ", $header[$i], "'s\n");
      @narray = ($header[$i]) x $genotypes[$i];
      push(@genoarray, @narray);
    }
  }
  
  #foreach(@genoarray){print($_, " ");}
  # sample off shuffled array mcmc times
  $i = 1;
  #print($tpop, "\n");
  for($i; $i <= $mcmcs; $i++){
    #print("\n", $i, "\n");

    my $pg = $g - 1;
    # shuffle
    @ind = shuffle(0..$#genoarray);
    @ind_picks = @ind[0..$pg];  
    @geno_picks = @genoarray[@ind_picks];
    
    #foreach $pick(@geno_picks){
    #  print($pick, "\t");
    #}
    #print("\n");
    
    # unique options
    @geno_picks = uniq(@geno_picks);
    
    
    
    if($#geno_picks > 0){
      # seg
      #print("seg!\n");
      $results{$tpop}[$i - 1]++;
    }
    else{
      $alleles = join(@geno_picks);
      @alleles = split(//, $alleles);
      @alleles = uniq(@alleles);
      if($#geno_picks > 0){
        #print("seg!\n");
        # seg
        $results{$tpop}[$i - 1]++;
      }
    }
  }
  
  $locus_tracker++;
}


close(FILE); close(FILE2);

foreach $key (keys(%results)){
  $pop = $key;
  $pop =~ s/\R//g;
  print OUTFILE $pop , "\t";
  foreach (@{$results{$key}}){
    print OUTFILE $_ , "\t";
  }
  print OUTFILE "\n";
}


close(OUTFILE);

sub uniq {
  my %seen;
  return grep { !$seen{$_}++ } @_;
}