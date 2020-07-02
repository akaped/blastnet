#!/usr/bin/env perl

use strict;
use warnings;

#Version_2.2

#PROGRAM NAME
my $wuBlastName = programName('blastp');

#HELP
my $help  = 0;
my $tmp;
($help = 1) if ($#ARGV < 1);
foreach my $field (0..$#ARGV){
  $help = 2 if (($ARGV[$field] eq '-h') or ($ARGV[$field] eq '-help') or ($ARGV[$field] eq '--help'));
  ($tmp = $ARGV[1+$field]) if ($ARGV[$field] eq '-tmp');
}

my $helpMessage = "\nNAME
blastR.pl - It runs WU-blastp using the BLOSUMR scoring matrix\n
SYNOPSIS
wu-blastR.pl database queryFile [blastp options][-pg]\n
DESCRIPTION
   * wu-blastR takes a nucleotidic FASTA file as query input.
   * It converts it into a dinucleotidic sequences.
   * It calls blastp on it using the BLOSUMR as scoring matrix.\n
OPTIONS
   * Use the usual blastp parameters (for instance \"wu-blastR.pl database queryFile W=7 mformat=2\").
   * wu-blastR uses as default -Q=10 -R=10 -T=9999999 -W=5.
   * By default wu-blastR.pl assumes that WU-blastp program name on your computer is just \'blastp\'. If this is not the case you must provide the field -pg with the blastp path.\n
TROUBLESHOOTING
   * Both AB-Blast and WU-Blast packages are provided with their own blastp. If needed it is therefore important that the user properly specifies, with -pg option, the WU/BlastP one.
";
if ($help == 1){
    print "$helpMessage\n\n\n";
    system "$wuBlastName --help\n";
    exit;
}
if ($help == 2){
    print "$helpMessage\n";
    exit;
}

#GLOBAL VARIABLES
my (%translation) = ('A'  => 'A',
		     'C'  => 'C',
		     'AA' => 'D',
		     'AC' => 'E',
		     'AG' => 'F',
		     'G'  => 'G',
		     'AT' => 'H',
		     'CA' => 'I',
		     'CC' => 'K',
		     'CG' => 'L',
		     'CT' => 'M',
		     'GA' => 'N',
		     'GC' => 'P',
		     'GG' => 'Q',
		     'GT' => 'R',
		     'TA' => 'S',
		     'T'  => 'T',
		     'TC' => 'V',
		     'TG' => 'W',
		     'A-' => 'X',
		     'C-' => 'X',
		     'G-' => 'X',
		     'T-' => 'X',
		     'AN' => 'X',
		     'CN' => 'X',
		     'GN' => 'X',
		     'TN' => 'X',
		     'NA' => 'X',
		     'NC' => 'X',
		     'NG' => 'X',
		     'NN' => 'X',
		     'NT' => 'X',
		     'N-' => 'X',
		     '-'  => 'X',
		     'N'  => 'X',
		     'TT' => 'Y',
		     '-A' => '*',
		     '-C' => '*',
		     '-G' => '*',
		     '-T' => '*',
		     '-N' => '*',
		     '--' => '*');
my (%lower) = ( 'a' => 1 ,'c' => 1 ,'g' => 1 ,'t' => 1 );
my %compl = ("A" => "T", "C"=>"G", "G"=>"C" , "T"=>"A" , "a" => "t", "c"=>"g", "g"=>"c","t"=>"a");


#SANITY CHECK
my $blastErrorMessage = `$wuBlastName 2>&1`;
if ($blastErrorMessage =~/BLAST query\/options error: Either a BLAST database or subject sequence\(s\) must be specified/){
  die "Error. You are probably calling NCBI-Blast. Please specify the WU-blastP path on your computer\n";
}


#CREATING BLOSUMR MATRIX (wu format) AND POINTING TO IT
my %matrixHash =(
    ' '       => '     A    R    N    D    C    Q    E    G    H    I    L    K    M    F    P    S    T    W    Y    V    B    Z    X    *' ,
    'A'       => '     0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0' ,
    'R'       => '     0    4   -4   -3    0   -4   -2    0    3   -1   -4   -1    3    0   -1   -1    0   -1    2   -1    0    0    0    0' ,
    'N'       => '     0   -4    3    3    0   -4   -4    0   -4    1   -4   -3   -3   -2   -5    3    0   -1   -3   -2    0    0    0    0' ,
    'D'       => '     0   -3    3    4    0   -2   -6    0   -4    3   -3   -3   -3   -4   -3    1    0   -3   -2   -4    0    0    0    0' ,
    'C'       => '     0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0' ,
    'Q'       => '     0   -4   -4   -2    0    3   -3    0   -3   -3    2    0   -3    3   -4   -4    0    2   -3   -3    0    0    0    0' ,
    'E'       => '     0   -2   -4   -6    0   -3    3    0   -3   -2   -4    2    1   -5    4   -2    0   -2   -2    3    0    0    0    0' ,
    'G'       => '     0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0' ,
    'H'       => '     0    3   -4   -4    0   -3   -3    0    2   -3   -4   -3    1   -6    0   -3    0   -3    1   -2    0    0    0    0' ,
    'I'       => '     0   -1    1    3    0   -3   -2    0   -3    3   -3   -5   -4   -2   -2    3    0   -2   -3   -2    0    0    0    0' ,
    'L'       => '     0   -4   -4   -3    0    2   -4    0   -4   -3    2   -6   -7    3   -2   -3    0    3   -4   -5    0    0    0    0' ,
    'K'       => '     0   -1   -3   -3    0    0    2    0   -3   -5   -6    3   -2   -2    4   -4    0   -3   -3    3    0    0    0    0' ,
    'M'       => '     0    3   -3   -3    0   -3    1    0    1   -4   -7   -2    5   -3    1   -4    0   -3    3   -1    0    0    0    0' ,
    'F'       => '     0    0   -2   -4    0    3   -5    0   -6   -2    3   -2   -3    4   -1   -4    0    3   -4   -3    0    0    0    0' ,
    'P'       => '     0   -1   -5   -3    0   -4    4    0    0   -2   -2    4    1   -1    4   -2    0   -2   -1    4    0    0    0    0' ,
    'S'       => '     0   -1    3    1    0   -4   -2    0   -3    3   -3   -4   -4   -4   -2    2    0   -4   -4   -4    0    0    0    0' ,
    'T'       => '     0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0' ,
    'W'       => '     0   -1   -1   -3    0    2   -2    0   -3   -2    3   -3   -3    3   -2   -4    0    3   -5   -6    0    0    0    0' ,
    'Y'       => '     0    2   -3   -2    0   -3   -2    0    1   -3   -4   -3    3   -4   -1   -4    0   -5    4   -2    0    0    0    0' ,
    'V'       => '     0   -1   -2   -4    0   -3    3    0   -2   -2   -5    3   -1   -3    4   -4    0   -6   -2    4    0    0    0    0' ,
    'B'       => '     0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0' ,
    'Z'       => '     0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0' ,
    'X'       => '     0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0   -1    0' ,
    '*'       => '     0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0' ,
);
my @rows = (' ','A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V','B','Z','X','*');
my $matrixFolderName          = matrixNameGenerator();
my $matrixFolderName_toRemove = $matrixFolderName;
mkdir $matrixFolderName or die "Cannot create blastR matrix directory " . $matrixFolderName ."$!\n";
$ENV{'WUBLASTMAT'} = $matrixFolderName;
$matrixFolderName .= 'aa/';
$ENV{'BLASTMAT'} = $matrixFolderName;
mkdir $matrixFolderName or die "Cannot create blastR matrix directory " . $matrixFolderName ."$!\n";
open (MAT, ">$matrixFolderName"."/BLOSUM62") or die "Cannot create the blastR matrix file $matrixFolderName"."/BLOSUM62 $!\n";
foreach my $row (@rows){
  print MAT $row . $matrixHash{$row} ."\n";
}
close MAT;


#PARSING the command line: BLASTP database queryfile [options]
my (@arrayCommandLine, $inputQueryFile , $databaseName);
my $spy_Q   = 1;
my $spy_R   = 1;
my $spy_T   = 1;
my $spy_W   = 1;
my $spy_tmp = 1;
foreach my $field (0..$#ARGV){
   #detecting the database
  if ($field == 0){
    $databaseName = $ARGV[0];
    next;
  }
  #the input queries
  if ($field == 1){
    $inputQueryFile = $ARGV[1];
    next;
  }
  #skipping the -tmp
  if ($ARGV[$field] eq '-tmp'){
    $spy_tmp = 2;
    next;
  }
  if ($spy_tmp == 2){
    $spy_tmp = 3;
    next;
  }
  #checking -Q , -R , -T
  $spy_Q = 2 if (($ARGV[$field] eq '-Q') or ($ARGV[$field] eq 'Q') or ($ARGV[$field]=~/^-Q=.+$/) or ($ARGV[$field]=~/^Q=.+$/));
  $spy_R = 2 if (($ARGV[$field] eq '-R') or ($ARGV[$field] eq 'R') or ($ARGV[$field]=~/^-R=.+$/) or ($ARGV[$field]=~/^R=.+$/));
  $spy_T = 2 if (($ARGV[$field] eq '-T') or ($ARGV[$field] eq 'T') or ($ARGV[$field]=~/^-T=.+$/) or ($ARGV[$field]=~/^T=.+$/));
  $spy_W = 2 if (($ARGV[$field] eq '-W') or ($ARGV[$field] eq 'W') or ($ARGV[$field]=~/^-W=.+$/) or ($ARGV[$field]=~/^W=.+$/));

  push (@arrayCommandLine,$ARGV[$field]);
 }

#default -Q , -R
push (@arrayCommandLine,"-Q=10")      if ($spy_Q == 1);
push (@arrayCommandLine,"-R=10")      if ($spy_R == 1);
push (@arrayCommandLine,"-T=9999999") if ($spy_T == 1);
push (@arrayCommandLine,"-W=5")       if ($spy_W == 1);
my $commandLine = join( " ", "@arrayCommandLine" );



#CONVERTING the input
my $input          = $inputQueryFile;
my $translatedFile = translateLoop($input , 'top' , 0 , 'tmp');



#WU-BLAST
my $command = "$wuBlastName $databaseName $translatedFile "." $commandLine";
(system "$command") == 0 or die "failed to execute the command $command $!\n";
(system "rm -rf $matrixFolderName_toRemove $translatedFile") == 0 or die ("unable to remove the BLOSUMR matrix folder $matrixFolderName_toRemove or the translated input $translatedFile after that the script is finished $!\n");







###
#FUNCTIONS
sub trad_dino_to_mono {
    my ($dino,$pos) = @_;
    my $lowerCaseSpy = 0;
    ($lowerCaseSpy = 1) if (defined $lower{substr $dino, 0, 1});
    $dino = uc $dino;
    if ($translation{"$dino"}){
      if ($lowerCaseSpy){
	return lc($translation{"$dino"});
      }
      else{
	return $translation{"$dino"};
      }
    }
    else {die "unknown dino \"$dino pos $pos\" - crash!\n"}
}


sub complement{
  my ($string) = @_;
  my ($symbol , $complement);
  my $position = 0;
  while ($string =~ /(\w|\-)/g) {
    $symbol = (substr $string, $position, 1);
    if (defined $compl{$symbol}){$complement .= $compl{$symbol};}
    else{$complement .= $symbol;}
    $position++;
  }
  return $complement;
}

sub fileNameGenerator{
  my ($nameRoot) = @_;
  my $tmp_name_counter = 0;
  my $tmp_name;
  while (!$tmp_name || -f $tmp_name) {
    $tmp_name_counter++;
    $tmp_name = "$nameRoot$$".".$tmp_name_counter";
  }
  return $tmp_name;
}


#TRANSLATE FUNCTION. It produce either the top or both strands translated. It consider either a file or the standard input
sub translateLoop{
  my ($inputFile , $strand , $fromStdin , $outDirection) = @_;

  #Sanity Check
  die "Error. Specify in the script \'top\' or \'both\' as strand for the function translateLoop \n" if (($strand ne 'top') and ($strand ne 'both'));

  #read from file. Define that file as Standard Input
  if ($fromStdin == 0){
    open (STDIN, "< $inputFile") or die ("cannot read the $inputFile file $!\n");
  }

  #choose where to produce the translated file
  my ($rootName, $translatedFileName);
  if (defined $tmp){
    $rootName = "${tmp}/blastRpackage_input.translated";
  }
  elsif ($outDirection eq 'tmp'){
    $rootName = '/tmp/blastRpackage_input.translated';
  }
  elsif ($outDirection eq 'local'){
    $rootName = "$inputFile".".translated";
  }
  $translatedFileName = fileNameGenerator ($rootName);
  open (OUT,">$translatedFileName") or die "cannot create the temporary file $translatedFileName $!\n";

  my $totrad          = "";
  my $totrad_revcompl = "";
  my $seqname         = "";
  my $line            = <STDIN>;
  my ($binucQ , $binucQ_revcompl);
  while ($line) {
    if ($line =~ /^>/) {
      $seqname=$line;
      $seqname=~s/\r\n/\n/g;
      $seqname=~s/\n\r/\n/g;
      $seqname=~s/\r\r/\n/g;
      $seqname=~s/\r/\n/g;
      chomp $seqname;
      print OUT "$seqname\n";
      while (defined ($line = <STDIN>) && ($line !~ /^>/)) {
	$line=~s/\r\n/\n/g;
	$line=~s/\n\r/\n/g;
	$line=~s/\r\r/\n/g;
	$line=~s/\r/\n/g;
	chomp $line;
	$line =~ tr/[RYSWKMBDHVryswkmbdhv]/N/;
	$line =~ tr/U/T/;
	$line =~ tr/u/t/;
	$line =~ s/ //g;
	$totrad .= $line;
	$totrad_revcompl .= complement ($line) if ($strand eq 'both');
      }
      my $pos           = 0;
      my $i             = 0;
      my $trad          = "";
      my $trad_revcompl = "";
      $totrad_revcompl=reverse $totrad_revcompl if ($strand eq 'both');
      while ($totrad =~ /(\w|\-)/g) {
	$i++;
	$binucQ = (substr $totrad, $pos, 2);
	$trad .= trad_dino_to_mono($binucQ,$pos);
	if ($strand eq 'both'){
	  $binucQ_revcompl = (substr $totrad_revcompl, $pos, 2);
	  $trad_revcompl .= trad_dino_to_mono($binucQ_revcompl,$pos) if ($strand eq 'both');
	}
	if ($i == 60) {
	  $trad .= "\n";
	  $trad_revcompl .= "\n" if ($strand eq 'both');
	  $i = 0;
	}
	$pos++;
      }
      print OUT "$trad\n";
      $totrad = "";
      #here print the rev compl
      if ($strand eq 'both'){
	print OUT "${seqname}__rev__\n" ;
	print OUT "$trad_revcompl\n" if ($strand eq 'both');
	$totrad_revcompl="";
      }
    }
    else {$line = <STDIN>}
  }
  close OUT;
  return $translatedFileName;
}

sub matrixNameGenerator{
  my $tmp_matrix_name_counter = 0;
  my $tmp_matrix_name;
  while (!$tmp_matrix_name || -d $tmp_matrix_name) {
    $tmp_matrix_name_counter++;
    if (defined $tmp){
	$tmp_matrix_name = "${tmp}/blastRpackage_matrix$$".".$tmp_matrix_name_counter"."/";
    }
    else{
	$tmp_matrix_name = "/tmp/blastRpackage_matrix$$".".$tmp_matrix_name_counter"."/";
    }
  }
  return $tmp_matrix_name;
}


sub programName {
  my ($default) = @_;
  my ($pgName , @cleanARGV);
  my $spyPg = 1;
  foreach my $field (0..$#ARGV){
    if ($ARGV[$field] eq '-pg'){
      $pgName = $ARGV[1+$field];
      $spyPg = 2;
      next;
    }
    if ($spyPg == 2){
      $spyPg = 3;
      next;
    }
    push (@cleanARGV , $ARGV[$field]);
  }
  if (! defined $pgName){
    $pgName = $default;
  }
  @ARGV = @cleanARGV;
  return $pgName;
}
