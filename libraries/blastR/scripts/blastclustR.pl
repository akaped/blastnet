#!/usr/bin/env perl

use strict;
use warnings;

#Version_2.2

#PROGRAM NAME
my $blastClustName = programName('blastclust');

#HELP
my $help  = 0;
my ($spy_i , $standardInput , $tmp);
foreach my $field (0..$#ARGV){
  $help = 2 if (($ARGV[$field] eq '-h') or ($ARGV[$field] eq '--help') or ($ARGV[$field] eq '-help'));
  ($spy_i = 1) if ($ARGV[$field] eq '-i');
  ($tmp = $ARGV[1+$field]) if ($ARGV[$field] eq '-tmp');
  if (($ARGV[$field-1] !~/^-/) and ($ARGV[$field]!~/^-/ )){
      my $element = $field + 1;
      die "[blastclustR.pl] ERROR: Arguments must start with \'-\' (the offending argument \#$element was: \'$ARGV[$field]\')\n";
  }
}
if (! defined $spy_i){
  $standardInput = 1;
}

$help = 1 if ($#ARGV == -1);

my $helpMessage = "\nNAME
blastclustR.pl - blastclust using the BLOSUMR scoring matrix\n
SYNOPSIS
blastclustR.pl (-i multiFASTAfile | STDIN) [blastclust options][-pg]\n
DESCRIPTION
   * BlastclustR takes a nucleotidic multi FASTA file as input (fileName or stdin).
   * It converts it into a dinucleotidic sequences.
   * It calls blastclust on it using the BLOSUMR as scoring matrix.\n
OPTIONS
   * Use the usual blastclust parameters (for instance \"blastclustR.pl -i multiFASTAfile -c configurationFile -o outputFile\").
   * Since blastclustR works on protein-like sequences the blastclust parameter \"-p\" is unnecessary and therefore skipped if given in the commandline.
   * The option \"-S 2\" in the configuration file (if any) is not supported.
   * The flag \"-U\" masks lower case nucleotides displaying them with an X.
   * BlastclustR uses as default -S 0 -L 0.
   * BlastclustR uses as default a temporary configuration file setting -G 6 -E 2.
   * By default blastclustR assumes that blastclust program name on your computer is just \'blastclust\'.
     If this is not the case you must provide the field -pg with the blastclust path.
";
if ($help == 1){
    print "$helpMessage\n\n\n";
    system "$blastClustName --help\n";
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



#CREATING BLOSUMR MATRIX (ab or ncbi format) AND POINTING TO IT
my %matrixHash =(
    ' '       => '   A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  J  Z  X  *' ,
    'A'       => '   0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0' ,
    'R'       => '   0  4 -4 -3  0 -4 -2  0  3 -1 -4 -1  3  0 -1 -1  0 -1  2 -1  0  0  0  0  0' ,
    'N'       => '   0 -4  3  3  0 -4 -4  0 -4  1 -4 -3 -3 -2 -5  3  0 -1 -3 -2  0  0  0  0  0' ,
    'D'       => '   0 -3  3  4  0 -2 -6  0 -4  3 -3 -3 -3 -4 -3  1  0 -3 -2 -4  0  0  0  0  0' ,
    'C'       => '   0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0' ,
    'Q'       => '   0 -4 -4 -2  0  3 -3  0 -3 -3  2  0 -3  3 -4 -4  0  2 -3 -3  0  0  0  0  0' ,
    'E'       => '   0 -2 -4 -6  0 -3  3  0 -3 -2 -4  2  1 -5  4 -2  0 -2 -2  3  0  0  0  0  0' ,
    'G'       => '   0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0' ,
    'H'       => '   0  3 -4 -4  0 -3 -3  0  2 -3 -4 -3  1 -6  0 -3  0 -3  1 -2  0  0  0  0  0' ,
    'I'       => '   0 -1  1  3  0 -3 -2  0 -3  3 -3 -5 -4 -2 -2  3  0 -2 -3 -2  0  0  0  0  0' ,
    'L'       => '   0 -4 -4 -3  0  2 -4  0 -4 -3  2 -6 -7  3 -2 -3  0  3 -4 -5  0  0  0  0  0' ,
    'K'       => '   0 -1 -3 -3  0  0  2  0 -3 -5 -6  3 -2 -2  4 -4  0 -3 -3  3  0  0  0  0  0' ,
    'M'       => '   0  3 -3 -3  0 -3  1  0  1 -4 -7 -2  5 -3  1 -4  0 -3  3 -1  0  0  0  0  0' ,
    'F'       => '   0  0 -2 -4  0  3 -5  0 -6 -2  3 -2 -3  4 -1 -4  0  3 -4 -3  0  0  0  0  0' ,
    'P'       => '   0 -1 -5 -3  0 -4  4  0  0 -2 -2  4  1 -1  4 -2  0 -2 -1  4  0  0  0  0  0' ,
    'S'       => '   0 -1  3  1  0 -4 -2  0 -3  3 -3 -4 -4 -4 -2  2  0 -4 -4 -4  0  0  0  0  0' ,
    'T'       => '   0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0' ,
    'W'       => '   0 -1 -1 -3  0  2 -2  0 -3 -2  3 -3 -3  3 -2 -4  0  3 -5 -6  0  0  0  0  0' ,
    'Y'       => '   0  2 -3 -2  0 -3 -2  0  1 -3 -4 -3  3 -4 -1 -4  0 -5  4 -2  0  0  0  0  0' ,
    'V'       => '   0 -1 -2 -4  0 -3  3  0 -2 -2 -5  3 -1 -3  4 -4  0 -6 -2  4  0  0  0  0  0' ,
    'B'       => '   0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0' ,
    'J'       => '   0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0' ,
    'Z'       => '   0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0' ,
    'X'       => '   0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 -1  0' ,
    '*'       => '   0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0' ,
);
my @rows = (' ','A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V','B','J','Z','X','*');
my $matrixFolderName = matrixNameGenerator();
mkdir $matrixFolderName or die "Cannot create blastclustr matrix directory " . $matrixFolderName ."$!\n";
open (MAT, ">$matrixFolderName"."/BLOSUM62") or die "Cannot create the blastclustr matrix file $matrixFolderName"."/BLOSUM62 $!\n";
foreach my $row (@rows){
  print MAT $row . $matrixHash{$row} ."\n";
}
close MAT;

$ENV{'BLASTMAT'} = $matrixFolderName;



#PARSING the command line
my (@arrayCommandLine , $inputMultiFastaFile , $tmp_configFile_name , $configurationFile , $customOut);
my $spyInput = 1;
my $spy_p    = 1;
my $spy_S    = 1;
my $spy_L    = 1;
my $spy_U    = 1;
my $spy_W    = 3;
my $spy_tmp  = 1;
my $spy_o    = 1;

foreach my $field (0..$#ARGV){
  #taking the input
  if ($ARGV[$field] eq '-i'){
    $inputMultiFastaFile = $ARGV[1+$field];
    $spyInput = 2;
    next;
  }
  if ($spyInput == 2){
    $spyInput = 3;
    next;
  }
  #skipping the -p
  if ($ARGV[$field] eq '-p'){
    $spy_p = 2;
    next;
  }
  if ($spy_p == 2){
    $spy_p = 3;
    next;
  }
  #skipping the -o
  if ($ARGV[$field] eq '-o'){
    $spy_o = 2;
    next;
  }
  if ($spy_o == 2){
    $spy_o = 3;
    $customOut = $ARGV[$field];
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
  #taking the -U
  if ($ARGV[$field] eq '-U'){
    $spy_U = 2;
    next;
  }
  #configuration file -c
  if ($ARGV[$field] eq '-c'){
    $configurationFile = $ARGV[1+$field];
  }
  #checking -S , -L , -W
  $spy_S = 2 if ($ARGV[$field] eq '-S');
  $spy_L = 2 if ($ARGV[$field] eq '-L');
  $spy_W = $ARGV[1+$field] if ($ARGV[$field] eq '-W');

  push (@arrayCommandLine,$ARGV[$field]);
}

#default -S , -L
push (@arrayCommandLine,"-S 0") if ($spy_S == 1);
push (@arrayCommandLine,"-L 0") if ($spy_L == 1);
my $commandLine = join( " ", "@arrayCommandLine" );


#create deafault configuration file
$tmp_configFile_name = tempConfigFileGenerator() unless (defined $configurationFile);


#STRAND CHOICE. CONFIGURATION FILE ANALYSIS
my $S_option;
if (defined $configurationFile){
  open (CONFIG, "<$configurationFile") or die "Cannot open the configuration file $configurationFile $!\n";
  foreach my $line (<CONFIG>){
    if ($line=~/-S\s+(\d+)/){
      $S_option = $1;
      die "Aborting. -S 2 option not supported. Choose -S 1 to search the top strands or -S 3 to search both strands.[default 3]\n" if ($S_option == 2);
    }
    if ($line=~/-W\s+(\d+)/){
      $spy_W = $1;
    }
  }
}


#CONVERTING the input
my ($input , $translatedFile , @filteredOut);
if (! defined $standardInput){$input = $inputMultiFastaFile;}
else{$input = 'input';}
#default
if (! defined $S_option){
  if (! defined $standardInput){
    $translatedFile = translateLoop($input , 'both' , 0 ,'tmp');
  }
  else{
    $translatedFile = translateLoop($input ,'both' , 1 ,'tmp');
  }
}
else{
  if ($S_option == 1){
    if (! defined $standardInput){
      $translatedFile = translateLoop($input ,'top' , 0 ,'tmp');
    }
    else{
      $translatedFile = translateLoop($input , 'top' , 1 ,'tmp');
    }
  }
  if ($S_option == 3){
    if (! defined $standardInput){
      $translatedFile = translateLoop($input , 'both' , 0 ,'tmp');
    }
    else{
      $translatedFile = translateLoop($input , 'both' , 1 ,'tmp');
    }
  }
}


#BLASTCLUST
my $command = "$blastClustName -i $translatedFile -p T"." $commandLine";
$command .= " -c $tmp_configFile_name" unless (defined $configurationFile);


print STDERR "#command: $command\n";
if ($spy_U == 2){
  print STDERR "#-U enabled lower case filtering\n";
}
if (defined $tmp_configFile_name){
  print STDERR "#$tmp_configFile_name: -G 6 -E 2\n";
}
else {
  open (CONFIG, "<$configurationFile") or die "Cannot open the configuration file $configurationFile $!\n";
  my $lineConfig;
  foreach my $line (<CONFIG>){
    chomp $line;
    $lineConfig .= $line;
  }
  if (defined $lineConfig){
    print STDERR "#$configurationFile: $lineConfig\n";
  }
  else{
    print STDERR "#$configurationFile: ";
  }
}

my @blastclustOut = `$command 2>&1`;
if ($?){
    foreach my $el (@blastclustOut){
	print "$el";
    }
    die ("blastclust didn't work for the command $command:$!\n");
}

#OUTPUT CLEANING
my %allLines;
push(@blastclustOut , @filteredOut);
if (defined $customOut){
    system ("rm $customOut") if (-e $customOut);
    open (CUSTOM_OUT, ">$customOut") or die ("cannot write on the $customOut file $!\n");
}
foreach my $line (@blastclustOut){
   my @lineArray = split(/\s+/, $line);

   foreach my $member (@lineArray){
     $member=~s/__rev__$//;
   }
   my $lineToPrint = $line;
   @lineArray    = sort(@lineArray);
   $line         = join(" ",@lineArray);

   if ((! defined $allLines{$line}) or (scalar(@lineArray) == 1)){
       if (defined $customOut){
	   print CUSTOM_OUT "$lineToPrint";
       }
       else{
	   print "$lineToPrint";
       }
   }
   $allLines{$line} = 1;
}
if (defined $customOut){
    close CUSTOM_OUT; 
}

(system "rm -rf $matrixFolderName $translatedFile") == 0 or die ("unable to remove the BLOSUMR matrix folder $matrixFolderName or the translated input $translatedFile after that the script is finished $!\n");
if (defined $tmp_configFile_name){
  (system "rm $tmp_configFile_name") == 0 or die ("unable to remove the temporary configuration file $tmp_configFile_name after that the script is finished $!\n");
}





####
#FUNCTIONS

sub trad_dino_to_mono {
    my ($dino,$pos) = @_;
    my $lowerCaseSpy = 0;
    ($lowerCaseSpy = 1) if (defined $lower{substr $dino, 0, 1});
    $dino = uc $dino;
    if ($translation{"$dino"}){
      if ($lowerCaseSpy){
	if ((defined $spy_U) and ($spy_U == 2)){   #WARNING: variable defined just for blastclustR
	  return 'X';
	}
	else{
	  return lc($translation{"$dino"});
	}
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
	next if ($line=~/^\s*$/);
	$totrad .= $line;
	$totrad_revcompl.= complement ($line) if ($strand eq 'both');
      }
      my $pos           = 0;
      my $trad          = "";
      my $trad_revcompl = "";
      $totrad_revcompl=reverse $totrad_revcompl if ($strand eq 'both');
      while ($totrad =~ /(\w|\-)/g) {
	$binucQ = (substr $totrad, $pos, 2);
	$trad .= trad_dino_to_mono($binucQ,$pos);
	if ($strand eq 'both'){
	  $binucQ_revcompl = (substr $totrad_revcompl, $pos, 2);
	  $trad_revcompl .= trad_dino_to_mono($binucQ_revcompl,$pos) if ($strand eq 'both');
	}
	$pos++;
      }
      if (seqeunceControl($trad) == 1){
	print OUT "$seqname\n";
	print OUT "$trad\n";
      }
      else{
	  if ($seqname=~/>(.*)/){
	      push(@filteredOut,$1."\n");
	  }
      }
      $totrad = "";
      #here print the rev compl
      if ($strand eq 'both'){
	if (seqeunceControl($trad_revcompl) == 1){
	  print OUT "${seqname}__rev__\n" ;
	  print OUT "$trad_revcompl\n";
	}
	else{
	  if ($seqname=~/>(.*)/){
	      push(@filteredOut,"$1"."__rev__\n");
	  }
	}
	$totrad_revcompl = "";
      }
    }
    else {$line = <STDIN>}
  }
  close OUT;
  return $translatedFileName;
}

sub tempConfigFileGenerator{
    my $tmp_name_counter = 0;
    my $tmp_name;
    while (!$tmp_name || -f $tmp_name) {
	$tmp_name_counter++;
      	$tmp_name = "config$$.$tmp_name_counter";
      }
    open (CONFIG, ">/tmp/$tmp_name") or die "cannot create the temporary blastclust configuration file. $!\n";
    print CONFIG "-G 6 -E 2";
    close CONFIG;
    $tmp_name = "/tmp/$tmp_name";
    return $tmp_name;
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

sub seqeunceControl{
    my ($seq) = @_;
    my $switch = 0;
    $switch = 1 if (($seq=~/^\S*[^X]\S*\S$/) && (length($seq) >= ($spy_W * 2)));
    return $switch;
}
