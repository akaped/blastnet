#!/usr/bin/env perl

use strict;
use warnings;

#Version_2.2

#PROGRAM NAME
my $formatdbName = programName('../blast/formatdb');

#HELP
my $help                   = 0;
my ($standardInput , $tmp);

foreach my $field (0..$#ARGV){
  $help = 2 if (($ARGV[$field] eq '-h') or ($ARGV[$field] eq '-help') or ($ARGV[$field] eq '--help'));
  ($tmp = $ARGV[1+$field]) if ($ARGV[$field] eq '-tmp');
}
$help = 1 if (scalar @ARGV < 1);

my $helpMessage = "\nNAME
formatdbR.pl - create databases for NCBI-BLAST from one or more input nucleotidic files in FASTA format;\n
SYNOPSIS
perl formatdbR.pl [formatdb options]\n
DESCRIPTION
   * It takes as input one (or more) nuleotidic multi FASTA files.
   * It converts it (them) into a dinucleotidic sequences apending the reverse complement of the input.
   * It creates an ncbi-blast formatted database running formatdb on the dinucleotidic database.\n
OPTIONS
   * Use the usual formatdb parameters (for instance \"formatdbR.pl -i inputDatabase -n outDatabaseName\")
   * Since formatdbR works on protein-like sequences the parameter \"-p\" is unnecessary and therefore skipped if given in the commandline
   * By default formatdb assumes that formatdb program name on your computer is just \'formatdb\'. If this is not the case you must provide the field -pg with the formatdb path.";

if ($help == 1){
    print "$helpMessage\n\n\n";
    system "$formatdbName --help\n";
    exit;
}
if ($help == 2){
    print "$helpMessage\n\n";
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


#PARSING THE COMMAND LINE
my (@arrayCommandLine, $inputDatabaseFile , $userDefinedOutput);
my $spyInput   = 1;
my $spy_p      = 1;
my $spy_tmp    = 1;
foreach my $field (0..$#ARGV){
  #the inputs
  if (($ARGV[$field] eq '-i') and ($spyInput == 1)){
    $spyInput = 2;
     next;
  }
  if ($spyInput == 2){
    $inputDatabaseFile = $ARGV[$field];
    $spyInput = 1;
    next;
  }
  #skip the type
  if ($ARGV[$field] eq '-p'){
    $spy_p = 2;
    next;
  }
  if ($spy_p == 2){
    $spy_p = 3;
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
  #user output name. Important to understand where to produce the translated temporary file to optimize the process
  if ($ARGV[$field] eq '-n'){
      $userDefinedOutput = 1;
  }

  push (@arrayCommandLine,$ARGV[$field]);
 }
my $commandLine;
if (@arrayCommandLine){
  $commandLine = join( " ", "@arrayCommandLine" );
}


#SANITY CHECK
die "[formatdb] ERROR: No database name was specified\n" unless ($inputDatabaseFile);
if ($inputDatabaseFile eq "stdin"){
  $standardInput = 1;
}



#TRANSLATING the multifasta database in protein-like sequences
my ($input , $translatedFile);
if (! defined $standardInput){
    $input = $inputDatabaseFile;
}
else{$input = 'stdin';}

#input file
if (! defined $standardInput){
    if (defined $userDefinedOutput){
	$translatedFile = translateLoop($input , 'both' , 0 , 'tmp');
    }
    else{
	$translatedFile = translateLoop($input , 'both' , 0 , 'local');
    }
}
#stdin
else{
    if (defined $userDefinedOutput){
	$translatedFile = translateLoop($input , 'both' , 1 , 'tmp');
    }
    else{
	$translatedFile = translateLoop($input , 'both', 1 , 'local');
    }
}



#RUN FORMATDB
my ($command, $inputName);
if (defined $commandLine){
  $command = "$formatdbName -i $translatedFile "."-p T" . " $commandLine";
}
else{
  $command = "$formatdbName -i $translatedFile "."-p T";
}
(system $command) == 0 or die ("formatdb did not work for command " . $command . ": ". $!);



#CLEANING
(system "rm $translatedFile") == 0 or die "cannot remove the temporary translated file $translatedFile $!\n";
outputRenamer($inputDatabaseFile);






###
#FUNCTIONS

sub outputRenamer{
  my ($nameRoot) = @_;
  my @allTranslatedFiles = `ls ${translatedFile}* 2>/dev/null`;
  if (@allTranslatedFiles){
    foreach my $currentTranslatedFile (@allTranslatedFiles){
      chomp $currentTranslatedFile;
      if ($currentTranslatedFile =~/(\.\w\w\w)$/){
	my $extension = $1;
	(system "mv $currentTranslatedFile ${nameRoot}$extension") == 0 or die "cannot rename the output file $currentTranslatedFile into ${nameRoot}$extension after that the script is finished$!\n";
      }
    }
  }
}

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

#TRANSLATE FUNCTION. It produce either the top or both strans translated. It consider either a file or the standard input
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
  else{
      if ($outDirection eq 'tmp'){
	  $rootName = '/tmp/blastRpackage_input.translated';
      }
      if ($outDirection eq 'local'){
	  $rootName = "$inputFile".".translated";
      }
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
	next if ($line=~/^\s*$/);
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
	print OUT "${seqname}_rev\n" ;
	print OUT "$trad_revcompl\n" if ($strand eq 'both');
	$totrad_revcompl="";
      }
    }
    else {$line = <STDIN>}
  }
  close OUT;
  return $translatedFileName;
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
