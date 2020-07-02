#!/usr/bin/env perl

use strict;
use warnings;

#Version_2.2

#PROGRAM NAME
my $xdformatName = programName('xdformat');

#HELP
my $help = 0;
my $tmp;
foreach my $field (0..$#ARGV){
  $help = 2 if (($ARGV[$field] eq '-h') or ($ARGV[$field] eq '-help') or ($ARGV[$field] eq '--help'));
  ($tmp = $ARGV[1+$field]) if ($ARGV[$field] eq '-tmp');
}
$help = 1 if (scalar @ARGV < 1);

my $helpMessage = "\nNAME
xdformatR.pl - create databases for blastR.pl in XDF (eXtended Database Format) from one or more input nucleotidic files in FASTA format;\n
SYNOPSIS
xdformatR.pl [options] database\n
DESCRIPTION
   * It takes as input one (or more) nuleotidic multi FASTA files
   * It converts it (them) into a dinucleotidic sequences appending the reverse complement of the input.
   * It creates an AB-blast or WU-blast formatted database running xdformat on the dinucleotidic database.\n
OPTIONS
   * Use the usual xdformat parameters (for instance \"xdformatR.pl -o xdbname [options] fadb\")
   * Since xdformatR works on protein-like sequences the parameter [-n|-p] is unnecessary and therefore skipped if given in the commandline.
   * The options: -a -r -i -V -X are not supported. Please use the normal xdformat to use them.
   * By default xdformatR assumes that xdformat program name on your computer is just \'xdformat\'. If this is not the case you must provide the field -pg with the xdformat path.
   * WARNING: Specify the proper xdformat, either the one included in AB-blast or WU-blast. WU-Blast is unable to search sequence databases that were created by AB-xdformat\n";
if ($help == 1){
    print "$helpMessage\n\n\n";
    system "$xdformatName --help\n";
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
my %lower = ( 'a' => 1 ,'c' => 1 ,'g' => 1 ,'t' => 1 );
my %compl = ("A" => "T", "C"=>"G", "G"=>"C" , "T"=>"A" , "a" => "t", "c"=>"g", "g"=>"c","t"=>"a");


#CHECK if STDIN
my ($standardInput , @cleanARGV);
foreach my $element (@ARGV){
  if ($element eq '-'){
    $standardInput = 1;
    next;
  }
  push (@cleanARGV , $element);
}
@ARGV = @cleanARGV;
if ((scalar @ARGV < 1) && (! defined $standardInput)){
    system ("$xdformatName");
    exit;
}


#PARSE THE COMMAND LINE: taking the input database(s)
my (@arrayCommandLine, @inputDatabaseFiles , $commandLineEnd);
my %flagsWithArgument = (
			 "-M"  => 1 ,
			 "-t"  => 1 ,
			 "-v"  => 1 ,
			 "-d"  => 1 ,
			 "-C"  => 1 ,
			 "-K"  => 1 ,
			 "-T"  => 1 ,
			 "-o"  => 1
			);
my $currentField = $#ARGV;
while ($ARGV[$currentField])  {
  if ((-e $ARGV[$currentField]) and (! defined $flagsWithArgument{$ARGV[$currentField-1]})) {
    push (@inputDatabaseFiles, $ARGV[$currentField]);
  }
  else{
    $commandLineEnd = $currentField;
    last;
  }
  $currentField--;
  last if ($currentField == -1);
}

#SANITY CHECK 1: RETURN XDFORMAT ERROR avoiding unizializated variables if not found inputDatabaseFiles
#(defined $flagsWithArgument{$ARGV[-2]})
if (((! -e $ARGV[-1]) && (! defined $standardInput)) or ( (defined $ARGV[-2]) && (! defined $standardInput) && (defined $flagsWithArgument{$ARGV[-2]}))) {
  system "$xdformatName";
  die "\nFATAL:  Missing input database name(s).\nSYNOPSIS: xdformat [-p|-n] [options] fadb\n";
}

#PARSE THE COMMAND LINE: taking the commandline
my $spy_o   = 1;
my $spy_tmp = 1;
my $outputName;
if (defined $commandLineEnd){
  foreach my $field (0..$commandLineEnd){
    #exiting the program for not suported options
    notSupportedOptions($ARGV[$field]);

    #skipping the -p and -n on the command line
    if (($ARGV[$field] eq '-p') or ($ARGV[$field] eq '-n')){
      next;
    }
    #-o option
    if ($ARGV[$field] eq '-o'){
      push (@arrayCommandLine, $ARGV[$field]);
      $spy_o = 2;
      next;
    }
    if ($spy_o == 2){
      push (@arrayCommandLine, $ARGV[$field]);
      $spy_o = 3;
      $outputName = $ARGV[$field];
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
    #taking the the options
    push (@arrayCommandLine, $ARGV[$field]);
  }
}
my $commandLine;
if (@arrayCommandLine){
  $commandLine = join( " ", "@arrayCommandLine" );
}
$outputName = $inputDatabaseFiles[0] unless (defined $outputName);


#SANITY CHECK 2: FOR STDIN or MULTIPLE INPUT MUST HAVE A COSTUMER OUTPUT NAME
if (((scalar(@inputDatabaseFiles) > 1) or (defined $standardInput)) and ($spy_o != 3)){
  system "$xdformatName";
  die "\nFATAL:  When multiple input files are to be processed, the basename for the
        output database files must be specified with the -o or -a options.\n";
}



#TRANSLATING the multifasta database in protein-like sequences
my $translatedFile;
if (scalar(@inputDatabaseFiles) > 1 or ((defined $standardInput) and (scalar(@inputDatabaseFiles) == 1))) {
    my $concatenatedInput = fileNameGenerator("${outputName}_concatenated");
    open (CONC,">$concatenatedInput") or die "cannot create the temporary concatenated input file $!\n";
    foreach my $inputDatabaseFile (@inputDatabaseFiles){
	open (IN_DB, "<$inputDatabaseFile") or die "cannot read $inputDatabaseFile $!\n";
	foreach my $line (<IN_DB>){
	    print CONC "$line";
	}
	close IN_DB;
    }
    if ($standardInput){
	foreach my $st_in (<STDIN>){
	    print CONC "$st_in";
	}
    }
    close CONC;
    $translatedFile = translateLoop("$concatenatedInput" , 'both'  , 0 , 'local');
    (system "rm $concatenatedInput") == 0 or die "unable to remove the concatenated input $concatenatedInput $!\n" if (-e "$concatenatedInput");
}
elsif ((! @inputDatabaseFiles) and (defined $standardInput)){
    my $input = $outputName;
    $translatedFile = translateLoop($input , 'both'  , 1 , 'local');
}
else{
    my $input = $inputDatabaseFiles[0];
    $translatedFile = translateLoop($input , 'both'  , 0 , 'local');
}



#FORMATTING the translated database
my ($command, $inputName);

if (defined $commandLine){
  $command = "$xdformatName -p $commandLine"." $translatedFile";
}
else{
  $command = "$xdformatName -p $translatedFile";
}
(system "$command") == 0 or die "failed to execute the command $command $!\n";



#CLEANING
(system "rm $translatedFile") == 0 or die "unable to remove the translated input $translatedFile after that the script is finished $!\n";
outputRenamer($inputDatabaseFiles[0]);







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
  my ($rootName , $translatedFileName);
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
      $totrad_revcompl = reverse $totrad_revcompl if ($strand eq 'both');
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
	print OUT "$trad_revcompl\n";
	$totrad_revcompl="";
      }
    }
    else {$line = <STDIN>}
  }
  close OUT;
  return $translatedFileName;
}



sub notSupportedOptions{
  my ($feature) = @_;
  if ($feature eq '-a'){
    die "Aborting. '-a' is not an option supported by xdformatR.pl. Please use the normal xdformat to do it:\n ".
      "   Append sequences to an existing database:
    xdformat [-p|-n] -a xdbname [options] fadb...";
  }
  if ($feature eq '-r'){
    die "Aborting. '-r' is not an option supported by xdformatR.pl. Please use the normal xdformat to do it:\n ".
      "   Report the contents of existing database(s) to stdout in FASTA format:
    xdformat [-p|-n] -r [options] xdbname...";
  }
  if ($feature eq '-i'){
    die "Aborting. '-i' is not an option supported by xdformatR.pl. Please use the normal xdformat to do it:\n ".
      "   Describe the contents of existing database(s):
    xdformat [-p|-n] -i xdbname...";
  }
  if ($feature eq '-V'){
    die "Aborting. '-V' is not an option supported by xdformatR.pl. Please use the normal xdformat to do it:\n ".
      "   Verify the integrity of existing database(s):
    xdformat [-p|-n] -V xdbname...";
  }
  if ($feature eq '-X'){
    die "Aborting. '-X' is not an option supported by xdformatR.pl. Please use the normal xdformat to do it:\n ".
      "   Index or re-index the sequence identifiers in existing database(s):
    xdformat [-p|-n] -X xdbname...";
  }
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
