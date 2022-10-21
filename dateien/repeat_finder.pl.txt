#!/usr/local/bin/perl

## Program Info:
#
# Name: repeat_finder
#
# Purpose:
#   Searches for DNA repeats in a sequence. It searches for every
#   permutation of an n-mer, and reports the number of hits.
#
# Author: John Nash
#
# Copyright (c) National Research Council of Canada, 2002,
#   all rights reserved.
#
# Licence: This script may be used freely as long as no fee is charged
#   for use, and as long as the author/copyright attributions
#   are not removed.
#
# History:
#   4-12 July 2002: Alpha versions.
#
#
# Thanks: To Ed Taboada for suggesting hashing the repeats
#
# To-do:  Put ranges in for command line switches
##

use strict;
use warnings;
use Text::Wrap;
$Text::Wrap::columns = 65;

my $time_zero = time;

# who am I and where am I?
my $title = "repeat_finder";
my $version = "1.13";
my $date = "17 July, 2002";
my $error_msg = "Type \"$title -h\" for help.";
print "$title $version (released $date)\n";

#  Command line switches:
#    $nmer_query - size of the DNA sequence to be searched for repeats.
#    $modifier - how many times larger than expected the nmer
#      has to be present in order to be flagged and reported.
#    $rank_hits - number of relevantly flagged hits to display.
#    $sort_key - sort by number of hits or more than expected.
#    $coords - do we display co_ordinates.
#    $circular - is the sequence circular.

# Get and process the command line params:
# Returns array of $fasta_file and $orf_file;
my ($nmer_query, $rank_hits, $sort_key, $coords, $circular) = 
	process_command_line();

##  Handle the input sequence.
## Does the input sequence exist:  handle errors
# If $ARGV[0] is not blank, test for file's existence:
if (defined $ARGV[0]) {
	unless (-e $ARGV[0]) {
		die("\nError: Sequence file \'$ARGV[0]\' does *not* exist. \n", 
				$error_msg, "\n");
	}
	open (FILE, $ARGV[0]);
}

# If it has come in from a redirection or pipe, 
#  check it is bigger than 0:
else {
	my $fh = *STDIN;
	unless ((-p $fh) or (-s $fh)) {
		die("\nError: Piped sequence file does *not* exist. \n", 
				$error_msg, "\n");
	}
	*FILE = *STDIN;
}


# Command-line range checks:
if ($nmer_query <=2) {
	die("\nError: Nmer query-size has to be greater than 2 \n",
			$error_msg, "\n");
}

$rank_hits = 2 if ($rank_hits < 2);

## Read in the sequence from a FASTA file:
#   For multiple sequences, concatenate them with ">".
#   There is no reason for multiple sequences to be thus analysed
#   unless they are contigs from a single project.

my ($seq_name, $seq_length, $seq_str1, $seq_str2, $count);
# read the header:
$count = 0;
while (<>) {
	s/\r\n/\n/g;
	chomp;
	if (/^>/)  {
		$seq_name = substr($_, 1, length $_);
		$seq_str1 .= ">" if ($count > 0);
		$count++;
	}
	else {
		$seq_str1 .= uc $_;
	}
}

# Some final sequence processing:
$seq_length = length $seq_str1;
$seq_str1 = uc $seq_str1;

# Eliminate all nonACGTs:
$seq_str1 =~ tr/XNATGCBDKRVHMYSW/XXATGCXXXXXXXXXX/;
$seq_str2 = reverse($seq_str1);
$seq_str2 =~ tr/XNATGCBDKRVHMYSW/XXTACGXXXXXXXXXX/;

## Hash search routine:
# Search for each nmer string by hashing it from the sequence:
print STDERR "** Hashing all ",$nmer_query,"-mers...";

# Holds the hits:
my %nmers_t;
my $k;

# First strand:
print STDERR "  First strand...";

# Account for circularity (add $nmer's of sequence from the 
# start to the end)
$seq_str1 .= substr $seq_str1, 0, $nmer_query if ($circular eq "yes");

for (1..((length $seq_str1) - $nmer_query)) {
	$k = substr $seq_str1, $_, $nmer_query;
# remove Ns, Xs and the "between sequence delimiter" '>':
	next if ($k =~ /X/);
	next if ($k =~ />/);
	$nmers_t{$k}++;
}

# Second strand:
# reverse-complement the sequence:
print STDERR "  Second strand...\n";
$seq_str2 .= substr $seq_str2, 0, $nmer_query if ($circular eq "yes");

for (1..((length $seq_str2) - $nmer_query)) {
	$k = substr $seq_str2, $_, $nmer_query;
# remove Ns, Xs and the "between sequence delimiter" '>':
	next if ($k =~ /X/);
	next if ($k =~ />/);
	$nmers_t{$k}++;
}

print STDERR "** Found elements: ", scalar (keys %nmers_t), ";";
print STDERR " Theoretical answer: ", 4**$nmer_query, "\n";

# Calculating incomings ACGT content:
my ($As, $Cs, $Gs, $Ts, $totalACGT, $fractA, $fractC, $fractG, $fractT);
$As = ($seq_str2 =~ tr/A//);
$Cs = ($seq_str2 =~ tr/C//);
$Gs = ($seq_str2 =~ tr/G//);
$Ts = ($seq_str2 =~ tr/T//);
$totalACGT = $As + $Cs + $Gs + $Ts;

# Convert to 2 decimal places:
$fractA = sprintf("%0.2f", $As/$totalACGT);
$fractC = sprintf("%0.2f", $Cs/$totalACGT);
$fractG = sprintf("%0.2f", $Gs/$totalACGT);
$fractT = sprintf("%0.2f", $Ts/$totalACGT);

# Expected hits calculations:
my (%expected_hits, %times_expected);
my $exp_hits;

print STDERR "** Expected hits and observed hits...\n";
foreach (keys %nmers_t) {
	$exp_hits = expected_hits($_) * $seq_length * 2;
	if ($exp_hits > 10000) { 
		$exp_hits = sprintf ("%0.1e", $exp_hits);
	}
	elsif ($exp_hits < 0.1)  { 
		$exp_hits = sprintf ("%0.1e", $exp_hits);
	}
	else {
		$exp_hits = sprintf ("%0.1f", $exp_hits);
	}
	$expected_hits{$_} = $exp_hits;
	
	$times_expected{$_} = $nmers_t{$_}/$exp_hits;
	if ($times_expected{$_} > 10000) { 
		$times_expected{$_} = sprintf ("%0.1e", $times_expected{$_}); 
	}
	elsif ($times_expected{$_} < 0.1)  { 
		$times_expected{$_} = sprintf ("%0.1e", $times_expected{$_}); 
	}
	else { 
		$times_expected{$_} = sprintf ("%0.1f", $times_expected{$_});
	}
}

# Flags zero hits to report "No hits found":
my $count_0 = 0;

# Expected hits are calculated based on the GC content of the
#   incoming sequence:

# Sort the hash by number of hits:
my @sorted_hits;
print STDERR "** Sorting by $sort_key...\n";
$sort_key eq "number_of_hits" ?
	(@sorted_hits = sort sort_by_number_of_hits(keys %nmers_t)) :
	(@sorted_hits = sort sort_by_most_overrepresented (keys %nmers_t));

# Select the number of hits desired by the user:
@sorted_hits = @sorted_hits[0..($rank_hits-1)];


# Output the data:
print STDERR "** Printing...\n";
# Basic information:
# Keep the following line, it's used as a regex start point!"

print wrap '', '', "Sequence: [",
	(($count > 1) ? "Multiple entries" : $seq_name),"] ($seq_length bp) ", 
	($circular eq "yes"? "circular": "linear"), "molecule\n";

# If the incoming sequence has some degenerate bases. Warn user.
if ($totalACGT < $seq_length) {
	print "The analysed sequence has ", $seq_length-$totalACGT, 
	" degenerate bases.  This may skew the results.\n";
}

printf "\nThe top %d %d-mers sorted by %s:\n",
	$rank_hits, $nmer_query, $sort_key;

print "Rank:\tHits:\tExp:\t\tX:\t\tSequence:\n";
print "----\t----\t---\t\t-\t\t--------\n";

foreach (@sorted_hits) {
	$count_0++;
	printf "%d\t%d\t%-8s\t%-8s\t%s\n",
	$count_0,
	$nmers_t{$_},
	exists $expected_hits{$_} ? $expected_hits{$_}:0,
	exists $times_expected{$_}? $times_expected{$_}:0,
	$_;
} # end of foreach

print "No hits found\n" if ($count_0 == 0);
print "\n";

# Map the position of each required hit:
if ($coords eq "yes") {
	print "Co-ordinates: \n";
	foreach (@sorted_hits) {
		my @coords;
# forward strand:
		print "$_ (forward_strand): ";
		while ($seq_str1 =~ m/(?=$_)/g) {
	    push @coords, (pos $seq_str1) + 1;
		}
		print scalar @coords, " hits\n";
		print wrap ('', '', "@coords\n");
		if (scalar @coords == 0) {
			print "0\n";
		}
		
# reverse strand:
		@coords = ();
		print "$_ (reverse_strand): ";
		while ($seq_str2 =~ m/(?=$_)/g) {
	    push @coords, (pos $seq_str2) + 1;
		}
		print scalar @coords, " hits\n";
		print wrap ('', '', "@coords\n");
		if (scalar @coords == 0) {
			print "0\n";
		}
		print "\n";
	}
}

# Final time stamp:
printf STDERR "** Analysis complete: %d seconds (%0.2f minutes) elapsed\n", 
	time - $time_zero, (time - $time_zero)/60;

#<--- end of main function

sub expected_hits {
## Function: Calculates predicted frequency of a sequence
#  based on ACGT content
#  Expects: a string of nucleotides
#  Returns: an expected value of hits.
	
	my $incoming = $_[0];
	my $As = ($incoming =~ tr/A//);
	my $Cs = ($incoming =~ tr/C//);
	my $Gs = ($incoming =~ tr/G//);
	my $Ts = ($incoming =~ tr/T//);
	
	my $outgoing = ($fractA**$As) * ($fractC**$Cs) * 
		($fractG**$Gs) * ($fractT**$Ts);
	return $outgoing;
}

sub sort_by_number_of_hits {
	$nmers_t{$b} <=> $nmers_t{$a}
	or $a cmp $b;
}

sub sort_by_most_overrepresented {
	$times_expected{$b} <=> $times_expected{$a}
	or $a cmp $b;
}

sub process_command_line {
# Variables:
	my %opts = ();    # command line params, as entered by user
	my @cmd_line;     # returned value
	my @list;         # %opts as an array for handling
	my $cmd_args;	    # return value for getopts()
	my $item;
	
# Holders for command line's files:
	my $nmer_query = 9;
	my $rank_hits = 20;
	my $sort_key = "most_overrepresented";
	my $coords = "no";
	my $circular = "yes";
	
# Get the command=line parameters:
	use vars qw($opt_n $opt_r $opt_S $opt_O $opt_L $opt_h);
	use Getopt::Std;
	$cmd_args = getopts('n:r:SOLh', \%opts);
	
# Die on illegal argument list:
	if ($cmd_args == 0) {
		die ("Error: Missing or incorrect command line parameter(s)!\n",
				 $error_msg, "\n");
	}
	
# Make the hashes into an array:
  @list = keys %opts;
	
# Do a quick check for "help" and the compulsory parameters:
#   If the compulsory files are not there, squawk and die:
	foreach $item (@list)  {
# Help:
		if ($item eq "h")  {
	    help();
		}
# nmer_query:
		elsif ($item eq "n") {
	    $nmer_query = $opts{$item};
		}
# rank_hits:
		elsif ($item eq "r") {
	    $rank_hits = $opts{$item};
		}
# sort_key
		elsif ($item eq "S") {
	    $sort_key = "number_of_hits";
		}
# coords:
		elsif ($item eq "O") {
	    $coords = "yes";
		}
# circular:
		elsif ($item eq "L") {
	    $circular = "no";
		}
	}
	
# Put it in an array:
	@cmd_line = ($nmer_query, $rank_hits, $sort_key, $coords, $circular);
	return @cmd_line;
	
} #end of sub process_command_line()

sub help {
	
print <<EOHelp;
	
Function:  Searches for oligonucleotides which are over-represented in a
   genomic sequence.  These could be \"DNA uptake sequences\" or other
   repeats.  N-mers of any size can be searched and the results can be
   sorted by number of hits or by over-representation.  Over-representation
   is calculated by dividing the number of observed Hits by the
   expected abundance (taking %GC into account).
    
   The results will be in "pairs" of complementary sequences because
   both DNA strands will be searched.  Co-ordinates are not really for human
   eyes to look at - rather, they should be passed on to other computer
   programs, or even Excel, for further analysis.  Co-ordinates are based upon
   the relevant strand (i.e. co-ordinate 2 from the reverse strand is the
   reverse complement of the second last base of the forward strand).
	
Syntax:  $title -nX -rX -L -O -S fasta.file
   or    $title -h for help

Command-line parameters:
   -nX where X is the number or Nmers to be searched (default: 9)
   -rX where X is the number of best hits to be returned (default: 20)

Command-line toggles:
   -L switches Linear sequence ON (default: assumes circular sequence)
   -O toggles display co-Ordinates (default: don\'t display co-ordinates)
   -S toggles Sort by number of hits (default: sort by over-representation)

Examples:
   \"$title -n9 -r22 -LOS < fasta_file\"
   will search \"fasta_file\", as a linear sequence, for 9-mers, and display
   the best 20 hits (and their co-ordinates), sorted by number of hits.

Warning:  No attempt is made to ensure that the input file is a valid
   FASTA file.  The program will accept a FASTA file containing multiple
   sequences, but it assumes that they are from the same project and
   concatenates them (with delimiters to keep them apart) for the repeat
   analysis, and reports on the total set.

\"$title\" will work as a pipe, i.e. will accept a redirected FASTA file,
   and output can be redirected to a file.  The \"status comments\"
   (starting with **) are sent to the computer screen and will not be
   redirected to the file.

Warning: $title is a serious memory hog, and should not be considered
   for searching for repeats of greater than 12 mers without at least 
   one megabtye of RAM, and preferable 4!!!

If this scrolls by too fast type: \"$title \| more\".
EOHelp
die ("\n");
} # end of sub help

