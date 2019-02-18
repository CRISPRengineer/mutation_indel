#!/usr/bin/perl -w
##########################################################################################
# name:			mutation.pl																 
# date:   		17-Oct-2014																 
# author:  		CRISPRengineer															 
# description: 																			
# usage: 		
# ./mutation.pl -i <indir> -b <left bound> -l length									 
# .mutation.pl -b 141 -l 20
# ./mutation.pl -b 119 -l 20
# ./mutation.pl -b 268 -l 20
# 2018-6-28 run
# ./mutation.pl -b 106 -l 25
# 2018-7-3 run
##########################################################################################

#########################################################
# Arguments:
# input directory to SAM files
# -i /path/to/sam/files 
#  
# default set to current working directory
# -5' boundary of reference to begin mutation analysis
# -b <Numeric>
#
# length of mutation analysis in bp's
# -l <Numeric>
#########################################################

our($opt_i, $opt_b, $opt_l); # set up global variables for Getopt

use strict;
use Cwd;

#---------------------------------------------------------------------------------------#
my $start=(times)[0];			# save the start time
#---------------------------------------------------------------------------------------#


#=======================================================#
#================== START Main Script ==================#
#=======================================================#



#---------------------------------------------------------------------------------------#
get_args();			# get the command line arguments 

my $lbound=$opt_b;				#5' boundary = -b input
my $rbound=$lbound+$opt_l;		#5' boundary = -b + -l

# $SAMsummary file handle for writing the summary
open(my $SAMsummary,">$opt_i/summary.txt");
#---------------------------------------------------------------------------------------#

################################################################
# @sam = all sam files from the input directory (default: cwd) #
################################################################


#---------------------------------------------------------------------------------------#
my @sam = @{get_sam($opt_i)};
# print the headers for the summary file
print "sample\ttotal reads\tmapped reads\tdeletion reads\tinsertion reads\tPE indel reads\tPE indel freq\n";
print $SAMsummary "sample\ttotal reads\tmapped reads\tdeletion reads\tinsertion reads\tPE indel reads\tPE indel freq\n";
#---------------------------------------------------------------------------------------#


####################################################
# Main loop through each SAM file in the directory #
####################################################

	
#---------------------------------------------------------------------------------------#		
foreach (@sam) {
	
	$_=~/(\S+)\.sam/;
	my $SAMfile=$1;
	
	# $SAM a filehandle for reading each input sam file
	open(my $SAM,"$opt_i/$_");
	
	# Filehandle: $SAMindel an output sam file of all the reads with indels
	open(my $SAMindel,">$opt_i/$SAMfile.indel.txt");
	
	# Filehandle: $SAMstats an output sam file of all the positions and sizes of all indels
	open(my $SAMstats,">$opt_i/$SAMfile.stats.txt");
	
	# Filehandle: $SAMpos an output sam file of all the positions and sizes of all indels
	open(my $SAMpos,">$opt_i/$SAMfile.pos.txt");
	
	# Filehandle: $SAMpos an output sam file of all the positions and sizes of all indels
	#open(SIZE, ">$opt_i/$SAMfile.size.txt");

	# main subroutine to parse the sam file and generate the outputs	
	read_sam($SAMfile,$SAM,$SAMindel,$lbound,$rbound,$SAMsummary,$SAMstats,$SAMpos);
	
	# close the files for current SAM file before the next loop
	close($SAMpos);
	close($SAMstats);	
	close($SAMindel);
	close($SAM);
}	
#---------------------------------------------------------------------------------------#

# summary file contains overall statistics of all the input sam files 		
close($SAMsummary);
											
#=====================================================#
#================== END Main Script ==================#
#=====================================================#


###############################################################
# get_args parses the commandline arguments using Getopt::Std #
###############################################################


#---------------------------------------------------------------------------------------#
sub get_args{
	use Getopt::Std;
	
	# set in main script:
	# our($opt_i, $opt_b, $opt_l); # set up global variables for Getopt
	# -i = input directory
	# -b = boundary of mutation detection window (5')
	# -l = length of mutation detection window

	getopt('ibl');		# return variables following -i -b -l
	
	if (!defined $opt_i) {$opt_i=getcwd()}
	if (!defined $opt_l) {$opt_l=30}	
}
#---------------------------------------------------------------------------------------#

############################################################
# Get all sam files in the input directory and return them #  
############################################################

#---------------------------------------------------------------------------------------#
sub get_sam {
	my @inputFiles;
	my $dir=shift @_;
	opendir(DIR, $dir) or die $!;
	while (my $file = readdir(DIR)) {
		next unless (-f "$dir/$file");
		if ($file !~ /^\./ & $file =~ /\.sam$/) {
    		push(@inputFiles, $file);
    	}
	}
	return \@inputFiles;
}
#---------------------------------------------------------------------------------------#

#########################
# Write out indel stats #  
#########################


#---------------------------------------------------------------------------------------#
sub indel_stats {
	#$SAMstats,\@TargetDel,\@TargetDelSize,\@TargetIns,\@TargetInsSize
	my $SAMstats=shift;
	my $SAMpos=shift;
	my $mappedreads=shift;
	my @TargetDel = @{$_[0]};
	my @TargetDelSize = @{$_[1]};
	my @TargetIns = @{$_[2]};
	my @TargetInsSize = @{$_[3]};
	
	my %TargetDel;
	my %TargetIns;
	
	my $delReads=scalar @TargetDel;
	my $insReads=scalar @TargetIns;
	
	foreach(@TargetDel){
		$TargetDel{$_}++;
	}
	
	my @keys = sort { $a <=> $b } keys(%TargetDel);
	
	foreach(@TargetIns){
		$TargetIns{$_}++;
	}

	my @Inskeys = sort { $a <=> $b } keys(%TargetIns);	

	for (my $c=0; $c <= $#keys; $c++) {
		#print $SAMpos "$keys[$c]\t$mappedreads\t$delReads\t$insReads\t$TargetDel{$keys[$c]}\t$TargetIns{$keys[$c]}\t".$TargetDel{$keys[$c]}/$delReads."\n";
	}

	#print $SAMpos scalar @TargetDel;
	
	for (my $i=0;$i<=$#TargetDel;$i++){
		if (defined $TargetIns[$i]){
			print $SAMstats "$TargetDel[$i]\t$TargetDelSize[$i]\t$TargetIns[$i]\t$TargetInsSize[$i]\n";
		}	
		
		else{print $SAMstats "$TargetDel[$i]\t$TargetDelSize[$i]\n";}
	}
	%TargetDel=();
	%TargetIns=();
	@keys=();
	@Inskeys=();
	#@vals=();
}
#---------------------------------------------------------------------------------------#

###############################
# Read in all input SAM files #  
###############################

#---------------------------------------------------------------------------------------#
sub read_sam {
	
	my ($SAMfile,$SAM,$SAMindel,$lbound,$rbound,$SAMsummary,$SAMstats,$SAMpos)=@_;
	
	##################################
	# Set up the temporary variables #  
	##################################
	
	my @CIGAR;
	my @tabs;
	my $position;
	my @CIGARlengths;
	my @CIGARoperations;
	
	my @del;
	my @ins;
	my @sizedel;
	my @sizeins;
	
	# if readprint == 1 print the current sam file line to $SAMindel
	my $readprint;
	
	# counts of the reads(total, mapped, deletion, insertions)
	my $totalreads=0;
	my $mappedreads=0;
	my $delreads=0;
	my $insreads=0;	
	my $prevread=0;
	my $PEindel=0;
	my $PEfreq=0;
	
	# arrays holding the position(TargetDel,TargetIns) and the size of all indels 
	my @TargetDel;
	my @TargetDelSize;
	my @TargetIns;
	my @TargetInsSize;
	
	
	############################################
	# loop through all lines of the SAM files: # 
	############################################
	
	while(<$SAM>) {
		if ($_=~/^@/) { 
			print $SAMindel "$_";
		}
		else {
					
			# parse each line of the sam file splitting tabs	
			@tabs = split("\t", $_);
			
			# increment total reads
			if ($prevread=$tabs[0]){
				$totalreads++;	
			}	
			
			if ($tabs[5]!~/\*/){ 		# if CIGARstring is defined the read is aligned
				# increment mapped reads
				if ($prevread=$tabs[0]){
					$mappedreads++;
				}	
				# assign position to tab 4
				$position = $tabs[3];
				# split the CIGAR string into numeric/character divider
				@CIGAR = split(/([MIDNSHP])/, $tabs[5]);
				#loop through each element of CIGAR
				# CIGAR strings consist of either:
				# a) a number of base pairs followed by
				# b) a character indicating alignment status: M(Match), I(Insertion), etc...
				
				foreach my $val (@CIGAR) {
					# if number add to lengths
					if ($val =~/\d+/){
						push(@CIGARlengths, $val);
					}
					# else push to operations
					else {
						push(@CIGARoperations, $val);
					}	
				}
				# loop through each pair of operations and lengths a
				for (my $i=0;$i<=$#CIGARoperations;$i++){
					# M = Match alignment
					if ($CIGARoperations[$i]=~/M/){
						# adjust the position tracker
						$position+=$CIGARlengths[$i];
					}
					# D = Deletion
					elsif ($CIGARoperations[$i]=~/D/){
						# push the position of each deletion
						push(@del, $position);
						# push the size of each deletion
						push(@sizedel, $CIGARlengths[$i]);
						# adjust the position tracker
						$position+=$CIGARlengths[$i];
					}
					# I = Insertion
					elsif ($CIGARoperations[$i]=~/I/){
						# push the location of each deletion
						push(@ins, $position);
						# push the size of each deletion
						push(@sizeins, $CIGARlengths[$i]);
					}	
				}
				# reset the temporary counters for the next round
				$position=0;
				@CIGARoperations = ();
				@CIGARlengths = ();
			
				# loop through all deletions in the read and push those in the target
				for (my $i=0;$i<=$#del;$i++){	
					# test if mutation is within target
					if ($del[$i]>$lbound and $del[$i]<$rbound) {
						# set readprint to output this line
						
						
						if ($prevread=$tabs[0]){
							$readprint=1;
							# increment the number of reads with deletions
							$delreads++;
						}	
						# push the position of the deletion
						push(@TargetDel, $del[$i]);
						# push the size of the deletion
						push(@TargetDelSize, $sizedel[$i]);
						$i=$#del;		# only count 1 deletion per target read
					}
				}				
				for (my $c=0;$c<=$#ins;$c++){
					# test if mutation is within target	
					if ($ins[$c]>$lbound and $ins[$c]<$rbound) {
						if ($prevread=$tabs[0]){
							# set readprint to output this line
							$readprint=1;
							# increment the number of reads with insertions
							$insreads++;
						}	
						# push the position of the insertion
						push(@TargetIns, $ins[$c]);
						# push the size of the insertion
						push(@TargetInsSize, $sizeins[$c]);
						$c=$#ins;		# only count 1 insertion per target read
					}
				}	
				
				if($readprint==1 and $prevread=$tabs[0]){
					print $SAMindel "$_";
					$PEindel++;
				}
				$readprint=0;
			}	
		}	
		# reset the temp counts for each line
		$prevread=$tabs[0];
		@del=();		
		@sizedel=();	
		@ins=();		
		@sizeins=();
		$readprint=0;
	}    
	
	# print the indel position and size data to a tab delimited file
	# indel_stats(\@TargetDel,\@TargetDelSize,\@TargetIns,\@TargetInsSize);
	
	# print the sam file and mapped read as comments to the stats file
	
	#####################
	# Print the headers #
	#####################

	print $SAMstats "## indel position and size of $SAMfile\n";
	print $SAMstats "## Number of mapped reads: $mappedreads\n";
	print $SAMstats "del.pos\tdel.size\tins.pos\tins.size\n";
	
	print $SAMpos "## indel position frequency $SAMfile\n";
	print $SAMpos "## Number of mapped reads: $mappedreads\n";
	print $SAMpos "pos\tmapped.reads\ttotal.del.reads\t\ttotal.ins.reads\tdel.reads\t\tins.reads\tdel.freq\n";
	
	
	indel_stats($SAMstats,$SAMpos,$mappedreads,\@TargetDel,\@TargetDelSize,\@TargetIns,\@TargetInsSize);
	
	######################
	# Display the counts #
	######################
	
	$PEfreq=100*$PEindel/$mappedreads;
	print "$SAMfile\t$totalreads\t$mappedreads\t$delreads\t$insreads\t$PEindel\t$PEfreq\n";
	print $SAMsummary "$SAMfile\t$totalreads\t$mappedreads\t$delreads\t$insreads\t$PEindel\t$PEfreq\n";

		
	$totalreads=0;
	$mappedreads=0;
	$delreads=0;
	$insreads=0;
	$PEindel=0;
	$PEfreq=0;
	
	@TargetInsSize=();
	@TargetDelSize=();
	@TargetDel=();
	@TargetIns=();
} 
#---------------------------------------------------------------------------------------#


################################################
# Calculate and print execution time of script #
################################################

my $end=(times)[0];			# save the end time
my $dt =$end-$start;		# difference is the execution time
print STDERR "Execution time = $dt seconds\n"; # outputs to STDERR