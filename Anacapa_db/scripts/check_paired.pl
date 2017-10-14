#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
########################
sub suffix_remover($);
################        ###########################        ############################        ###########################
=cut
Should you decide to make this script publicly available, copy the suffix_remover() subroutine into here.
This script takes as input two fastq files with the forward and reverse unmatched mates for a paired-end read data set.
Input: provide on the command line, the paired-end read files (forward and reverse) containing reads that need to be matched.
Output: the forward and reverse mate files will be matched and written to READNAME_sorted.fastq, where
'READNAME' is simply the name of the file provided. The READNAME_singletons.fastq files contain the singleton reads
(reads with no matching mates) for both the input files.

by: Armin PEYMANN
=cut
################        ###########################        ############################        ###########################
sub sequence_separator ($);
sub hash_maker_using_fastq_sequences (\@);
############################################
my $path_read1 = shift @ARGV;
my $path_read2 = shift @ARGV;
my $qx_output_read1 = qx(wc -l $path_read1);
my ($number_of_lines_read1) = split(/\s/, $qx_output_read1);
my $qx_output_read2 = qx(wc -l $path_read2);
my ($number_of_lines_read2) = split(/\s/, $qx_output_read2);
if ($number_of_lines_read1 % 2 != 0 || $number_of_lines_read2 % 2 != 0){
die "The number of lines in one of the files is not dividable by 2.\nFor each sequence in your fastq files, you must have"
        . " the following lines:\n"
        . "\@HEADER\nSEQUENCE\n+QUALITY-HEADER\nQUALITY VALUES\n"
        . "Also make sure there is no empty line at the end of your file.\n";
}


 my @read1 = @{sequence_separator($path_read1)};
 my %read1 = %{hash_maker_using_fastq_sequences(@read1)};
 undef @read1;
 my @read2 = @{sequence_separator($path_read2)};
 my %read2 = %{hash_maker_using_fastq_sequences(@read2)};
 undef @read2;
 my $dir_read1 = dirname($path_read1);
 my $read1_name_without_suffix = suffix_remover($path_read1);
 my $path_out_read1 = $dir_read1 . "/" . $read1_name_without_suffix . "_sorted.fastq";
 open(FH_OUT1, ">$path_out_read1");
 my $dir_read2 = dirname($path_read2);
 my $read2_name_without_suffix = suffix_remover($path_read2);
 my $path_out_read2 = $dir_read2 . "/" . $read2_name_without_suffix . "_sorted.fastq";
 open(FH_OUT2, ">$path_out_read2");    
 my $path_out_read1_singletons = $dir_read1 . "/" . $read1_name_without_suffix . "_singletons.fastq";
 open(FH_OUT3, ">$path_out_read1_singletons");    
 foreach my $key (sort keys %read1){    
if (exists $read2{$key}){    
    print FH_OUT1 $read1{$key};
    print FH_OUT2 $read2{$key};
}else{
    print FH_OUT3 $read1{$key};
}
 }
 close FH_OUT1;
 close FH_OUT2;
 close FH_OUT3;
 print "sorted reads were written to:\n";
 print "Check out $path_out_read1" . "\n";
 print "Check out $path_out_read2". "\n" . "\n";
 my $path_out_read2_singletons = $dir_read2 . "/" . $read2_name_without_suffix . "_singletons.fastq";
 open(FH_OUT4, ">$path_out_read2_singletons");
 foreach my $key (sort keys %read2){
     unless (exists $read1{$key}){
         print FH_OUT4 $read2{$key};
     }
 }
 close FH_OUT4;
 print "single reads with no mate were written to:\n";
 print "Check out $path_out_read1_singletons" . "\n";
 print "Check out $path_out_read2_singletons" . "\n";
 
 sub hash_maker_using_fastq_sequences (\@){
     my $array_ref = shift;
     my @read = @{$array_ref};
     my %read;
    foreach my $sequence (@read){
        my $copy_sequence = $sequence;
        my ($first_header) = split(/\n|\r/, $copy_sequence);
        $first_header =~ /(.+)[# ]/;    # one single space or '#' should cover most of fastq files.
        my $pair_id = $1;
        $read{$pair_id} = $sequence;
        
        }    
    return \%read;
 }
     
sub sequence_separator ($){
    my $path_read = shift;
my $line_counter = 0;
my $event;
my @events;
open(FH_IN, $path_read) or die "unable to open FH_IN1\n";
while(my $line = <FH_IN>){
    $line_counter++;
        
    if ($line_counter <= 4){        
        $event .= $line;
    }
    if ($line_counter == 4){
    push(@events, $event);
    undef $event;
    $line_counter = 0    
    }
  }
  close FH_IN;
  return \@events;
}

sub suffix_remover($){                            #get rid of the .<format> in the file name. (If there are more
    my $pathOfDesiredFile = shift;                    #than one dot in the file name, the shortest part of the file name will be taken!)
    $pathOfDesiredFile =~ /([^\/]+)(?!\/)$/;
    my $desiredFileName = $1 if $1;
    $desiredFileName =~ s/(\..+)(?!\.)// if $desiredFileName =~ /\./ ;                                                         
    return $desiredFileName;    
}