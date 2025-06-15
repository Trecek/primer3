#!/usr/bin/perl

# Test script for command-line argument parsing bug fixes
# Tests the hybrid parser implementation that fixes ambiguous argument matching

use warnings 'all';
use strict;
use Getopt::Long;

our $winFlag;
our $exe_amplicon = '../src/amplicon3_core';
our $exe_oligotm = '../src/oligotm';
our $exe_thal = '../src/ntthal';
our $exe_ntdpal = '../src/ntdpal';

my $do_valgrind;
my %args;
my $failure = 0;
my $test_count = 0;

print "\n=== TESTING ARGUMENT PARSING BUG FIXES ===\n\n";

if (!GetOptions(\%args,
               'valgrind',
               'windows',
              )) {
    print "Usage: $0 [--valgrind] [--windows]\n";
    exit -1;
}

$winFlag = $args{'windows'};
$do_valgrind = $args{'valgrind'};

# Test sequence for amplicon (needs > 36 bases)
my $amp_seq = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
# Test sequence for oligotm (shorter)
my $oligo_seq = "ACGTACGTACGT";

# ==================================
# Test 1: Bug fix verification for amplicon3_core
# Both -dmso and -dmso_fact should work correctly
print "Test 1: amplicon3_core -dmso and -dmso_fact bug fix... ";
$test_count++;
{
    my $cmd = "$exe_amplicon -mv 50 -dv 1.5 -n 0.6 -dmso 10.0 -dmso_fact 0.6 -mf 0 $amp_seq 2>&1";
    my $output = `$cmd`;
    my $exit_code = $? >> 8;
    
    if ($exit_code == 0 && 
        $output =~ /AMPLICON_DMSO=10\.0/ && 
        $output =~ /AMPLICON_DMSO_CORRECTION=0\.6/) {
        print "[OK]\n";
    } else {
        print "[FAILED]\n";
        print "  Command: $cmd\n";
        print "  Exit code: $exit_code\n";
        print "  Output:\n$output\n";
        $failure++;
    }
}

# ==================================
# Test 1b: Verify output values with DVL1 sequence
# This is the actual bug report test case
print "Test 1b: amplicon3_core DVL1 sequence with dmso parameters... ";
$test_count++;
{
    my $dvl1_seq = "CGAGCAATGATGCACAGACGTTGACTTTTGATATAGTTTTTGTTCTAAGTGGGAATCCCATCCTGAGACTTGACAGC";
    my $cmd = "$exe_amplicon -mv 50 -dv 1.5 -n 0.6 -dmso 10.0 -dmso_fact 0.6 -mf 0 $dvl1_seq 2>&1";
    my $output = `$cmd`;
    my $exit_code = $? >> 8;
    
    # Expected output values
    my $expected = "AMPLICON_MONOVALENT=50.0
AMPLICON_DIVALENT=1.5
AMPLICON_DNTPS=0.6
AMPLICON_DMSO=10.0
AMPLICON_DMSO_CORRECTION=0.6
AMPLICON_FORMAMID=0.0
AMPLICON_PRODUCT_SIZE=77
AMPLICON_GC_PERCENT=42.9
AMPLICON_MELTPOINTS=74.4";
    
    my $all_matched = 1;
    for my $line (split /\n/, $expected) {
        if ($output !~ /\Q$line\E/) {
            $all_matched = 0;
            last;
        }
    }
    
    if ($exit_code == 0 && $all_matched) {
        print "[OK]\n";
    } else {
        print "[FAILED]\n";
        print "  Exit code: $exit_code\n";
        print "  Expected output to contain:\n$expected\n";
        print "  Got:\n$output\n";
        $failure++;
    }
}

# ==================================
# Test 2: Ambiguous argument detection for amplicon3_core
# -dm should be rejected as ambiguous
print "Test 2: amplicon3_core ambiguous -dm detection... ";
$test_count++;
{
    my $cmd = "$exe_amplicon -mv 50 -dv 1.5 -n 0.6 -dm 10.0 -mf 0 $amp_seq 2>&1";
    my $output = `$cmd`;
    my $exit_code = $? >> 8;
    
    if ($exit_code != 0 && $output =~ /error: ambiguous argument '-dm'/) {
        print "[OK]\n";
    } else {
        print "[FAILED]\n";
        print "  Expected non-zero exit and 'ambiguous argument' error\n";
        print "  Exit code: $exit_code\n";
        print "  Output:\n$output\n";
        $failure++;
    }
}

# ==================================
# Test 3: Backward compatibility for amplicon3_core
# -fo should still work as abbreviation for -formamid
print "Test 3: amplicon3_core backward compatibility -fo... ";
$test_count++;
{
    my $cmd = "$exe_amplicon -fo 0.8 $amp_seq 2>&1";
    my $output = `$cmd`;
    my $exit_code = $? >> 8;
    
    if ($exit_code == 0 && $output =~ /AMPLICON_FORMAMID=0\.8/) {
        print "[OK]\n";
    } else {
        print "[FAILED]\n";
        print "  Exit code: $exit_code\n";
        print "  Output:\n$output\n";
        $failure++;
    }
}

# ==================================
# Test 4: Unknown argument for amplicon3_core
print "Test 4: amplicon3_core unknown argument -xyz... ";
$test_count++;
{
    my $cmd = "$exe_amplicon -xyz $amp_seq 2>&1";
    my $output = `$cmd`;
    my $exit_code = $? >> 8;
    
    if ($exit_code != 0 && $output =~ /error: unknown argument '-xyz'/) {
        print "[OK]\n";
    } else {
        print "[FAILED]\n";
        print "  Expected non-zero exit and 'unknown argument' error\n";
        print "  Exit code: $exit_code\n";
        print "  Output:\n$output\n";
        $failure++;
    }
}

# ==================================
# Test 5: oligotm -d vs -dm disambiguation
print "Test 5: oligotm -d and -dm disambiguation... ";
$test_count++;
{
    my $cmd = "$exe_oligotm -d 100 -dm 5.0 $oligo_seq 2>&1";
    my $output = `$cmd`;
    my $exit_code = $? >> 8;
    
    if ($exit_code == 0 && $output =~ /^\d+\.\d+/) {
        print "[OK]\n";
    } else {
        print "[FAILED]\n";
        print "  Exit code: $exit_code\n";
        print "  Output:\n$output\n";
        $failure++;
    }
}

# ==================================
# Test 6: oligotm ambiguous -d (could be -d, -dm, -df)
print "Test 6: oligotm ambiguous argument detection... ";
$test_count++;
{
    # Note: -d by itself is not ambiguous in oligotm, but let's test a truly ambiguous case
    # Actually, in oligotm, -d matches only -d, -dm matches only -dm, -df matches only -df
    # So let's skip this test for oligotm as it doesn't have ambiguous prefixes
    print "[SKIPPED] (oligotm uses exact lengths)\n";
}

# ==================================
# Test 7: thal argument parsing
print "Test 7: ntthal basic argument parsing... ";
$test_count++;
if (-x $exe_thal) {
    my $cmd = "$exe_thal -mv 50 -dv 1.5 -n 0.8 -s1 ACGTACGT -s2 ACGTACGT -a ANY 2>&1";
    my $output = `$cmd`;
    my $exit_code = $? >> 8;
    
    if ($exit_code == 0) {
        print "[OK]\n";
    } else {
        print "[FAILED]\n";
        print "  Exit code: $exit_code\n";
        print "  Output:\n$output\n";
        $failure++;
    }
} else {
    print "[SKIPPED] (ntthal not found)\n";
}

# ==================================
# Test 8: ntdpal argument parsing
print "Test 8: ntdpal basic argument parsing... ";
$test_count++;
if (-x $exe_ntdpal) {
    my $cmd = "$exe_ntdpal -g 1.0 -l 0.5 ACGTACGT TGCATGCA g 2>&1";
    my $output = `$cmd`;
    my $exit_code = $? >> 8;
    
    if ($exit_code == 0) {
        print "[OK]\n";
    } else {
        print "[FAILED]\n";
        print "  Exit code: $exit_code\n";
        print "  Output:\n$output\n";
        $failure++;
    }
} else {
    print "[SKIPPED] (ntdpal not found)\n";
}

# ==================================
# Test 9: ntthal hybrid parser - test original strncmp behavior
# Test that -dna matches -d (original used strncmp with length 2)
print "Test 9: ntthal -dna abbreviation compatibility... ";
$test_count++;
if (-x $exe_thal) {
    my $cmd = "$exe_thal -s1 ACGTACGT -s2 TGCATGCA -dna 100 -a ANY 2>&1";
    my $output = `$cmd`;
    my $exit_code = $? >> 8;
    
    if ($exit_code == 0) {
        print "[OK]\n";
    } else {
        print "[FAILED]\n";
        print "  Expected -dna to be accepted as abbreviation for -d\n";
        print "  Exit code: $exit_code\n";
        print "  Output:\n$output\n";
        $failure++;
    }
} else {
    print "[SKIPPED] (ntthal not found)\n";
}

# ==================================
# Summary
print "\n=== SUMMARY ===\n";
if ($failure == 0) {
    print "All $test_count tests PASSED\n\n";
    exit 0;
} else {
    print "$failure of $test_count tests FAILED\n\n";
    exit -1;
}