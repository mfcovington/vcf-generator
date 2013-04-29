#!/usr/bin/env perl
# parental-vcf-writer.pl
# Mike Covington
# created: 2013-04-25
#
# Description:
#
use strict;
use warnings;
use autodie;
use feature 'say';
use Vcf;

my $par1_id = "R500";
my $par2_id = "IMB211";

my $vcf_out = Vcf->new();
$vcf_out->add_columns($par1_id);
$vcf_out->add_columns($par2_id);
$vcf_out->add_header_line(
    {
        key         => 'FORMAT',
        ID          => 'GT',
        Number      => '1',
        Type        => 'String',
        Description => "Genotype"
    }
);
$vcf_out->add_header_line(
    {
        key         => 'FORMAT',
        ID          => 'DP',
        Number      => '1',
        Type        => 'String',
        Description => "Depth"
    }
);
$vcf_out->add_header_line(
    {
        key         => 'INFO',
        ID          => 'DP',
        Number      => 1,
        Type        => 'Integer',
        Description => "Total Depth"
    }
);
open my $vcf_fh, ">", "output.vcf";
print $vcf_fh $vcf_out->format_header();

my %db;
my @chromosomes = qw(A01 A02 A03 A04 A05 A06 A07 A08 A09 A10);
for my $chr (@chromosomes) {
    open my $polydb_fh, "<", "snp_master/polyDB.$chr.nr";
    my $header = <$polydb_fh>;
    while (<$polydb_fh>) {
        my ( $chr_id, $pos, $ref, $alt, $genotype ) = split /\t/, $_;
        if ( $ref eq 'INS' ) {
            next;    # skipping indels for now
        }
        elsif ( $alt eq 'del' ) {
            next;    # skipping indels for now
        }
        else {
            $db{$chr}{$pos} = {
                'chr' => $chr,
                'ref' => $ref,
                'alt' => $alt,
                'gen' => $genotype
            };
        }
    }

    my $par1_file = "genotyped/R500.$chr.genotyped.nr";
    my $par2_file = "genotyped/IMB211.$chr.genotyped.nr";
    open my $par1_fh, "<", $par1_file;
    open my $par2_fh, "<", $par2_file;

    add_sample( $par1_id, $par1_fh );
    add_sample( $par2_id, $par2_fh );

    for my $pos ( sort { $a <=> $b } keys $db{$chr} ) {
        my %out;
        $out{CHROM}  = $chr;
        $out{POS}    = $pos;
        $out{ID}     = "$db{$chr}{$pos}{'gen'}.$pos.snp";
        $out{REF}    = $db{$chr}{$pos}{'ref'};
        $out{QUAL}   = '.';
        $out{FILTER} = ['.'];
        my $par1_depth  = $db{$chr}{$pos}{"$par1_id.par1"};
        my $par2_depth  = $db{$chr}{$pos}{"$par2_id.par2"};
        my $total_depth = $par1_depth + $par2_depth;
        $out{INFO} = { DP => $total_depth };
        $out{FORMAT} = [ 'GT', 'DP' ];
        $out{gtypes}{$par1_id}{GT} = get_genotype( $par1_id, $chr, $pos );
        $out{gtypes}{$par2_id}{GT} = get_genotype( $par2_id, $chr, $pos );
        $out{gtypes}{$par1_id}{DP} = $par1_depth;
        $out{gtypes}{$par2_id}{DP} = $par2_depth;

        $vcf_out->format_genotype_strings( \%out );
        print $vcf_fh $vcf_out->format_line( \%out );
    }
}

sub add_sample {
    my ( $sample, $fh ) = @_;
    my $count;
    while (<$fh>) {
        chomp $_;
        my ( $chr, $pos, $par1, $par2, $tot ) = split /\t/, $_;
        next unless exists $db{$chr}{$pos};
        $count++;
        $db{$chr}{$pos}{"$sample.par1"} = $par1;
        $db{$chr}{$pos}{"$sample.par2"} = $par2;
        $db{$chr}{$pos}{"$sample.tot"}  = $tot;
    }
}

sub get_genotype {
    my ( $id, $chr, $pos ) = @_;
    my $ref = $db{$chr}{$pos}{'ref'};
    my $alt = $db{$chr}{$pos}{'alt'};
    return $db{$chr}{$pos}{'gen'} eq $id ? "$alt/$alt" : "$ref/$ref";
}
