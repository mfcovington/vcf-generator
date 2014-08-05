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
use Log::Reproducible;
use feature 'say';
use Vcf;
use Getopt::Long;

my $par1_id   = "R500";
my $par2_id   = "IMB211";
my $fa        = "B.rapa_genome_sequence_0830.fa";
my $cov_min   = 4;
my $directory = ".";

my $options = GetOptions(
    "par1_id=s"     => \$par1_id,
    "par2_id=s"     => \$par2_id,
    "fa=s"          => \$fa,
    "cov_min=i"     => \$cov_min,
    "directory=s"   => \$directory,
    "single_parent" => \$single_parent,
);

my @snp_files = @ARGV;

my $vcf_out = Vcf->new();
$vcf_out->add_columns($par1_id);
$vcf_out->add_columns($par2_id) unless $single_parent;
$vcf_out->add_header_line(
    {
        key   => 'source',
        value => $0,
    }
);
$vcf_out->add_header_line(
    {
        key   => 'reference',
        value => $fa,
    }
);
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
open my $vcf_fh, ">", "$directory/output.vcf";
print $vcf_fh $vcf_out->format_header();

my %snps;
for my $file (@snp_files) {
    open my $polydb_fh, "<", "$file";
    my $header = <$polydb_fh>;
    while (<$polydb_fh>) {
        my ( $chr, $pos, $ref, $alt, $genotype ) = split /\t/, $_;
        if ( $ref eq 'INS' ) {
            next;    # skipping indels for now
        }
        elsif ( $alt eq 'del' ) {
            next;    # skipping indels for now
        }
        else {
            $snps{$chr}{$pos} = {
                'chr' => $chr,
                'ref' => $ref,
                'alt' => $alt,
                'gen' => $genotype
            };
        }
    }
}

for my $chr ( sort keys %snps ) {
    my $par1_file = "$directory/genotyped/$par1_id.$chr.genotyped.nr";
    my $par2_file = "$directory/genotyped/$par2_id.$chr.genotyped.nr"
        unless $single_parent;
    open my $par1_fh, "<", $par1_file;
    open my $par2_fh, "<", $par2_file unless $single_parent;

    add_sample( $par1_id, $par1_fh );
    add_sample( $par2_id, $par2_fh ) unless $single_parent;

    for my $pos ( sort { $a <=> $b } keys $snps{$chr} ) {
        my %out;
        $out{CHROM}  = $chr;
        $out{POS}    = $pos;
        $out{ID}     = "$snps{$chr}{$pos}{'gen'}.$pos.snp";
        $out{REF}    = $snps{$chr}{$pos}{'ref'};
        $out{QUAL}   = '.';
        $out{FILTER} = ['.'];
        my $par1_depth  = $snps{$chr}{$pos}{"$par1_id.par1"};
        my $par2_depth = $snps{$chr}{$pos}{"$par2_id.par2"}
            unless $single_parent;
        next if $par1_depth < $cov_min;
        next if $par2_depth < $cov_min && !$single_parent;
        my $total_depth
            = $single_parent ? $par1_depth : $par1_depth + $par2_depth;
        $out{INFO} = { DP => $total_depth };
        $out{FORMAT} = [ 'GT', 'DP' ];
        $out{gtypes}{$par1_id}{GT} = get_genotype( $par1_id, $chr, $pos );
        $out{gtypes}{$par2_id}{GT} = get_genotype( $par2_id, $chr, $pos )
            unless $single_parent;
        $out{gtypes}{$par1_id}{DP} = $par1_depth;
        $out{gtypes}{$par2_id}{DP} = $par2_depth unless $single_parent;

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
        next unless exists $snps{$chr}{$pos};
        $count++;
        $snps{$chr}{$pos}{"$sample.par1"} = $par1;
        $snps{$chr}{$pos}{"$sample.par2"} = $par2;
        $snps{$chr}{$pos}{"$sample.tot"}  = $tot;
    }
}

sub get_genotype {
    my ( $id, $chr, $pos ) = @_;
    my $ref = $snps{$chr}{$pos}{'ref'};
    my $alt = $snps{$chr}{$pos}{'alt'};
    return $snps{$chr}{$pos}{'gen'} eq $id ? "$alt/$alt" : "$ref/$ref";
}
