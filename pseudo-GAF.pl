#!/usr/bin/env perl 
use strict;
use warnings;

## Change this to whatever taxon you are working with
my $taxon = 'taxon:1000';
chomp(my $date = `date +%Y%M%d`);

my (%aspect, %gos);
## Read the GO.terms_and_ids file to get the aspect (sub ontology)
## of each GO term. 
open(my $fh, $ARGV[0]) or die "Need a GO.terms_and_ids file as 1st arg: $!\n";
while (<$fh>) {
    next if /^!/;
    chomp;
    my @fields = split(/\t/);
    ## $aspect{GO:0000001} = 'P'
    $aspect{$fields[0]} = $fields[2];
}
close($fh);

## Read the list of gene annotations
open($fh, $ARGV[1]) or die "Need a list of gene annotattions as 2nd arg: $!\n";
while (<$fh>) {
    chomp;
    my ($gene, @terms) = split(/[\s;]+/);
    ## $gos{geneA} = (go1, go2 ... goN)
    $gos{$gene} = [ @terms ];
}
close($fh);

foreach my $gene (keys(%gos)) {
    foreach my $term (@{$gos{$gene}}) {
        ## Warn and skip if there is no aspect for this term
        if (!$aspect{$term}) {
            print STDERR "Unknown GO term ($term) for gene $gene\n";
            next;
        }
        ## Build a pseudo GAF line 
        my @out = ('DB', $gene, $gene, ' ', $term, 'PMID:foo', 'TAS', ' ', $aspect{$term},
                             $gene, ' ', 'protein', $taxon, $date, 'DB', ' ', ' ');
        print join("\t", @out). "\n";
    }
}
