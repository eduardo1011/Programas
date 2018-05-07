#!/usr/bin/perl

######################################################
#
#  new improved script to generate GO.terms_and_ids
#  using the gene_ontology.obo file
#
#  8-15 March 2004
######################################################

use strict;

my $outfile = "";

# get rid of old temporary file if there is one

if (-e "Terms_IDs")  { unlink "Terms_IDs" }

# open up new temporary holding file and actual file

open (TEMP, ">>Terms_IDs") or die "can't open temporary holding file \n";

# read in path to OBO format file
#  (for my use, provided as argument passed by cronjob)
#   value is /Users/midori/go/ontology/gene_ontology.obo

open (ONTFILE, "$ARGV[0]") or die "\n";

# now process the file contents ...
# want to bypass first few lines, which are essentially header info
# read each record into an array, then go thru the array to find id, name, namespace
#  (cos they're not necessarily in predictable positions)

while (<ONTFILE>){

    $/ = "\n[T";       # $/ is input record separator;
                       #  partial so both [Term] and [Typedef] will separate

    my $record = $_;
    my @contents = split(/\n/, $record);  # read the info into an array, splitting on newlines

    my $id = '';         # variables for what gets printed to temp file
    my $termname = '';   #   i.e. the bits i want out of the file
    my $namespace = '';

    foreach my $line (@contents) { # go thru the array
    chomp($line);
    if ($line =~ "^id: GO:") {     # find the ID line
	($id = $line) =~ s/id: //;
        }
    elsif ($line =~ "name: ") {    # find the term name
        ($termname = $line) =~ s/name: //;
        }
    elsif ($line =~ "namespace: ") { # ... and process/function/component
        if ($line =~ "function") { $namespace = "F"; }
        elsif ($line =~ "process") { $namespace = "P"; }
        elsif ($line =~ "component") { $namespace = "C"; }
        }
    }
    if ($id =~ "GO:" ) {                # don't print the junk from the header
        unless ($id =~ "GO:0003673") {  # don't include GO:0003673 Gene_Ontology
	    print TEMP "$id\t$termname\t$namespace\n";  ## original
	   #print TEMP "$id\t$termname\n";  ### change
        }
    }
}


# now we have the file of terms and primary IDs

# print header (lots of print commands just for legibility)

# get rid of old output file

if (-e $outfile)  { unlink $outfile }

open (OUT, ">>$outfile") or die " $outfile \n";
print OUT "! version: \$Revision: 1.3 $\n";
print OUT "! date: \$Date: 2008/01/26 00:34:29 $\n";
print OUT "!\n";
print OUT "! GO IDs and text strings (primary only)\n";
print OUT "! GO:0000000 [tab] text string	[tab] F|P|C  \n";
print OUT "! where F = molecular function, P = biological process, C = cellular component\n";
print OUT "!\n";

# sort the terms & IDs and plunk into the real file

system ("sort Terms_IDs >> $outfile");

# now the tricky part : cvs commit!!
#  must change to ssh eventually - can do so once keys are set up
#  mar 12 note: keys ok; prob have to remove passphrase - see man ssh-keygen

# system ("cvs -d :pserver:midori\@gocvs.geneontology.org:/share/go/cvs ci -m 'regular update' go/doc/GO.terms_and_ids");

# system ("cvs -d :ext:midori\@ext.geneontology.org:/share/go/cvs commit -m 'regular update' go/doc/GO.terms_and_ids");

__END__
