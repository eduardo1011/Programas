#!/usr/bin/perl -w
#
#########################################################################################
#########################################################################################
#######CON ESTE SELLO IDENTIFICO QUE ES EL MODIFICADO Y EL QUE GENERA GRAFICOS###########
#########################################################################################
#########################################################################################
#
#    phobius.pl 1.01
#     A hidden Markov Model capable of predicting both Transmembrane Topology
#     and Signal peptides
#
#     (c) 2004 Lukas Kall, Anders Krogh, Erik Sonnhammer.
#
#     Please cite:
#
#     Lukas Käll, Anders Krogh and Erik L. L. Sonnhammer.
#     A Combined Transmembrane Topology and Signal Peptide Prediction Method.
#     Journal of Molecular Biology, 338(5):1027-1036, May 2004.
#

use IPC::Open2;
use Fcntl qw(:DEFAULT :flock);
use POSIX;

use FindBin;
#use lib $FindBin::RealBin;
#use seq_tools;

my $PHOBIUS_DIR =      "$FindBin::RealBin";
my $DECODEANHMM =      "$PHOBIUS_DIR/decodeanhmm";
my $PHOBIUS_OPT_FILE = "$PHOBIUS_DIR/phobius.options";
my $MODELFILE =        "$PHOBIUS_DIR/phobius.model";
my $PHOBIUS =          "$DECODEANHMM -f $PHOBIUS_OPT_FILE $PHOBIUS_DIR/phobius.model";
my $PHOBIUS_PLP =      "$DECODEANHMM -plp $PHOBIUS_DIR/phobius.model";
my $PHOBIUS_VER =      "\nPhobius ver 1.01\n\n";  #(c) 2004 Lukas Kall, Anders Krogh, Erik Sonnhammer\n\n";
my $GNUPLOT =          "gnuplot";
my $GNUPLOT_TERMINAL = "postscript eps color solid";


###### debido a que presentó errores para generar plots en formato .png, lo he modificado de tal forma que genere los gráficos sin problema
###### esto por si se desea adquirir el gráfico de alguna proteína de interés.
###### sustituí png por eps.
################################################################################################

sub tempname {
    my $suffix = $_[0];
    my $dir = $_[1];
    foreach $prefix ("aaaa".."zzzz") {
        $name = "$dir/$prefix$suffix";
        if (! (-e $name)) {
            open(TMP,">$name") || die; # Touch it to reserve it.
            close(TMP);
            return "$name";
        }
    }
}

# Print an entry with a label
sub print_entry {
    if (defined($_[1])) { $fh = $_[1]; }
    else {$fh = \*STDOUT;}

    $linelen = 70;
    @entry = @{$_[0]};
    print $fh ">" if ( $entry[0]!~/^>/ );
    print $fh "$entry[0]\n";
    $movein = "";
    if (defined($entry[2])) {
	$movein = '  ';
	$labpref = '# ' 
    }
    if ($#entry>2) {
	$movein = '   ';
	$labpref = '#  ';
    }
    $l = length($entry[1]);
    for ($i=0; $i<$l; $i += $linelen) {
	print $fh $movein . substr($entry[1],$i,$linelen) . "\n";
	if (defined($entry[2])) {
	    print $fh $labpref . substr($entry[2],$i,$linelen) . "\n";
	}
	for ($k=3; $k<$#entry; $k += 2) {
	    if ($entry[$k+1]) { 
		if (defined($entry[$k])) { print $fh "?$entry[$k] ";}
		else { print $fh '?  '; }
		print $fh substr($entry[$k+1],$i,$linelen) . "\n";
	    }
	}
    }
}


# This sub reads an entry in fasta-like format
# It returns a list containing the id, the sequence,the correct labels (#)
# and then alternating labels and strings.
$lastline = "";
sub read_entry {
    if (defined($_[0])) { $fh = $_[0]; }
    else {$fh = \*STDIN;}
    @entry = ();
    %labels = ();

    $nr=0;
    $nl = 3;
    $entry[0] = "";
    $nid = 0;
    if ($lastline) {
	$entry[0]=$lastline;
	$lastline = "";
	$nid=1;
    }
    else {$nid=0;}
    $entry[1] = "";

    while (<$fh>) {
	chop $_;
	if ($_ =~ /^>/) {
	    if ($nid==1) {
		$lastline = $_;
		return @entry;
	    }
	    $entry[0] = $_;
	    ++$nid;
	    ++$nr;
	}
	elsif ($nid && $_ =~ /^%/) {
	    # Append comments to the id line
	    $entry[0] = $entry[0] . "\n" . $_;
	}
	elsif ($nid) {
	    if ($_ =~ /^#/) {
		if ( !defined($entry[2]) ) { $entry[2] = ""; }
		$_ =~ s/^#//;
		$k = 2;
	    }
	    elsif ($_ =~ /^\?/) {
		$lab = substr($_,1,1);
		if (defined($labels{$lab})) {
		    $k=$labels{$lab};
		}
		else {
		    $entry[$nl] = $lab;
		    $labels{$lab} = $nl+1;
		    $k = $nl+1;
		    $nl +=2;
		    $entry[$k] = "";
		    if ( !defined($entry[2]) ) { $entry[2] = ""; }
		}
		$_ =~ s/^\?.//;
	    }
	    else { $k = 1; }
	    $_ =~ s/ //g;
	    $entry[$k] .= $_;
	    ++$nr;
	}
    }
    if (!$nr) { @entry = () ;}
    return @entry;
}


sub read_predict {
  my @entry = @{$_[0]};
  my ($rawpred,$line,$predstr);
  my ($label,$start,$stop);
  my @preds=();
  my $reclen = 1000;
  my $no_splits = ceil(length($entry[4])/$reclen);

  $stop=0;
  $label='';
  foreach $line (unpack "a$reclen" x $no_splits ,$entry[4] ) { # Take one part a time (titin protection)
    while ($line =~ /((\w)\2*)/g) {
      if ($label eq $2) {
	$stop += length($1);	  
      } else {
	push @preds,[$label,$start,$stop] if $label;
	$label = $2;
	$start = $stop + 1;
	$stop = $start + length($1) - 1;
      }
    }
  }
  push @preds,[$label,$start,$stop] if $label;
  return @preds;
}

sub printLongN1{
  my $seqid = $_[1];
  my $predstr;
  my ($label,$start,$stop);
  my ($tm,$sp);
  my $hanger = 0;
  my $template="FT   %-8s %6i %6i       %s\n";
  my $arrowtempl="set arrow from %i,-0.02 to %i,-0.02 lt %1i lw 30\n";
  $arrows = "";
  
  my @preds=&read_predict($_[0]);
  $tm=0;
  $sp=0;
  $stop=0;
  foreach my $ref (@preds) {
    ($label,$start,$stop) = @{$ref};
    if ($label eq "n") {
      $predstr=$predstr . sprintf($template,("DOMAIN",$start,$stop,"N-REGION."));
    }
    elsif ($label eq "h") {
      $predstr=$predstr . sprintf($template,("DOMAIN",$start,$stop,"H-REGION."));
      $sp = 1;
    }
    elsif ($label eq "c") {
      $predstr=$predstr . sprintf($template,("DOMAIN",$start,$stop,"C-REGION."));
      $predstr=sprintf($template,("SIGNAL",1,$stop,"")) . $predstr;
      $arrows=sprintf($arrowtempl,(1,$stop,6));
      $arrows=$arrows.sprintf($arrowtempl,($stop,1,6));
    }
    elsif ($label eq "C") {
      $hanger=$start;
    }
    elsif ($label eq "O"  || $label eq "o") {
      if ($hanger>0) {
	$start=$hanger;
	$hanger=0;
      }
      $predstr=$predstr . sprintf($template,("DOMAIN",$start,$stop,"NON CYTOPLASMIC."));
      $arrows=$arrows.sprintf($arrowtempl,($start,$stop,3));
      $arrows=$arrows.sprintf($arrowtempl,($stop,$start,3));
    }
    elsif ($label eq "M" ) {
      $predstr=$predstr . sprintf($template,("TRANSMEM",$start,$stop,""));
      $arrows=$arrows.sprintf($arrowtempl,($start,$stop,1));
      $arrows=$arrows.sprintf($arrowtempl,($stop,$start,1));
      $tm = $tm + 1;
    }
    elsif ($label eq "i") {
      $predstr=$predstr . sprintf($template,("DOMAIN",$start,$stop,"CYTOPLASMIC."));
      $arrows=$arrows.sprintf($arrowtempl,($start,$stop,2));
      $arrows=$arrows.sprintf($arrowtempl,($stop,$start,2));
    }
  }
  $predstr="ID   $seqid\n$predstr\/\/\n";
  print $predstr;
}

sub printShortN1{
  my $seqid = $_[1];
  my ($waist,$line,$predstr,$rawpred);
  my ($label,$start,$stop);
  my ($tm,$sp);
  my @preds=&read_predict($_[0]);

  $tm=0;
  $sp=0;
  foreach my $ref (@preds) {
    ($label,$start,$stop) = @{$ref};
    if ($label eq "" || $label eq "" || $label eq "") { ## eliminé la etiqueta "n" para que no aparezca en los resultados, y también la etiqueta "o" e "i"
      $predstr=$predstr . $label;
    }
    elsif ($label eq "h") {
      $predstr=$predstr .""; ## sustituí "$start-$stop" por "" para que no me muestre las posiciones de H-REGIONS
      $sp ="Signal_peptide";
    }
    elsif ($label eq "c") {
      $predstr=$predstr ."1\t$stop"; ## agregué el número 1, y espacios antes y después de la diagonal
    }
    elsif ($label eq "C") {
      $predstr=$predstr ."";
    }
    elsif ($label eq "O" ) {
      $predstr=$predstr .""; #$predstr=$predstr . "o";                ============modificado
    }
    elsif ($label eq "M" ) {
      $predstr=$predstr .""; #$predstr=$predstr . " $start-$stop ";              ===========modificado
      $tm += 1;
    }
  }
  printf "%s\t%s\t%s\t%s\n",($seqid,$sp,$tm,$predstr); ## cambié la posición de $sp y $tm
  #printf "%s\t%s\t%s\t%s\n",($seqid,$sp,$tm,$predstr); ## cambié la posición de $sp y $tm
}


sub printPosterior{
    my @seq = @{$_[0]};
    my $seqid = $_[1];
    my $options = $_[2];
    my $seqlen = length($seq[1]);
    local (*PI,*PO,*FPLP);
    delete $seq[4];
    delete $seq[2];
#    my $pid = open2(\*PO,\*PI,"$PHOBIUS -plp $options");
    my $pid = open2(\*PO,\*PI,"$PHOBIUS -plp $options 2>/dev/null");
    &print_entry(\@seq,\*PI);
    close(PI);
    open (\*FPLP,">$nplp");
    $lnum = 0;
    while($line=<PO>) {
      if ($line=~/^#\s+i\s+o\s+O/) {
	$line = "#pos\taa\ti\to\tM\tS\n";
      }
      elsif ($line=~/^\w/) {
	@myfield = split(/\s+/,$line);
	if ($myfield[1] !~ /[0-9,.]+/) {splice(@myfield,1,1)};
	$lnum = $lnum + 1;
	$line = "$lnum\t$myfield[0]\t$myfield[1]\t" . ($myfield[2]+$myfield[3]) . "\t$myfield[4]\t" . ($myfield[5]+$myfield[6]+$myfield[7]+$myfield[8]) . "\n" ;
      }
      print FPLP $line;
    }
    close(PO);
    close(FPLP);
    waitpid($pid,0);
    open(FGP,">$ngp");
    print FGP $arrows."set title \"Phobius posterior probabilities for $seqid \"
set key below
set yrange [-0.04:1.0]
set ylabel \"Posterior label probability\"
set xrange [1:$seqlen]
set terminal $GNUPLOT_TERMINAL
set output \"$neps\"
plot \"$nplp\" using 1:5 title \"transmembrane\" with impulses lt 1 lw 2, \\
\"\" using 1:3 title \"cytoplasmic\" with line lt 2 lw 3, \\
\"\" using 1:4 title \"non cytoplasmic\" with line lt 3 lw 3, \\
\"\" using 1:6 title \"signal peptide\" with line lt 6 lw 3
exit";
    close(FGP);
    system("$GNUPLOT $ngp");
    print "<img src=\"$neps\" alt=\"Posterior label probability plot\"></img>" if ($WEB_SERVER);
    print "<p><font size=\"-1\">The data to generate the plot is found <a href=$nplp>here</a>,
	and the gnuplot script is <a href=$ngp>here</a></font>.\n" if ($WEB_SERVER);
}

# main routine
my $format;
# OPTION PARSING ##########################################
use Getopt::Long;

$opt_short = 0;
$opt_long = 1;
$opt_plp = "";
$opt_eps = "";
$opt_gp = "";
$opt_raw = "";
$opt_h = 0;
$opt_help = 0;

# Process options
$result = GetOptions ('short!','long!','raw!','plp=s','eps=s','gp=s','h!','help!');
die ("Error on command line") unless $result;

print STDERR $PHOBIUS_VER;
if ($#ARGV>0 || $opt_h || $opt_help ) {
    print "usage: phobius.pl [options] [infile]\
\
  infile            A fasta file with the query protein sequences.\
                    If not present input will be read from stdin.\
\
options:\
  -h, -help         Show this help message and exit\
  -short            Short output. One line per protein.\
  -long             Long output. This is the default.\
  -raw              decodeanhmm's raw output. This overrides other options.\
  -eps FILE         If given a posterior label probability plot will be\
                    returned in the file FILE in eps format.\
  -gp FILE          Will in combination with the eps option write the\
                    gnuplot file to produce the plot into the file FILE.\
  -plp FILE         If given a posterior label probability values will be \
                    returned in the file FILE in ascii format.\n\n";
    exit;
}
local *SOUT;

if ($#ARGV==0)
{
    my $fname = shift @ARGV;
    if ($opt_raw) {
	system("$PHOBIUS $fname 2> /dev/null");
	exit;
    }
    open(*SOUT,"$PHOBIUS $fname 2> /dev/null |") || die "Can't open fasta file $fname: $!";
}
else
{
    if ($opt_raw) {
	system("$PHOBIUS 2> /dev/null");
	exit;
    }
    open(*SOUT,"$PHOBIUS 2> /dev/null |") || die "Can't reopen stdin: $!;";
}
$format = ($opt_long?"long":$format);    
$format = ($opt_short?"short":$format);    
$format = ($opt_eps?"plp":$format);    
if ($opt_eps || $opt_plp) {
    $neps = $opt_eps || tempname(".eps",".");
    $nplp = $opt_plp || tempname(".plp",".");
    $ngp =  $opt_gp || tempname(".gnuplot",".");	
}
#
# Get going
#
my $entries_done=0;
printf "%s\t%s\t%s\t%s\n",("qacc","Domain","TM","Start\tEnd") if ( $format eq "short" ); 
##### printf "%s\t%s\t%s\t%s\n",("ID","SP","TM","PREDICTION") if ( $format eq "short" ); 
##### En la linea anterior: cambié la posición de SP y TM, también sustituí los espacios por tabuladores, para mejor manejo del resultado
while ( @seq = &read_entry(\*SOUT)) {
  my @wrds = split(" ",$seq[0]);
  my $sid = $wrds[0];
  $sid =~ s/>//;
  if ( $format eq "short" ) {
    &printShortN1(\@seq,$sid,"");
  }
  elsif ($format eq "plp"){
    &printLongN1(\@seq,$sid,"");
    &printPosterior(\@seq,$sid,"") if ($entries_done<=0);
  }
  else {
    &printLongN1(\@seq,$sid,"");
  }
  $entries_done++;
}
close(*SOUT);
die ("Could not read provided fasta sequence") if ($entries_done==0);
if ($format eq "plp") {
  system ("rm $neps") if (!$opt_eps);
  system ("rm $nplp") if (!$opt_plp);
  system ("rm $ngp") if (!$opt_gp);
}
