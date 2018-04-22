#!/usr/bin/perl -w
use strict;
use HTTP::Request::Common;
use LWP::UserAgent;
use Pod::Usage;
use Getopt::Long;
use Bio::SeqIO;

#run this script with --help to see the options

=pod

=head1 NAME

SMART_batch  -  submit sequences from a FASTA file to SMART

=head1 SYNOPSIS

<SMART_batch.pl> <options>

=head1 NAME

options

=head1 SYNOPSIS

<--inputFile> <--outputDirectory> <--includePfam> <--includeSignalP> <--includeRepeats> <--includeDISEMBL> <--includeSchnipsel>

=head1 DESCRIPTION

Use <SMART_batch.pl> to submit multiple protein sequences from a FASTA file into the SMART analysis queue. Results are saved into plain text files.

=head1 GENERAL OPTIONS


=over 4

=item <--help>

display this message

=item <--inputFile>

FASTA file with sequences to submit

=item <--outputDirectory>

Directory which will be used to store the results. Will be created if it doesn't exist. Defaults to 'SMART_results'.


=back

=head1 ANALYSIS OPTIONS 

=over 4

=item <--includePfam>

Include Pfam domains in the search. (http://pfam.sanger.ac.uk/)

=item <--includeSignalP>

Include signal peptide predictions. (http://www.cbs.dtu.dk/services/SignalP/)

=item <--includeRepeats>

Include internal repeat predictions. (http://www.well.ox.ac.uk/rmott/ARIADNE/)

=item <--includeDISEMBL>

Include predictions of internal protein disorder. (http://dis.embl.de/)

=item <--includeSchnipsel>

Include predictions of outlier homologues and homologues of known structures. (http://smart.embl.de/help/smart_glossary.shtml#outlier)

=back


=head1 SEE ALSO

 SMART Home page : http://smart.embl.de
 SMART FAQ       : http://smart.embl.de/help/FAQ.shtml

=head1 AUTHOR

 Ivica Letunic <ivica@letunic.com>
 Contact me if you have any questions or comments.

=cut

my $submit_url = "http://smart.embl.de/smart/show_motifs.pl";
my $job_status_url = "http://smart.embl.de/smart/job_status.pl";
my ($show_help, $input_file, $output_directory, $do_pfam, $do_signalp, $do_rep, $do_disembl, $do_schnipsel);
my $op_r = GetOptions (
                       "help" => \$show_help,
                       "inputFile=s"   => \$input_file,
                       "outputDirectory=s"   => \$output_directory,
                       "includePfam" => \$do_pfam,
                       "includeSignalP" => \$do_signalp,
                       "includeRepeats" => \$do_rep,
                       "includeDISEMBL" => \$do_disembl,
                       "includeSchnipsel" => \$do_schnipsel,
                      ); 

unless ($input_file) { $show_help = 1; }

pod2usage(VERBOSE => 2) if ( $show_help );

my $ua  = LWP::UserAgent->new();
$ua->agent("SMARTbatch1.0");


print "\nSMART batch analysis\n======================\n";

unless (defined $output_directory) { $output_directory = 'SMART_results'; }
unless (-d $output_directory) { mkdir $output_directory; }
unless (-e $input_file) { print STDERR "Input file does not exist."; exit;}

my $io = new Bio::SeqIO(-format=> 'fasta', -file=> $input_file);

#process sequences one by one. ALWAYS wait for the results before submitting the next sequence.

while (my $seq = $io->next_seq) {
  my $seq_id = $seq->display_id;
  my $output_file = $output_directory . "/" . $seq_id . ".txt";
  if (-e $output_file) {
    my @s = stat($output_file);
    if ($s[7] == 0) {
      print "Removing empty results file $output_file.\n";
      unlink $output_file;
    } else {
      print "Skipping sequence $seq_id because the results file already exists.\n";
      next;
    }
  }
  print "Submitting sequence $seq_id...\n";
  #prepare the basic POST data
  my %post_content;
  $post_content{'SEQUENCE'} = $seq->seq;
  $post_content{'TEXTONLY'} = 1;
  if ($do_pfam) { $post_content{'DO_PFAM'} = 'DO_PFAM'; }
  if ($do_signalp) { $post_content{'INCLUDE_SIGNALP'} = 'INCLUDE_SIGNALP'; } 
  if ($do_rep) { $post_content{'DO_PROSPERO'} = 'DO_PROSPERO'; } 
  if ($do_disembl) { $post_content{'DO_DISEMBL'} = 'DO_DISEMBL'; } 
  if ($do_schnipsel) { $post_content{'INCLUDE_BLAST'} = 'INCLUDE_BLAST'; } 
  my $req = POST $submit_url, Content_Type => 'form-data', Content => [ %post_content ];
  my $response = $ua->request($req);
  if ($response->is_success()) {
    my @res = split(/\n/, $response->content);
    #check if we got the results directly (precomputed results)
    shift @res if ($res[1] =~ /^(--|[A-Z])/); ## he modificado este comando (original:   /^--\ SMART\ RESULT/   ) 
    if ($res[0] =~ /^(--|[A-Z])/) {  ## he modificado este comando (original:   /^--\ SMART\ RESULT/   )
      open (OUT, ">$output_file") or die "Cannot write to $output_file";
      print OUT $response->content;
      close OUT;
      print "Results saved to '$output_file'\n";
    } else {
      #we're in the queue, or there was an error
      my $job_id;
      for (my $i = 0; $i <= $#res; $i++) {
        if ($res[$i] =~ /job_status\.pl\?jobid=(\d+.+?)'/) {
          $job_id = $1;
          last;
        }
      }
      unless (length $job_id) {
        #there is no job ID, so an error occured
        my $error_file = "$output_directory/$seq_id\_SMART_error.html";
        open (ERR, ">$error_file") or die "Cannot write to $error_file";
        print ERR $response->content;
        close ERR;
        print "SMART returned an error page, which was saved into '$error_file'.\nPlease check the file for details. Aborting further submissions.\n";
        exit;
      } else {
        #we have a jobID, check every 10 seconds until we get the results
        print "Job entered the queue with ID $job_id. Waiting for results.\n";
        my $job_status_req = GET "$job_status_url?jobid=$job_id";
        sleep 5;
        while (1) {
          my $job_status_response = $ua->request($job_status_req);
          if ($job_status_response->is_success) {
            #check if we got the results
            my @job_status_res = split(/\n/, $job_status_response->content);
	     shift @job_status_res if ($job_status_res[1] =~ /^(--|[A-Z])/);  ## he modificado este comando (original:   /^--\ SMART\ RESULT/   )
            if ($job_status_res[0] =~ /^(--|[A-Z])/m) {   ## he modificado este comando (original:   /^--\ SMART\ RESULT/   )
              open (OUT, ">$output_file") or die "Cannot write to $output_file";
              print OUT $job_status_response->content;
              close OUT;
              print "Results saved to '$output_file'\n";
              last;
            } else {
              #still in queue
              sleep 10;
            }
          } else {
            print "SMART returned a web server error. Full message follows:\n\n";
            print $response->as_string;
            die;
          }
        }
      }
    }
  } else {
    print "SMART returned a web server error. Full message follows:\n\n";
    print $response->as_string;
    die;
  }
  #be nice to other users
  sleep 5;
}
