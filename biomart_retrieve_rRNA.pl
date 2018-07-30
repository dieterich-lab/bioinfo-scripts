#!/usr/bin/env perl

# An example script demonstrating the use of BioMart API.
# This perl API representation is only available for configuration versions >=  0.5
use strict;
use warnings;
use Getopt::Std;
use BioMart::Initializer;
use BioMart::Query;
use BioMart::QueryRunner;


my %opts;

getopts( 's:', \%opts );

sub usage {
    my $msg = shift;
    print "usage: $0 -s species_string \n\n";
    die "$msg\n";
}

usage("ERROR: No species string provided\n") unless ( $opts{'s'} );

my $full_name = $opts{'s'};
my @species = split('_',$opts{'s'});
my $file_name = substr($species[0], 0, 1)."".$species[1];
my $file_mt = $file_name.".mt-rRNA.fasta";
my $mt_name =  $file_name.".mt-rRNA";
my $file_rrna = $file_name.".rRNA.fasta";
my $rrna_name =  $file_name.".rRNA";

system("mkdir -v $full_name;");
#system("mkdir -v $full_name/rRNA");
#system("mkdir -v $full_name/mt-rRNA;");
#system("mkdir -v $full_name/mt-tRNA;");

my $confFile = "/biosw/biomart-perl/0.7/conf/martURLLocation.xml";

my $db_name = $file_name."_gene_ensembl";

# NB: change action to 'clean' if you wish to start a fresh configuration
# and to 'cached' if you want to skip configuration step on subsequent runs from the same registry

my $action='cached';
my $initializer = BioMart::Initializer->new('registryFile'=>$confFile, 'action'=>$action);
my $registry = $initializer->getRegistry;

print "Using DB $db_name\n";

my $query_mtRNA = BioMart::Query->new('registry'=>$registry,'virtualSchemaName'=>'default');

        $query_mtRNA->setDataset($db_name);
	$query_mtRNA->addFilter("biotype", ["Mt_rRNA"]);
        $query_mtRNA->addAttribute("ensembl_gene_id");
        $query_mtRNA->addAttribute("ensembl_transcript_id");
        $query_mtRNA->addAttribute("gene_exon");

$query_mtRNA->formatter("FASTA");

my $query_runner = BioMart::QueryRunner->new();

$query_runner->uniqueRowsOnly(1);
$query_runner->execute($query_mtRNA);

open( MT , ">", $full_name."/$file_name.mt-rRNA.fasta" ) || die "Can't open file 1\n";
$query_runner->printHeader(*MT);
$query_runner->printResults(*MT);
$query_runner->printFooter(*MT);
close(MT);


$query_mtRNA = BioMart::Query->new('registry'=>$registry,'virtualSchemaName'=>'default');

        $query_mtRNA->setDataset($db_name);
	$query_mtRNA->addFilter("biotype", ["Mt_tRNA"]);
        $query_mtRNA->addAttribute("ensembl_gene_id");
        $query_mtRNA->addAttribute("ensembl_transcript_id");
        $query_mtRNA->addAttribute("gene_exon");

$query_mtRNA->formatter("FASTA");

$query_runner = BioMart::QueryRunner->new();

$query_runner->uniqueRowsOnly(1);
$query_runner->execute($query_mtRNA);

open( MT , ">", $opts{'s'}."/$file_name.mt-tRNA.fasta" ) || die "Can't open file 2\n";
$query_runner->printHeader(*MT);
$query_runner->printResults(*MT);
$query_runner->printFooter(*MT);
close(MT);

my $query_rRNA = BioMart::Query->new('registry'=>$registry,'virtualSchemaName'=>'default');

	$query_rRNA->setDataset($db_name);
	$query_rRNA->addFilter("biotype", ["rRNA"]);
	$query_rRNA->addAttribute("ensembl_gene_id");
	$query_rRNA->addAttribute("ensembl_transcript_id");
	$query_rRNA->addAttribute("gene_exon");

$query_rRNA->formatter("FASTA");

$query_runner = BioMart::QueryRunner->new();

$query_runner->uniqueRowsOnly(1);
$query_runner->execute($query_rRNA);

open( RNA , ">", $opts{'s'}."/$file_name.rRNA.fasta" ) || die "Can't open file 3\n";
$query_runner->printHeader(*RNA);
$query_runner->printResults(*RNA);
$query_runner->printFooter(*RNA);
close(RNA);

#system("module load bowtie2; cd $full_name/rRNA_and_mt-rRNA/; bowtie2-build $file_mt $mt_name;");
#system("module load bowtie2; cd $full_name/rRNA/; bowtie2-build $file_rrna $rrna_name;");
