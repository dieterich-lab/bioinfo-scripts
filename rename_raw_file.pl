#!/usr/bin/env perl

# @Author: Tobias Jakobi <tjakobi>
# @Date:   03/05/2016
# @Email:  tobias.jakobi@med.uni-heidelberg.de
# @Project: University Hospital Heidelberg, Section of Bioinformatics and Systems Cardiology
# @Last modified by:   tjakobi
# @Last modified time: Wednesday, May 4, 2016 10:29 AM
# @License: CC BY-NC-SA

use Getopt::Std;
use Cwd;
use strict;
use warnings;

### print usage information
sub usage() {
    print  <<EOF;

  This script links/renames raw sequencing files received from sequencing centers
  to a more readable and informative scheme:

  e.g.:
  1557_AC_run251_ATGAGCAT_L006_R1_001_f.fastq.gz
  -> Human_brain_R2_run251_002.fastq.gz

  usage: $0
            -o [list of original sample identifier, e.g. AC]
            -n [list of new IDs, e.g. HomeSapiens_WholeBrain_RNaseR_minus]
            -r rename files (instead of linking)
            -e [file extension, e.g. gz or bz2]
            -s [in-name seperator e.g. _ or -]
            -t [target directory for links e.g. /home/juser/prj_1/links/]
            -h this help

Contact: tobias.jakobi\@med.Uni-Heidelberg.DE
         https://github.com/tjakobi/bioinfo-scripts


EOF
    exit;
}


### set options
my %options = ();

getopts( "s:o:n:e:t:drlh", \%options ) or usage;

usage() if $options{h};

## check for parameters

die("please specify list of original sample identifier, e.g. AC,AD,AF\n")
if ( !defined $options{o} );

die("please specify list of new IDs, e.g. Control,Treatment,Treatment2\n")
if ( !defined $options{n} );

die("-o and -n list have to have the same number of entries\n")
if ( split(',', $options{o}) != split(',', $options{n}) );

die("please specify the file extension\n")
if ( !defined $options{e} );

die("please specify the in-name seperator\n")
if ( !defined $options{s} );

die("please specify the target directory\n")
if ( !defined $options{t} );

my @old_tags = split(',', $options{o});
my @new_tags = split(',', $options{n});


### does the actual work
sub do_rename {

  if ($_[0] =~ /.*$options{s}$_[1]$options{s}run(\d{1,}).*R([1,2]).*(\d{3})/) {

  my $new_name = $_[2];
     $new_name .= $options{s};
     $new_name .= "R$2";
     $new_name .= $options{s};
     $new_name .= "run$1";
     $new_name .= $options{s};
     $new_name .= "$3";
     $new_name .= ".";
     $new_name .= $options{e};

  return $new_name;
 }

}

### scan for files

opendir (DIR, '.') or die "Cannot list current directory";

while (my $file = readdir(DIR)) {

    # ignore files beginning with a period
    next if ($file =~ m/^\./);

    # ignore files with wrong extension
    next if ($file !~ m/$options{e}$/);

    ## look for suitable files
    for(my $i=0;$i<=$#old_tags;$i++){

      # ignore files with wrong extension
      if ($file =~ m/.*$options{s}$old_tags[$i]$options{s}.*/){

        # launch processing
        my $name = do_rename($file,$old_tags[$i],$new_tags[$i]);

        if ($name){

          my $here = getcwd();
          my $cmd = "ln -v -s $here/$file $options{t}/$name";

          if ($options{d}){
            print "$cmd\n";
            #print "$file => $name \n";
          } else {
            system($cmd);
          }
        }
      }

    }
}
