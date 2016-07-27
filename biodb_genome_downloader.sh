#!/bin/bash

## original script by Jorge Boucas (@jorgeboucas / GitHub)
## http://mpg-age-bioinformatics.github.io/

# @Author: Tobias Jakobi <tjakobi>
# @Date:   Tuesday, July 26, 2016 2:00 PM
# @Email:  tobias.jakobi@med.uni-heidelberg.de
# @Project: University Hospital Heidelberg, Section of Bioinformatics and Systems Cardiology
# @Last modified by:   tjakobi
# @Last modified time: Wednesday, July 27, 2016 8:00 PM
# @License: CC BY-NC-SA

## updated for HD cluster by Tobias / 14.04.2016

# This script builds the references and annotations folder for each selected organism and selected release.
# eg. caenorhabditis_elegans release: WBcell235.78
# Level1 folder = caenorhabditis_elegans
# Level2 folder = WBcell235_78
# reference genome (fasta and index) = WBcel235.dna.toplevel
# reference GTF = WBcel235.78.gtf
# reference GTF index = caenorhabditis_elegans/WBcell235_78/GTF_index/
# cuffcompare fixed GTF = cuffcmp_GTF.WBcel235.78.gtf
# other cuffcompare output files = cuffcmp.results.tar.bz2
# cuffcompare fixed GTF index = caenorhabditis_elegans/WBcell235_78/GTF_cuffcmp_GTF_index/
# chromosomes fasta = caenorhabditis_elegans/WBcell235_78/fasta

# example tophat options:
# --transcriptome-index $references_directory/caenorhabditis_elegans/WBcell235_78/GTF_cuffcmp_GTF_index
# for the reference indexed genome use $references_directory/caenorhabditis_elegans/WBcell235_78/bowtie2/WBcel235.dna.toplevel

# example cufflinks options:
# -g $references_directory/caenorhabditis_elegans/WBcell235_78/cuffcmp_GTF.WBcel235.78.gtf

# example cuffcompare options:
# -s $references_directory/caenorhabditis_elegans/WBcell235_78/chromosomes


## check if organsim directory already exists, if not create it
organism=$1

if [ -d "$organism" ]; then
  printf "Directory $organism already exists, skipping creation\n"
else
  mkdir $organism
fi

cd $organism

## Choice: either we get the latest release (default) or the user specifies
## the release

release_force=$2

if [ $release_force ] ; then

  release="release-$2"

  printf "ENSEMBL release $release manually selected\n"

  # check the latest release
  gtf_import=$(curl -s -l ftp://ftp.ensembl.org/pub/$release/gtf/$organism/ | grep gtf.gz)
  releaseA=$(echo $gtf_import | cut -f2 -d.)
  releaseB=$(echo $gtf_import | cut -f3 -d.)
  sep=_
  release_code=$releaseA$sep$releaseB
  printf "Genome ID is $release_code\n"

else
  printf "Get genome releases...\n"

  # check the latest release
  gtf_import=$(curl -s -l ftp://ftp.ensembl.org/pub/current_gtf/$organism/ | grep gtf.gz)
  releaseA=$(echo $gtf_import | cut -f2 -d.)
  releaseB=$(echo $gtf_import | cut -f3 -d.)
  release="release-$releaseB"
  sep=_
  release_code=$releaseA$sep$releaseB

  printf "Genome ID is $release_code\n"
fi

mkdir $release_code
cd $release_code



printf "Downloading data...\n"

wget -nc ftp://ftp.ensembl.org/pub/$release/fasta/$organism/dna/*.dna.*
wget -nc ftp://ftp.ensembl.org/pub/$release/gtf/$organism/*.gtf.gz

rm -v *abinitio.gtf.gz
rm -v *.chr.gtf.gz

printf "Uncompressing...\n"

parallel -v -j 40 srun gunzip {} ::: *.gz

printf "Housekeeping...\n"

namep1=$(echo $gtf_import | cut -f1 -d.)

for file in $(ls $namep1.*); do
nname=${file#"$namep1."}
nname=${nname%"$namep1."}
mv -v $file $nname; done


mkdir -v chromosomes
for file in $(ls *chromosome*); do
if [ ${file} != chromosomes: ]; then
chrp1=$(echo $file | cut -f4 -d.)
chrp2=$(echo $file | cut -f5 -d.)
sep=.
chr=$chrp1$sep$chrp2

mv -v ${file} chromosomes/${chr}; fi; done
mv -v *nonchromosomal* chromosomes/nonchromosomal.fa

original=$(ls *dna.toplevel.fa)
toplevel=${original#".fa"}
toplevel=${toplevel%".fa"}

# BWA index creation

printf "Generating BWA index\n"

module load bwa
mkdir bwa
cd bwa
full_path=$(pwd)
ln -s ../${original} ${original}
echo "#!/bin/bash
bwa index -a bwtsw -p ${full_path}/${original::(-3)} ${original}
which bwa" > bwa.sh
chmod 770 bwa.sh;
sbatch -J "bwa index $organism" -o bwa.log bwa.sh
cd ..


printf "Toplevel: $toplevel\n"

echo "Building bowtie2 index...\n"
# bowtie2 index
module load bowtie2
mkdir bowtie2
cd bowtie2
ln -vs ../${original} ${original}
srun -J "bowtie2 index $organism" -e bowtie2.err  -o bowtie2.log -c 2 --mem 40000 bowtie2-build-s $original $toplevel

cd ..
pwd
gtf=$(ls -C *.gtf)

# STAR index creation

printf "Generating STAR index\n"

module load star
mkdir star
cd star
full_path=$(pwd)
ln -s ../${original} ${original}
ln -s ../${gtf} ${gtf}
echo "#!/bin/bash
STAR --runMode genomeGenerate --genomeDir ${full_path} --genomeFastaFiles ${original} --runThreadN 40 --sjdbGTFfile ${gtf} --sjdbOverhang 100 --limitGenomeGenerateRAM 240000000000
which STAR" > star.sh;
chmod 770 star.sh;
sbatch -J "STAR index: $organism"  --cpus-per-task=40 --mem=256GB -o star.log star.sh
cd ..


# Fix GTF with cuffcompare

printf "Fixing GTF $gtf...\n"

module load cufflinks
which cuffcompare >> cuffcompare.log
srun -J "cuffcompare fix $organism" -c 20 --mem=40000 -o cuffcompare.log cuffcompare -V -CG -s chromosomes -r $gtf $gtf
mv -v cuffcmp.combined.gtf cuffcmp_GTF.$gtf
srun -o tar.log -J "cuffcompare taring" tar -cvf cuffcmp.results.tar cuffcmp.* --remove-files
srun -o pbzip2.log -J "cuffcompare bzipping" -c 20 --mem=4000 pbzip2 -p20 cuffcmp.results.tar

# Generate TOPHAT transcriptome indexes

printf "Indexing cuffcompare GTF\n"

module load tophat
mkdir cuffcmp_GTF_index
srun -J "TopHat2 $organism"  -e tophat2.err -o tophat2.log -c 40 --mem=40GB tophat2 -p 40 -G cuffcmp_GTF.$gtf --transcriptome-index cuffcmp_GTF_index bowtie2/$toplevel
which tophat2 >> tophat2.log
mv tophat2.log tophat2.err cuffcmp_GTF_index/
rm -r tophat_out


printf "Indexing GTF\n"

mkdir GTF_index
srun -J "TopHat2 $organism"  -e tophat2_gtf_index.err -o tophat2_gtf_index.log -c 40 --mem=40GB tophat2 -p 40 -G $gtf --transcriptome-index GTF_index bowtie2/$toplevel
which tophat2 >> tophat2_gtf_index.log
mv tophat2_gtf_index.log tophat2_gtf_index.err cuffcmp_GTF_index/
rm -r tophat_out


exit

else

cd $release

printf "
You already have this release in:
"
pwd
printf "
Press enter to list components.
"
read
ls
printf "
Do you wish to add more components? (y/n)
"
read answer

# If the user wants to add more components then they should be implemented into the script.

if [ $answer == 'n' ]; then exit;

elif [ $answer == 'y' ]; then

printf "
Please contact me: tobias.jakobi@med.uni-heidelberg.de | (+49) 6221 56 39126
";

else

printf "
Exiting.. you can only type y or n
"
exit

fi; fi

exit
