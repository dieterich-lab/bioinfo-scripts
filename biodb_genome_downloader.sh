#!/bin/bash

## original script by Jorge Boucas (@jorgeboucas / GitHub)
## http://mpg-age-bioinformatics.github.io/

# @Author: Tobias Jakobi <tjakobi>
# @Email:  tobias.jakobi@med.uni-heidelberg.de
# @Project: University Hospital Heidelberg, Section of Bioinformatics and Systems Cardiology
# @License: CC BY-NC-SA

## updated for HD cluster by Tobias / 14.04.2016

# This script builds the references and annotations folder for each selected organism and selected release.
# eg. caenorhabditis_elegans release: WBcell235.78
# Level1 folder = caenorhabditis_elegans
# Level2 folder = WBcell235_78
# reference genome (fasta and index) = WBcel235.dna.toplevel
# reference GTF = WBcel235.78.gtf
# reference GTF index = caenorhabditis_elegans/WBcell235_78/gtf_index/
# cuffcompare fixed GTF = cuffcmp_GTF.WBcel235.78.gtf
# other cuffcompare output files = cuffcmp.results.tar.bz2
# cuffcompare fixed GTF index = caenorhabditis_elegans/WBcell235_78/GTF_cuffcmp_gtf_index/
# chromosomes fasta = caenorhabditis_elegans/WBcell235_78/fasta

# example tophat options:
# --transcriptome-index $references_directory/caenorhabditis_elegans/WBcell235_78/GTF_cuffcmp_gtf_index
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

latest_ensembl_release=$(curl -s ftp://ftp.ensembl.org/pub/current_README | grep Release | cut -f3 -d " ")

## Choice: either we get the latest release (default) or the user specifies
## the release
release=0
release_force=$2
echo $2
alias sbatch='sbatch --mail-type=END,FAIL,TIME_LIMIT_80  --mail-user=tobias.jakobi@med.uni-heidelberg.de'
alias srun='srun --mail-type=END,FAIL,TIME_LIMIT_80  --mail-user=tobias.jakobi@med.uni-heidelberg.de'

if [ -n "$release_force" ] && [ -z "$3" ] && [ "$release_force" == "pre" ] ; then # get the pre release

  release=$2
  release_str="$release"

  printf "ENSEMBL release $release manually selected ($latest_ensembl_release is current)\n"

  # check the latest release
  gtf_import=$(curl -s -l ftp://ftp.ensembl.org/pub/pre/gtf/$organism/ | grep gtf)
  releaseA=$(echo $gtf_import | cut -f2 -d.)
  releaseB=$(echo $gtf_import | cut -f3 -d.)
  sep=_
  release_code=$releaseA$sep$releaseB
  if [ "${release_code/$release}" != "$release_code" ]
  then
    echo "Genome release includes ENSEMBL build";
  else
    echo "Genome release does not include ENSEMBL build, adding";
    release_code=$release_code"_"$release
  fi

  release=$release
  printf "Genome ID is $release_code\n"
  unset $release_force

elif [ -n "$release_force" ]  && [ -z "$3" ] && [ "$release_force" != "pre" ] ; then # get the selected release

  release=$2
  release_str="release-$release"

  printf "ENSEMBL release $release manually selected ($latest_ensembl_release is current)\n"

  # check the latest release
  gtf_import=$(curl -s -l ftp://ftp.ensembl.org/pub/release-$release/gtf/$organism/ | grep gtf.gz)
  releaseA=$(echo $gtf_import | cut -f2 -d.)
  releaseB=$(echo $gtf_import | cut -f3 -d.)
  sep=_
  release_code=$releaseA$sep$releaseB
  if [ "${release_code/$release}" != "$release_code" ]
  then
    echo "Genome release includes ENSEMBL build";
  else
    echo "Genome release does not include ENSEMBL build, adding";
    release_code=$release_code"_"$release
  fi

  release=release-$release
  printf "Genome ID is $release_code\n"

elif [ -z "$3" ] ; then # just get the last ensembl release
  printf "Get genome releases...\n"

  printf "Latest ENSEMBL release is: $latest_ensembl_release\n"
  gtf_import=$(curl -s -l ftp://ftp.ensembl.org/pub/release-$latest_ensembl_release/gtf/$organism/ | grep gtf.gz)
  releaseA=$(echo $gtf_import | cut -f2 -d.)
  releaseB=$(echo $gtf_import | cut -f3 -d.)
  release=release-$latest_ensembl_release
  sep=_
  release_code=$releaseA$sep$releaseB
  if [ "${release_code/$latest_ensembl_release}" != "$release_code" ]
  then
    echo "Genome release includes ENSEMBL build";
  else
    echo "Genome release does not include ENSEMBL build, adding";
    release_code=$release_code"_"$latest_ensembl_release
  fi

  printf "Genome ID is $release_code\n"
elif [ -n "$3" ] ;then # no ensembl, take the supplied fasta and gtf file

	printf "Manual node selected.\n"
	printf "Organism name: $1\n"
	printf "Release name: $2\n"
	printf "FASTA file: $3\n"
	printf "GTF file: $4\n"

	organism=$1
	release_code=$2
	orginal=$4
	
	if [ -d "$release_code" ]; then
	  printf "Directory $release_code already exists, exiting program\n"
	  exit
	fi

	mkdir $release_code
	cd $release_code
	cp $4 . -av
	cp $3 . -av

        original=$(ls *.fa)
        toplevel=${original#".fa"}
        toplevel=${toplevel%".fa"}
fi







if  [ -z "$3" ] ; then

	if [ -d "$release_code" ]; then
	  printf "Directory $release_code already exists, exiting program\n"
	  exit
	fi

	printf "Downloading data...\n"
	
	mkdir $release_code
	cd $release_code



	if [ "$release_force" = "pre" ] ; then
		# we have to test if we have .dna. files, otherwise we take the _RM files
		files=$(curl -s -l ftp://ftp.ensembl.org/pub/pre/fasta/dna/$organism/ | grep ".dna\.")
		
		if [ -n "$files" ]
		then 
			wget -nc -nv ftp://ftp.ensembl.org/pub/pre/fasta/dna/$organism/*.dna.toplevel.fa.gz
		else	
			wget -nc -nv ftp://ftp.ensembl.org/pub/pre/fasta/dna/$organism/*.dna_rm.toplevel.fa.gz
		fi
	else
		wget -nc -nv ftp://ftp.ensembl.org/pub/$release/fasta/$organism/dna/*.dna.toplevel.fa.gz
		wget -nc -nv ftp://ftp.ensembl.org/pub/$release/fasta/$organism/dna/*.dna.primary_assembly.fa.gz
                wget -nc -nv ftp://ftp.ensembl.org/pub/$release/fasta/$organism/cdna/*.cdna.all.fa.gz
	fi

	wget -nc -nv ftp://ftp.ensembl.org/pub/$release/gtf/$organism/*.gtf.gz

	rm -v *abinitio.gtf.gz
	rm -v *.chr.gtf.gz
	rm -v *.patch_hapl_scaff.gtf.gz	

	printf "Uncompressing...\n"

	parallel -v -j 40 srun gunzip {} ::: *.gz

	printf "Housekeeping...\n"

	namep1=$(echo $gtf_import | cut -f1 -d.)

	for file in $(ls $namep1.*); do
	nname=${file#"$namep1."}
	nname=${nname%"$namep1."}
	mv -v $file $nname; done


	#mkdir -v chromosomes
	#for file in $(ls *chromosome*); do
	#if [ ${file} != chromosomes: ]; then
	#chrp1=$(echo $file | cut -f4 -d.)
	#chrp2=$(echo $file | cut -f5 -d.)
	#sep=.

	#chr=$(echo $file | egrep -o '[[:alnum:]]*.fa')

	#mv -v ${file} chromosomes/${chr}; fi; done
	#mv -v *nonchromosomal* chromosomes/nonchromosomal.fa

	original=$(ls *dna.primary_assembly.fa)
	# since 09/2016 we use primary assembly instead of toplevel
	#original=$(ls *dna.toplevel.fa)

	# but we still use top level as fallback if no assembly exists, i.e. pre builds
	if [ -z "$original" ] ; then
		original=$(ls *dna.toplevel.fa)
	fi


	mv -v ${original} ${release_code}.fa 

	original=$(ls ${release_code}.fa)

	toplevel=${original#".fa"}
	toplevel=${toplevel%".fa"}
fi 



printf "Fasta indexing\n"

module load samtools
echo "#!/bin/bash
samtools faidx ${original}
which samtools" > samtools.sh
chmod 770 samtools.sh;
sbatch -J "Fasta index: $organism / $release" -o samtools.log samtools.sh



printf "Running RepeatModeler & RepeatMasker\n"

module load repeatmasker
module load repeatmodeler

mkdir -v repeat_masker
cd repeat_masker
full_path=$(pwd)
ln -vs ../${original} ${original}
echo "#!/bin/bash
BuildDatabase -name $organism -engine ncbi ${original}
RepeatModeler -pa 10 -engine ncbi -database $organism
MODEL_PATH=\$(grep consensi.fa.classified repeat_masker.log | cut -f3 -d ' ')
RepeatMasker -s -u -ace -xm -html -norna --gff -pa 10 -lib \$MODEL_PATH/consensi.fa.classified ${original}  
# RepeatMasker -s -u -ace -xm -html -norna --gff -pa 10 -lib consensi.fa.classified ${original} 
# RepeatMasker -s -u -ace -xm -html -norna --gff -pa 10 -species "${organism/_/ }" ${original}
ln -vs repeat_masker/${original}.out.gff ../${original::(-3)}_repeatmasker.gtf
tar -cvf repeat_masker.tar \$MODEL_PATH
pbzip2 -v -p40 repeat_masker.tar
rm -v RM_* -rf
which RepeatMasker" > repeat_masker.sh
chmod 770 repeat_masker.sh;
sbatch -p long --cpus-per-task=40 --mem=250GB -J "RepeatModeler: $organism / $release" -o repeat_masker.log repeat_masker.sh
cd ..



module load blast
mkdir blast
cd blast
full_path=$(pwd)
ln -s ../${original} ${original}
echo "#!/bin/bash
makeblastdb -in ${original} -out ${original/\.fa/}  -dbtype nucl -logfile ../blast.log " > blast.sh
chmod 770 blast.sh;
sbatch -p general --cpus-per-task=2 --mem=10GB -J "makebelastdb: $organism / $release" -o /dev/null blast.sh
cd ..

module load hisat2
mkdir hisat2
cd hisat2
full_path=$(pwd)
ln -s ../${original} ${original}
echo "#!/bin/bash
hisat2-build $original $toplevel
which hisat2" > hisat2.sh
chmod 770 hisat2.sh;
sbatch -p general --cpus-per-task=2 --mem=50GB -J "HISAT2 index: $organism / $release" -o hisat2.log hisat2.sh
cd ..


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
sbatch -p general --cpus-per-task=2 --mem=50GB -J "BWA index: $organism / $release" -o bwa.log bwa.sh
cd ..


printf "Toplevel: $toplevel\n"

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
STAR --runMode genomeGenerate --genomeDir ${full_path} --genomeFastaFiles ${original} --runThreadN 40 --limitGenomeGenerateRAM 240000000000
which STAR" > star.sh;
chmod 770 star.sh;
sbatch -p long -J "STAR index: $organism / $release"  --cpus-per-task=40 --mem=250G -o star.log star.sh
cd ..


# salmon index creation

printf "Generating salmon index\n"

module load salmon
mkdir salmon
cd salmon
full_path=$(pwd)
ln -s ../*cdna*fa .
transcript=$(ls -C *.fa)
echo "#!/bin/bash
# k31 -> 75bp reads or longer
salmon index -k31 -t ${transcript} -i ${full_path} -p 40 
which salmon" > salmon.sh;
chmod 770 salmon.sh;
sbatch -p long -J "salmon index: $organism / $release"  --cpus-per-task=40 --mem=250G -o salmon.log salmon.sh
cd ..

# Fix GTF with cuffcompare

printf "Fixing GTF $gtf...\n"

module load cufflinks
which cuffcompare >> cuffcompare.log
srun -J "Cuffcompare: $organism / $release" -c 20 --mem=40G -o /dev/null cuffcompare -V -CG -s ${original} -r $gtf $gtf
mv -v cuffcmp.combined.gtf tmp.gtf
rm -v cuffcmp.*
mv -v tmp.gtf cuffcmp.$gtf
#srun -o tar.log -J "cuffcompare taring" tar -cvf cuffcmp.results.tar cuffcmp.* --remove-files
#srun -o pbzip2.log -J "cuffcompare bzipping" -c 20 --mem=4000 pbzip2 -p20 cuffcmp.results.tar


echo "Building bowtie2 index...\n"
# bowtie2 index
module load bowtie2
mkdir bowtie2
cd bowtie2
ln -vs ../${original} ${original}
srun -J "Bowtie2 index: $organism / $release" -e bowtie2.err  -o bowtie2.log -c 2 --mem 50G bowtie2-build-s $original $toplevel

cd ..
pwd

# Generate TOPHAT transcriptome indexes

printf "Indexing cuffcompare GTF\n"

module load tophat
mkdir cuffcmp_gtf_index
srun -J "TopHat2: $organism / $release"  -e tophat2.err -o tophat2.log -c 40 --mem=40GB tophat2 -p 40 -G cuffcmp.$gtf --transcriptome-index cuffcmp_gtf_index bowtie2/$toplevel
which tophat2 >> tophat2.log
mv tophat2.log tophat2.err cuffcmp_gtf_index/
rm -r tophat_out


printf "Indexing GTF\n"

mkdir gtf_index
srun -J "TopHat2 GTF index: $organism / $release"  -e tophat2_gtf_index.err -o tophat2_gtf_index.log -c 40 --mem=40GB tophat2 -p 40 -G $gtf --transcriptome-index gtf_index bowtie2/$toplevel
which tophat2 >> tophat2_gtf_index.log
mv tophat2_gtf_index.log tophat2_gtf_index.err cuffcmp_gtf_index/
rm -r tophat_out


