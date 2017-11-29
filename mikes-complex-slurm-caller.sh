#!/bin/bash
# Example "mapping/SID1120_NoIndex_HISATmapping/Aligned.out.bam"

### USAGE ######################################################################
read -r -d '' USAGE << ENDOFUSAGE

Configure the Makefile.snps for the specific files found in the mapping
directory (created previously in pipeline).

  -h | --help
        Print this help statement.

  -E | --example-parameter
        Put a switche that takes a parameter.
        (DEFAULT: "default_num_or_string")

  -e | --example-flag
        Put a switche that triggers a state.
        (DEFAULT: "default_num_or_string")

  -L | --logfile
        The path to the logfile here.

  -v | --verbose
        Verbose output.

ENDOFUSAGE
### USAGE ######################################################################

#--- SET GETOPTS --------------------------------------------------------------
OPTS=`getopt -o vheE: --longoptions mapping-dir,verbose,help,example-flag,example-parameter,logfile: -n 'parse-options' -- "$@"`
if [ $? != 0 ] ; then
    echo "Failed parsing options." >&2
    echo "$USAGE" >&2
    exit 1
fi
# echo "$OPTS"
eval set -- "$OPTS"

#--- SET DEFAULTS --------------------------------------------------------------
ROOT_DIR="./"
EXAMPLE_PARAM="default_num_or_string"
VERBOSE=false
LOGFILE=""
DRY_RUN=false

#--- RUN GETOPTS --------------------------------------------------------------
while true; do
  case "$1" in
    -h | --help ) 			  	echo "$USAGE" >&2;		shift ;;
    -v | --verbose ) 		  	VERBOSE=true; 			shift ;;
    -L | --logfile ) 			LOGFILE="$2"; 			shift; shift ;;
    -e | --example-flag ) 		FLAGEXAMPLE=true; 			shift ;;
    -d | --dry-run ) 		    DRY_RUN=true; 			shift ;;
    -E | --example-parameter ) 	VALUEEXAMPLE="$2"; 			shift; shift ;;

    -- ) shift; break ;;
    * ) break ;;
  esac
done

#--- PARAM MODS --------------------------------------------------------------
#LOGFILE should always come first
if [ ${#LOGFILE} -lt 1 ]; then # NOT SET
  LOGFILE="./MyLog.$(date +%y%m%d.%H%M%S).log"

elif [[ ! $LOGFILE =~ .*/.* ]]; then # SET, but NOT a Full path
  LOGFILE="${ROOT_DIR}/${LOGFILE}"

# else # No else, just go with whats set
fi
# And redirect for logfile
exec > >(tee -a $LOGFILE)


#--- FUNCTIONS --------------------------------------------------------------
function some_function_here() {
        # This is actually the return line
        echo $1 # $1 is what was input to the function
}

# MAIN() ---------------------------------------------------------------------
if [ ${VERBOSE} == true ]; then
  echo "ROOT_DIR        = ${ROOT_DIR}"
  echo "EXAMPLE_PARAM   = ${EXAMPLE_PARAM}"
  echo "VERBOSE         = ${VERBOSE}"
  echo "LOGFILE         = ${LOGFILE}"
fi

# RUN CODE TO BUILD SBATCH FILE'S COMMAND(S) ##################################

for line in some_function_here ${EXAMPLE_PARAM}; do
    command="sometool ${line} | grep -v nonsense > ${EXAMPLE_PARAM}.out"
    command2="sometool ${line} | grep -v nonsense > ${EXAMPLE_PARAM}.out"

# === sbatch file =========================================================
cat << EOF > $tmpfile
#!/bin/bash
#SBATCH --job-name=bam-calls-${patient_name}
#SBATCH --output=${slurmoutput}
#SBATCH --mem=10G
#SBATCH -c 2

# FASTA_REFERENCE using ... ${GENOME_REF_PATH}
# based on "workflow_DNAseq_PE_MPIZ_homoSapiens.pl"

module load sometool
echo "Running command: ${command}"
${command} && echo "OK" || "FAILED!"
echo "Running command: ${command2}"
${command2} && echo "OK" || "FAILED!"
EOF
# === END sbatch file =========================================================

# Run the temp script
  count=$(($count + 1))
  chmod 777 ${tmpfile}
  echo "input_file = ${line}"
  if [ ${DRY_RUN} == true ]; then
    echo "Dry run with slurm batch file..."
    echo "---${tmpfile}-----------------------"
    cat ${tmpfile}
    echo "------------------------------------"
  else
    echo -n "Adding ${tmpfile} to slurm as job #${count} ... "
    sbatch ${tmpfile} && echo "OK" || echo "FAILED!"
  fi
  # once started via slurm, file can be deleted
  echo -n "Removing ${tmpfile} ..."
  rm -rf ${tmpfile} && echo "OK" || echo "FAILED!"
### END Make and run the sbatch file

done
