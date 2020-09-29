#!/bin/bash

helpMessage() {
	echo "$0 usage -i <trinity output directory> -h [show this message]"
	echo ""
	exit 1
}

# cmd line args
while getopts h:i: option
do
	case "${option}"
		in
		h) helpMessage;;
		i) inDir=${OPTARG};;
	esac
done
shift $((OPTIND - 1))

# loop
for x in ${inDir}/recursive_trinity.cmds.chunk.??
do
	# file basenames
	y=$(basename ${x}).completed
	#z=$(basename ${x})
	
	a=${inDir}/${y}

	cmds=$(wc -l ${x})
	num_cmds=$(sed 's/\s.*$//' <<< ${cmds})
	
	comp=$(wc -l ${a})
	num_comp=$(sed 's/\s.*$//' <<< ${comp})
	
	incomplete=$((num_cmds - num_comp))
	#echo $cmds
	#echo $comp
	echo "${x}  ...  ${incomplete}/${num_cmds} incomplete assemblies."
done
