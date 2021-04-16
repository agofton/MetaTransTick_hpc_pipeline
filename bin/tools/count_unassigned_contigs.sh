#!/bin/bash

helpMes() {
	"""
	$0 usage: -i <input tab separated lca.txt file from MEGAN6> -o <stats.txt> -h [show this message]
	 < > = required, [ ] = optional
	 Alexander W. Gofton 2021.
	 """
}

errorExit() {
	if [[ $? -ne 0 ]]; then
		echo $1; date; exit 1
	fi
}

# Command line arguments
while getopts hi:o: option
do
	case "${option}"
		in
		h) helpMessage;;
	  	i) input=${OPTARG};;
		o) stats_out=${OPTARG};;
	esac
done
shift $((OPTIND - 1))

# Count
total=$(cat ${input} | wc -l)

NCBI=$(grep -c NCBI$ ${input})
CO=$(grep -c cellular_organisms$ ${input})
NA=$(grep -c Not_assigned$ ${input})
tot_unass=$(expr $NCBI + $CO + $NA)

Bac=$(awk -F "\t" '{print $4}' ${input} | grep -c Bacteria)
Euk=$(awk -F "\t" '{print $4}' ${input} | grep -c Eukaryota)
Vir=$(awk -F "\t" '{print $3}' ${input} | grep -c Viruses)

echo ${input}
echo "Total contigs: ${total}"
echo "Unassigned contigs: ${tot_unass}"
echo "Bacterial contigs: ${Bac}"
echo "Eukaryotic contigs: ${Euk}"
echo "Viral contig: ${Vir}"
echo ""
echo "${tot_unass} ${Bac} ${Euk} ${Vir}"
echo ""

# write to file
echo ${input} > ${stats_out}
echo "Total contigs: ${total}" >> ${stats_out}
echo "Unassigned contigs: ${tot_unass}" >> ${stats_out}
echo "Bacterial contigs: ${Bac}" >> ${stats_out}
echo "Eukaryotic contigs: ${Euk}" >> ${stats_out}
echo "Viral contig: ${Vir}" >> ${stats_out}






