#!/bin/sh

UB="N"
GFILE=""

usage()
{
    echo "
Usage: qsEq $0 [options] [BAM]

Options:
    -g CHROMINFO    List of chromosomes for ordering
    -U              Use corrected UMI (UB tag)

" >&2    
    exit 1
}

ARGS=$(getopt -o "g:U" -- "$@")
if [ $? -ne 0 ]; then
    echo "Invalid command-line parameters" >&2
    usage
fi

eval set -- "$ARGS"

while [ $# -gt 0 ];
do
    case "$1" in
	-g) GFILE="$2"
	    shift 2
	    ;;
	-U) UB="Y"
	    shift
	    ;;
	--) shift
	    break
	    ;;
    esac
done

if [ -z "$1" ]; then
    usage
fi

SCRIPTDIR="../mod_oliver_latest_scripts"
BAM="$1"
BASE=$(basename ${BAM} \.bam)

## DEBUG mode, remove if statement for production
if [ "${UB}" == "Y" ]; then
    perl "$SCRIPTDIR/extract_BC_UMI.pl" -u UB ${BAM} | sort -k1,1 -S 40G --parallel=13 -T $PWD > ${BASE}_BCUMI.txt &
else
    perl "$SCRIPTDIR/extract_BC_UMI.pl" ${BAM} | sort -k1,1 -S 40G --parallel=13 -T $PWD > ${BASE}_BCUMI.txt &
fi

if [ -z "$GFILE" ]; then
    sh "$SCRIPTDIR/annotate.sh" ${BAM} &
else
    sh "$SCRIPTDIR/annotate.sh" ${BAM} ${GFILE} &
fi

wait

join -t "	" -j 1 ${BASE}_BCUMI.txt ${BASE}_fullAnnot.txt > ${BASE}_annotation_output.txt

echo "Done generating annotation output"

expandCols -i ${BASE}_annotation_output.txt -c 4 | sort -k2,2 -k3,3 -s -S 45G -T $PWD --parallel=13 | groupBy -g 2,3 -c 4,1 -o distinct | awk -F "	" -v OFS="	" '{print $3,$1,$2,$4}' | sort -k1,1 -T $PWD -S 45G --parallel=13 -s > "${BASE}_BCUMIannot.txt"

if [ $? -ne 0 ]; then
    echo "Error in sorting BC & UMI" >&2
    exit 1
else
    echo "Done sorting BC & UMI"
fi

perl $SCRIPTDIR/make_EqC_Bfh.pl ${BASE}_BCUMIannot.txt

if [ $? -ne 0 ]; then
    echo "Error in EqC generating step" >&2
    exit 1
else
    echo "Done EqC generation"
fi
