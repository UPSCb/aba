#!/bin/bash -l

usage(){
echo >&2 \
"
	Usage: $0 <genome> <host>

	Argument:
		genome: one of 'Potra01', 'Potrs01', 'Potri03', 'Potrx01'

  Note:
    the default node is watson. If you want to load on another computer provide the <host> argument
        Note:
             You need to set the UPSCb env. variable to your UPSCb git checkout directory.
"
	exit 1
}

## check we are set
if [ -z $UPSCb ]; then
    echo "You need to set the UPSCb environment variable to point to your UPSCb Git checkout"
    usage
fi

case $1 in
    Potra01)
	genome=/mnt/picea/storage/reference/Populus-tremula/v1.0/indices/STAR/Potra01
    ;;
    Potrs01)
	genome=/mnt/picea/storage/reference/Populus-tremuloides/v1.0/indices/STAR/Potrs01
	;;
    Potri03)
	genome=/mnt/picea/storage/reference/Populus-trichocarpa/v3.0/indices/STAR/Potri03
	;;
	  Potrx01)
	genome=/mnt/picea/storage/reference/Populus-tremula_X_Populus-tremuloides/v1.0/indices/STAR/Potrx01
	;;
    *)
	echo "ERROR: unknown genome"
	usage;;
esac

OPTION=watson
if [ ! -z $2 ]; then
	OPTION=$2
fi

## submit
sbatch -o ${genome}-${OPTION}-load.out -e ${genome}-${OPTION}-load.err -w $OPTION --mem=100G --mail-user="pal.miskolczi@slu.se" $UPSCb/pipeline/runSTARGenomeLoad.sh $genome

