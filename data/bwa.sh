if [ $# -ne 6 ];
then
	echo "sh $0
<ref.fa>
<r1,r2 or se>
<outdir>
<prefix>
<cpu for bwa>
<sam or bam for output>
"
	exit 1
fi
source /home/myth/miniconda3/install/etc/profile.d/conda.sh
conda activate  offtarget
set -e
ref=$1
qs=`echo $2|sed 's/,/ /g'`
outdir=$3
prefix=$4
cpu=$5
format=$6

mkdir -p $outdir
if [ ! -f "$ref.bwt" ];
then
	echo "start build index for $ref"
	bwa index $ref
fi

echo "get ref $ref"
echo "get query $qs"
out="$outdir/$prefix.$format"
if [ "$format" == "sam" ];
then	
	echo start: bwa mem $ref $qs -t $cpu -o $out
	bwa mem $ref $qs -t $cpu -o $out
elif [ "$format" == "bam" ];
then
	echo "start: bwa mem $ref $qs -t $cpu |samtools view -bS - >$out"
	bwa mem $ref $qs -t $cpu |samtools view -bS - >$out
else
	echo "error: not support output $format, only sam or bam"
	exit 1
fi
echo "output $out"
