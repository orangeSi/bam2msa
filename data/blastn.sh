if [ $# -ne 5 ];
then
	echo "sh $0
<ref.fa>
<query.fa>
<outdir>
<prefix>
<cpu for blastn>
"
	exit 1
fi
source /home/myth/miniconda3/install/etc/profile.d/conda.sh
conda activate  offtarget
set -e
ref=$1
q=$2
outdir=$3
prefix=$4
cpu=$5

mkdir -p $outdir
if [ ! -f "$ref.nsq" ];
then
	makeblastdb -in $ref  -dbtype nucl
fi
echo start: blastn -query $q -db $ref -outfmt 6 -num_threads $cpu -task blastn -out $outdir/$prefix.m6
blastn -query $q -db $ref -outfmt 6 -num_threads $cpu -task blastn -out $outdir/$prefix.m6
echo output $outdir/$prefix.m6
