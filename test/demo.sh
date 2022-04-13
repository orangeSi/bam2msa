set -e
./bam2msa ref.fa out.bwa.bam NC_045512.2_1bp_to_1680bp:1-1680 --span-whole-region-read-only 0 --display-read-boundary 0 >out.bam2msa
cat out.bam2msa |grep -v '^#'|awk -F '\t' '{print ">"$4"\n"$2"\n>"$5"\n"$1}'|sed 's/:.*//' >out.bam2msa.msa

ln -sf  out.bam2msa.msa case
sed 's/^>_R_/>/'  ../data/out.mafft.msa > control
diff control case
if [ $? -eq 0 ];
then
	echo "ok, test passed!"
else
	echo "sorry, test failed! give a issue to me please~~~"
fi
