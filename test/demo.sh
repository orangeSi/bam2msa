set -e
./bam2msa ref.fa out.bwa.bam >out.bam2msa
cat out.bam2msa |grep -v '^#'|awk '{print ">"$1"\n"$3"\n>"$4"\n"$5}' >out.bam2msa.msa

ln -sf  out.bam2msa.msa case
ln -sf  ../data/out.mafft.msa control
diff control case
if [ $? -eq 0 ];
then
	echo "ok, test passed!"
else
	echo "sorry, test failed! give a issue to me please~~~"
fi
