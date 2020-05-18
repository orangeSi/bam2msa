set -vex
./bam2msa ref.fa out.bam >out.bam2msa
cat out.bam2msa |awk '{if(NR==2){print ">"$1"\n"$2"\n>"$3"\n"$4};if(NR==3){print ">"$3"\n"$4}}' >out.bam2msa.msa

ln -sf  out.bam2msa.msa case
ln -sf  ../data/ref.read.fa.mafft.reformated.msa control
diff control case
echo done
