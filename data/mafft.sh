set -vex
cat read.se.fa|head -2 >tmp
rr=tmp
cat ref.fa $rr >tmp1
mafft --adjustdirection tmp1 >tmp1.mafft.msa
perl -e 'open MSA,"$ARGV[0]";$/=">";<MSA>;while(<MSA>){chomp;my ($id,$seq)=split(/\n/,$_,2);$seq=~ s/\n//g; $seq=uc($seq);print ">$id\n$seq\n"}close MSA'  tmp1.mafft.msa >tmp1.reformated.msa


cat read.se.fa|tail -2 >tmp
rr=tmp
cat ref.fa $rr >tmp2
mafft --adjustdirection tmp2 >tmp2.mafft.msa
perl -e 'open MSA,"$ARGV[0]";$/=">";<MSA>;while(<MSA>){chomp;my ($id,$seq)=split(/\n/,$_,2);$seq=~ s/\n//g; $seq=uc($seq);print ">$id\n$seq\n"}close MSA'  tmp2.mafft.msa >tmp2.reformated.msa
cat tmp1.reformated.msa tmp2.reformated.msa >out.mafft.msa
rm tmp1* tmp2*

echo done
