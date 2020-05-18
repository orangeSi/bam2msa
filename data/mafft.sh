cat ref.fa read.se.fa >ref.read.fa
mafft --adjustdirection ref.read.fa >ref.read.fa.mafft.msa
perl -e 'open MSA,"$ARGV[0]";$/=">";<MSA>;while(<MSA>){chomp;my ($id,$seq)=split(/\n/,$_,2);$seq=~ s/\n//g; $seq=uc($seq);print ">$id\n$seq\n"}close MSA'  ref.read.fa.mafft.msa >ref.read.fa.mafft.reformated.msa
