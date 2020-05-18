set -vex
ref=ref.fa
q=read.se.fa
blastn.sh $q $q . out.blastn.toself 2
blastn.sh $ref $q . out.blastn 2
bwa.sh $ref $q . out.bwa 2 bam
