# bam2msa
convert alignment file bam to multiple sequence alignment(msa)

```
$./src/bam2msa
Contact: ilikeorangeapple@gmail.com or go to https://github.com/orangeSi/grepfile/issues
Usage:
  ./src/bam2msa [flags...] <ref> <bam> [arg...]

convert bam to msa format

Flags:
  --help                       # Displays help for the current command.
  --primary-only (default: 1)  # only for primary alignment. 0 mean all alingment, 1 is only primary alignment
  --regions (default: "")      # default display the msa of whole ref for every read. If not, set the specific regsion, ex: chr1:100-300,chr3:500-800
  --version                    # Displays the version of the current application.

Arguments:
  ref (required)               # ref fasta file
  bam (required)               # bam alignemnt file


$cd test && cat demo.sh 
set -vex
./bam2msa ref.fa out.bwa.bam >out.bam2msa
cat out.bam2msa |grep -v '^#'|awk '{print ">"$1"\n"$2"\n>"$3"\n"$4}' >out.bam2msa.msa

ln -sf  out.bam2msa.msa case
ln -sf  ../data/out.mafft.msa control
diff control case
echo done

$ sh demo.sh

```
