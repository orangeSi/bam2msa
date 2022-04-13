# bam2msa
convert alignment bam file to multiple sequence alignment(msa) file

```
$ cd test/

$ ./bam2msa ref.fa out.bwa.bam NC_045512.2_1bp_to_1680bp:1-20|cut -f 1-7|c1
1                     2                     3                     4                               5         6     7
#query_msa            ref_msa               consensus_msa         ref_cut_region                  query_id  FLAG  POS_in_Bam
ATTAAAGGTTTATACC--CC  ATTAAAGGTTTATACCTTCC  ================DD==  NC_045512.2_1bp_to_1680bp:1-20  clone1    0     1

```


```
$ ../src/bam2msa
Usage:
  ./bam2msa [flags...] <ref> <bam> <regions> [arg...]

convert bam to msa format for alignment file

Flags:
  --help                                      # Displays help for the current command.
  --primary-only (default: 1)                 # only for primary alignment. 0 mean all alingment, 1 is only primary alignment
  --span-whole-region-read-only (default: 1)  # only for read which span the whole region. 0 mean all read which overlap with the region, 1 mean is read which span the whole region
  --version                                   # Displays the version of the current application.

Arguments:
  ref (required)                              # ref fasta file
  bam (required)                              # bam alignemnt file
  regions (required)                          # display read and ref msa alignment in these regions, example: chr1:1000-1200,chr2:2000-2300

```

```
$ cd test && cat demo.sh 
set -e
./bam2msa ref.fa out.bwa.bam NC_045512.2_1bp_to_1680bp:1-1680 --span-whole-region-read-only 0 >out.bam2msa
cat out.bam2msa |grep -v '^#'|awk -F '\t' '{print ">"$4"\n"$2"\n>"$5"\n"$1}'|sed 's/:.*//' >out.bam2msa.msa

ln -sf  out.bam2msa.msa case
ln -sf  ../data/out.mafft.msa control
diff control case
if [ $? -eq 0 ];
then
	echo "ok, test passed!"
else
	echo "sorry, test failed! give a issue to me please~~~"
fi

$ sh demo.sh
ok, test passed!

```
## for wasm support
```
samtools view data/out.bwa.bam|wasmer  --dir=data/ src/bam2msa_final.wasm -- data/ref.fa STDIN NC_045512.2_1bp_to_1680bp:1-1680 
```
## colorsize snp/indel of output with --colorize-snp-indel 1
```
bam2msa test/ref.fa test/out.bwa.bam NC_045512.2_1bp_to_1680bp:1-86 --span-whole-region-read-only 0 --colorize-snp-indel 1 |column -ts $'\t'|less -RS
```
