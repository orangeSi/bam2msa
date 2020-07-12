# bam2msa
convert alignment bam/sam file to multiple sequence alignment(msa) file

```
$ cd test/
$ ./bam2msa ref.fa out.bwa.bam --regions NC_045512.2_1bp_to_1680bp:1-10 --display-left-softclip 0 |cut -f 1-6|c1
1                          2               3           4          5           6
#refid                     ref_cut_region  ref_msa     query_id   query_msa   consensus_msa
NC_045512.2_1bp_to_1680bp  1-10            ATTAAAGGTT  clone1     ATTAAAGGTT  ==========
NC_045512.2_1bp_to_1680bp  1-10            ATTAAAGGTT  _R_clone2  ----------  DDDDDDDDDD

$ ./bam2msa ref.fa out.bwa.bam --regions NC_045512.2_1bp_to_1680bp:3-10 --display-left-softclip 0 |cut -f 1-6|c1
1                          2               3         4          5          6
#refid                     ref_cut_region  ref_msa   query_id   query_msa  consensus_msa
NC_045512.2_1bp_to_1680bp  3-10            TAAAGGTT  clone1     TAAAGGTT   ========
NC_045512.2_1bp_to_1680bp  3-10            TAAAGGTT  _R_clone2  --------   DDDDDDDD
```


```
$ ./src/bam2msa
Contact: https://github.com/orangeSi/bam2msa/issues
Usage:
  ./src/bam2msa [flags...] <ref> <bam> [arg...]

convert bam to msa format for alignment file

Flags:
  --display-left-softclip (default: 1)   # display softclip in the start of ref. 0 mean not display, 1 mean display
  --display-right-softclip (default: 1)  # display softclip in the end of ref. 0 mean not display, 1 mean display
  --help                                 # Displays help for the current command.
  --primary-only (default: 1)            # only for primary alignment. 0 mean all alingment, 1 is only primary alignment
  --regions (default: "")                # default display the msa of whole ref for every read. If not, set the specific regsion, ex: chr1:1000-3000
  --version                              # Displays the version of the current application.

Arguments:
  ref (required)                         # ref fasta file
  bam (required)                         # bam alignemnt file

```

```
$ cd test && cat demo.sh 
set -vex
./bam2msa ref.fa out.bwa.bam >out.bam2msa
cat out.bam2msa |grep -v '^#'|awk '{print ">"$1"\n"$2"\n>"$3"\n"$4}' >out.bam2msa.msa

ln -sf  out.bam2msa.msa case
ln -sf  ../data/out.mafft.msa control
diff control case
echo done

$ sh demo.sh
ok, test passed!


```
