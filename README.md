# doSameSpeciesLiftOver_nextflow

# **!!Deprecated!! see note below**

I’d recommend using either the “Flo” pipeline https://github.com/wurmlab/flo or Liftoff https://github.com/agshumate/Liftoff 

This Nextflow pipeline was never really that robust. 

See here for more details if interested: https://genomic.social/@photocyte/112255774455268103

---

A NextFlow pipeline to lift over GFF files using the UCSC liftover tools. 

Unlike many other liftOver pipelines, which use pre-computed liftover files (e.g. from UCSC), this script generates a custom liftover file by performing blat alignment of the provided "old" and "new" FASTA files.

- Inspired by: [doSameSpeciesLiftOver.pl](https://genome-source.gi.ucsc.edu/gitlist/kent.git/raw/master/src/hg/utils/automation/doSameSpeciesLiftOver.pl)
- And this: [using-liftover-to-convert-genome-assembly-coordinates/](https://iamphioxus.org/2013/06/25/using-liftover-to-convert-genome-assembly-coordinates/)
- Also this: [flo](https://github.com/wurmlab/flo)
- This is a nice general overview: [Griffith Lab - Introduction to liftover tools](http://www.genviz.org//module-01-intro/0001/06/02/liftoverTools/])


## Installation

1. Install [Miniconda3](https://conda.io/en/latest/miniconda.html)
2. Setup conda environment 

```
conda create --name doSameSpeciesLiftOver
conda activate doSameSpeciesLiftOver
```
3. Install conda dependencies:  

```
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority false
conda install nextflow graphviz
```
4. The `doSameSpeciesLiftOver.nf` script will dynamically install the rest of the conda dependencies as needed, but the dependicies can be preinstalled if you'd like. Install using the below line and simply delete the conda directives from the `doSameSpeciesLiftOver.nf` script, or set the `totalCondaEnvPath` parameter to an environment with the dependencies.

```
conda install ucsc-fatotwobit blat ucsc-fasplit ucsc-liftup \
ucsc-axtchain ucsc-chainmergesort ucsc-chainsplit ucsc-chainsort \
seqkit ucsc-chainnet ucsc-netchainsubset ucsc-liftover \
genometools-genometools gffutils
```

## Example

```
nextflow run doSameSpeciesLiftOver.nf \
-resume \
--gff examples/Homo_sapiens.GRCh38_chr6-subset.84.gff3 \
--oldFasta examples/hg38-chr6-subseq1.fa \
--newFasta examples/hg38-chr6-subseq2.fa \
-with-trace examples/trace.txt \
-with-report examples/report.html \
-with-dag examples/flowchart.svg
```
(Currently the example takes about 2 mins to run, with 50% of the computational time being conda installing things and the other 50% the blat alignment)

## Known issues
1. Splitting of FASTA files to a sub-record level (e.g. splitting a record every 4999 bp using 1500 bp overlaps by using non-default values for `params.splitDepth` `params.splitSize` `params.recordSplit` and `params.extraBases`), leads to an incorrect liftover file. I believe the NextFlow logic is correct, so the problem is something I don't quite understand about blat, chainfiles, and/or liftover files.
2. blat isn't multithreaded. Could use [pblat](https://github.com/icebert/pblat) instead.
3. Without sub-record splitting, and without multithreading, the script can't take advantage of parallel computing resources very well (only parallelizes at 1 blat alignment process per FASTA record). So, it takes a looong time for a whole genome liftover. The script is probably best used in a genome assembly tweaking context, where two versions of a single scaffold can be compared and lifted over. The script runs pretty fast in that context.
 - It's possible that the script could only compare/liftover the scaffolds which have a matching record ID (e.g., in different genome assembly versions where the record IDs are the same between the two versions), and also use a hash up front to confirm that the two scaffolds differ before using blat to align them, but not implemented.
4. Features which transverse "N" gaps, are oftentimes not lifted over properly. This is even with no sub-record splitting, and the `-extendThroughN` parameter for blat. The `rescue_unlifted_features` node attempts to fix some of the more trivial edge cases, but doesn't work super well.
5. Could use [CrossMap](http://crossmap.sourceforge.net), to liftover things like .bam, .vcf files, but was quite buggy for me & couldn't get it to work. There is still a vestigal node `crossmap_liftover`, partially implemeting this.
6. blat supports repeat-aware alignment (using the `-mask=` parameter). Including such repeat information would presumably make the alignment faster, but not implemented. Noted that the `constructOocFile` node and the `-ooc` parameter in blat do a sort of simple repeat annotation.


## Workflow flowchart

![Directed acyclic graph (DAG) for doSameSpeciesLiftOver_nextflow program execution](./examples/flowchart.svg)
