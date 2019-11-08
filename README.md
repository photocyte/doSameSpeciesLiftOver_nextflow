# doSameSpeciesLiftOver_nextflow
A NextFlow pipeline to lift over GFF files using the USCS liftover tools

### Installation

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
4. The `doSameSpeciesLiftOver.nf` script will dynamically install the rest of the conda dependencies

### Example

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

### Workflow flowchart

![Directed acyclic graph (DAG) for doSameSpeciesLiftOver_nextflow program execution](./examples/flowchart.svg)