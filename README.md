# doSameSpeciesLiftOver_nextflow
A NextFlow pipeline to lift over GFF files using the USCS liftover tools

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