//Inspired by: https://genome-source.gi.ucsc.edu/gitlist/kent.git/raw/master/src/hg/utils/automation/doSameSpeciesLiftOver.pl
//And this: 
//Also this: https://iamphioxus.org/2013/06/25/using-liftover-to-convert-genome-assembly-coordinates/
params.old = "old.fasta"
params.new = "new2.fasta"
params.gff = "Phengodes_nigromaculata.gff3"
params.splitDepth = 100000000 //number of sections per blat invocation. 100 for typical invocation (if chainMerge worked properly). Smaller for more parallelization. Larger to disable.
params.splitSize = 100000000 //length of a section, in base pairs. 5000 bp is the maximum size allowed for blat -fastmap. Interacts with param.extra. Larger to disable.
params.recordSplit = 1 //Split MultiFasta into files with this many fasta records. Default = 1 .  Higher for less parallelization. 
params.extra = 0 //Extra bases for splitFA for overlaps.  E.g. if splitsize is 2500, and extra is 2500, splitFA will make 5000 sections that overlap 2500 bp. 0 to disable.

oldGenome = Channel.fromPath(params.old)
newGenome = Channel.fromPath(params.new)
gffFile = Channel.fromPath(params.gff)

gffFile.into{ gffFile_1 ; gffFile_2 }

oldGenome.into {oldGenome_1 ; oldGenome_2 ; oldGenome_3 ; oldGenome_4 }
newGenome.into { newGenome_1 ; newGenome_2 ; newGenome_3 }

//Split multi-FASTA file into muliple files with typically one FASTA record per file
newGenome_2.splitFasta(by:params.recordSplit,file:true).set{fastaChunks}

process convertFAto2bit_old {
    conda "ucsc-fatotwobit"
    input:
    file fasta from oldGenome_1

    output:
    file "${fasta}.2bit" into old_2bit 
    script:
    """
    faToTwoBit ${fasta} ${fasta}.2bit    
    """
}

old_2bit.into{ old_2bit_1 ; old_2bit_2 }

process constructOocFile {
    conda "blat"
    input:
      file old_2bit from old_2bit_2
    output:
      file "${old_2bit}.ooc" into ooc
    script:
    """
    ##TODO Should follow protocol for repMatch using the "Construct ooc file" instructions from http://genomewiki.ucsc.edu/index.php/DoSameSpeciesLiftOver.pl
    blat ${old_2bit} /dev/null /dev/null -tileSize=11 -makeOoc=${old_2bit}.ooc -repMatch=300
    """
}

process evenSmallerChunks {
conda "ucsc-fasplit"
input:
 file fastaChunk from fastaChunks
output:
 set file("${fastaChunk}"), file("${fastaChunk}.lft"),file("${fastaChunk}.subsplit.fa") into subsplitFasta_liftUp
tag "${fastaChunk}"
script:
"""
faSplit size ${fastaChunk} ${params.splitSize} ${fastaChunk}.subsplit -lift=${fastaChunk}.lft -oneFile -extra=${params.extra}
"""
}

//Take the 1-per-scaffold FASTA files, which were subsplit internally into ~3000-4000 bp chunks, and separate them into 1-per-section chunks
//In memory splitting.

//Transpose is really the secret sauce below for getting the files setup properly.
subsplitFasta_liftUp.map{ values ->
 subChunks = values[2].splitFasta(by:params.splitDepth,file:true)
 return tuple(values[0],values[1],subChunks)}.transpose().set{ subFastaChunks }

subFastaChunks.combine(old_2bit_1).combine(ooc).set{blatCmds}

process blat_align {
conda "blat ucsc-fasplit ucsc-liftup"
memory '4 GB'
input:
 set file(originalFasta),file(liftupFile),file(fastaSubChunk),file(old_2bit),file(ooc) from blatCmds
output:
 file "${fastaSubChunk}.psl" into axtChainCmds
tag "${fastaSubChunk}"
script:
"""
if [ "${params.splitSize}" -lt "4000" ]; then
  blat ${old_2bit} ${fastaSubChunk} -ooc=${ooc} -tileSize=11 -minIdentity=98 -noHead -minScore=100 -fastMap -extendThroughN ${fastaSubChunk}.subsplit.psl
else
  blat ${old_2bit} ${fastaSubChunk} -ooc=${ooc} -tileSize=11 -minIdentity=98 -noHead -minScore=100 -extendThroughN ${fastaSubChunk}.subsplit.psl
fi

liftUp -pslQ ${fastaSubChunk}.psl ${liftupFile} warn ${fastaSubChunk}.subsplit.psl

##cleanup some temporary files
##rm -f ${fastaSubChunk}.subsplit.psl ${fastaSubChunk}.subsplit.fa
"""
}

process axtChain {
conda "ucsc-axtchain ucsc-fatotwobit"
input:
 file pslFile from axtChainCmds.collectFile(name:"merged.psl",keepHeader:true,skip:5)
 file oldFasta from oldGenome_3
 file newFasta from newGenome_3
output:
 file "${pslFile}.chain" into chains
tag "${pslFile}"
script:
"""
axtChain -linearGap=medium -faQ -faT -psl ${pslFile} ${oldFasta} ${newFasta} ${pslFile}.chain
"""

}

process chainMerge {
conda "ucsc-chainmergesort ucsc-chainsplit"
input:
 file chainFile from chains.collect()

output:
 file "chainMerge/*.chain" into sortMergedChains

script:
"""
##chainMerge is the output directory
chainMergeSort ${chainFile} | chainSplit chainMerge stdin -lump=50
"""

}

process chainSort {
conda "ucsc-chainsort"
input:
 file chainFile from sortMergedChains.collectFile(name: 'all.chain')

output:
 file "all.sorted.chain" into allSortedChain_1, allSortedChain_2
script:
"""
chainSort all.chain all.sorted.chain
"""

}

process calculateChromInfo {
conda "seqkit"
input:
 file oldGenome from oldGenome_2
 file newGenome from newGenome_1

output:
 set file("${oldGenome}.chromInfo"),file("${newGenome}.chromInfo") into chromInfos
script:
"""
 ##Equivalent command that can be run on FASTA files:
 seqkit fx2tab --only-id -nl ${oldGenome} | tr -s "\t" | sort -k2,2nr > ${oldGenome}.chromInfo
 seqkit fx2tab --only-id -nl ${newGenome} | tr -s "\t" | sort -k2,2nr > ${newGenome}.chromInfo

 ##Old way that used ucsc-twobitinfo from a 2bit file.
 ##twoBitInfo new.2bit new.2bit.chromInfo
 ##twoBitInfo old.2bit old.2bit.chromInfo
"""

}

process chainNet {
conda "ucsc-chainnet"
input:
 file allSortedChain from allSortedChain_1
 set file(oldInfo),file(newInfo) from chromInfos
output:
 file "all.net" into netFile
script:
"""
chainNet ${allSortedChain} ${oldInfo} ${newInfo} all.net /dev/null
"""
}

process produceLiftOverFile {
conda "ucsc-netchainsubset"
publishDir './liftover_output/',mode:'copy',overwrite:true
input:
 file netFile
 file allSortedChain_2
output:
 file "final.liftOver" into liftOverFile_2
script:
"""
netChainSubset ${netFile} ${allSortedChain_2} final.liftOver
"""
}

//process crossmap_liftover {
//conda "crossmap=0.3.6"
//input:
// file gffFile from gffFile_1
// file liftOverFile from liftOverFile_1
//output:
// file "crossmap-lifted_${gffFile}" into crossmap_lifted_gff
//script:
//"""
//Crossmap.py -v
//CrossMap.py gff ${liftOverFile} ${gffFile} lifted_unsorted_unfixed_${gffFile}
//###Below line is to fix a bug in crossmap where it outputs coordinates as floats rather than integers
//##Also bug where certain scores for features were set to null?
//cat lifted_unsorted_unfixed_${gffFile} | sed \$'s/.0\t/\t/g' | sed \$'s/\t\t/\t0\t/g'> crossmap-lifted_${gffFile}
//###
//"""
//}

gffFile_1.combine(oldGenome_4).set{normalizeCmds}

process normalizeGff {
publishDir './liftover_output/',mode:'copy',overwrite:true
input:
 set file(gff),file(fasta) from normalizeCmds
output:
 file "target.${gff}.gff3" into normalizedGff
 file "ignored.${gff}.gff3"
script:
"""
seqkit fx2tab --only-id -n ${fasta} | tr -s "\t" > target_scaffolds.txt
echo "##gff-version 3" >> target_scaffolds.txt
gt gff3 -tidy -sort -retainids -fixregionboundaries ${gff} > normalized.${gff}.gff3 
grep -f target_scaffolds.txt normalized.${gff}.gff3 > target.${gff}.gff3
grep -v -f target_scaffolds.txt normalized.${gff}.gff3 > ignored.${gff}.gff3
"""
}

//process ucsc_gff3togenepred {
//conda "ucsc-gff3togenepred genometools"
//input:
// file gffFile from gffFile_1
//output:
// file "${gffFile}.gp" into gpFile
//script:
//"""
//
//gt gff3 -tidy -sort -retainids -fixregionboundaries ${gffFile} > tmp.gff3
//gff3ToGenePred tmp.gff3 ${gffFile}.gp
//"""
//}

process ucsc_liftover {
conda "ucsc-liftover"
input:
 file gffFile from normalizedGff
 file liftOverFile from liftOverFile_2
output:
 file "ucsc-lifted_${gffFile}" into ucsc_lifted_gff
 file "unmapped_${gffFile}" into unmapped_gff
script:
"""
liftOver -gff ${gffFile} ${liftOverFile} ucsc-lifted_${gffFile} unmapped_${gffFile}
"""
}

//ucsc_lifted_gff.mix(crossmap_lifted_gff).set{liftedGffUnsorted}

process sort_gff {
publishDir './liftover_output/',mode:'copy',overwrite:true
conda "genometools"
input:
 file gff from ucsc_lifted_gff
 file unmapped from unmapped_gff
output:
 file "srt_${gff}" into final_gff
 file "${unmapped}"
tag "${gff}" 
script:
"""
 THENAME=${gff}
 NEWNAME=lifted_\${THENAME#lifted_unsorted_}
 cat ${gff} | grep -v "#" | gt gff3 -tidy -sort -retainids | grep -v "###" > srt_${gff}
"""
}

process compare_gffs {
echo true
input:
 file ogff from gffFile_2
 file fgff from final_gff
script:
"""
cat ${ogff} | grep -v "#" | cut -f 3 | sort | uniq > feature_types.txt
"""
}
