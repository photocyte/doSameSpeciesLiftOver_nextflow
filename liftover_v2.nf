params.old = "old.Pnm1.1_ph1.fasta"
params.new = "new.Pnm1.1_ph1.fasta"
params.gff = "Pnm1.1s2685p1_gmap_align.out.gff3"
oldGenome = Channel.fromPath(params.old)
newGenome = Channel.fromPath(params.new)
gffFile = Channel.fromPath(params.gff)

oldGenome.into {oldGenome_1 ; oldGenome_2}
newGenome.into {newGenome_1 ; newGenome_2 ; newGenome_3 }
gffFile.into{ gffFile_1 ; gffFile_2 }

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

process convertFAto2bit_new {
    conda "ucsc-fatotwobit"
    input:
    file fasta from newGenome_1

    output:
    file "${fasta}.2bit" into new_2bit

    script:
    """
    faToTwoBit ${fasta} ${fasta}.2bit    
    """
}

old_2bit.into{ old_2bit_1 ; old_2bit_2 ; old_2bit_3}
new_2bit.into{ new_2bit_1 ; new_2bit_2 }

process constructOocFile {
    conda "blat"
    input:
      file old_2bit from old_2bit_3
    output:
      file "${old_2bit}.ooc" into ooc
    script:
    """
    blat ${old_2bit} /dev/null /dev/null -tileSize=11 -makeOoc=${old_2bit}.ooc -repMatch=300
    """
}

process blat_align {
conda "blat"
input:
 file old_2bit from old_2bit_1
 file newFasta from newGenome_3
 file ooc
output:
file "*.psl" into psls
file "${newFasta}" into blatFasta
file "${old_2bit}" into old2bitToaxtChain
script:
"""
blat ${old_2bit} ${newFasta} -ooc=${ooc} -tileSize=11 -minIdentity=98 ${newFasta}.psl -noHead -minScore=100
"""
}

process axtChain {
conda "ucsc-axtchain ucsc-fatotwobit"
input:
 file pslFile from psls
 file old_2bit from old2bitToaxtChain
 file newFasta from blatFasta
output:
 file "*.chain" into chains

script:
"""
faToTwoBit ${newFasta} ${newFasta}.2bit
axtChain -linearGap=medium -psl ${pslFile} ${old_2bit} ${newFasta}.2bit ${pslFile}.chain
"""

}

process chainMerge {
conda "ucsc-chainmergesort ucsc-chainsplit"
input:
 file chainFile from chains

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
 file chainFile from sortMergedChains

output:
 file "all.sorted.chain" into allSortedChain_1, allSortedChain_2
script:
"""
cat *.chain > all.chain
chainSort all.chain all.sorted.chain
"""

}

process calculateChromInfo {
conda "seqkit"
input:
 file oldGenome from oldGenome_2
 file newGenome from newGenome_2

output:
 file "${oldGenome}.chromInfo" into oldInfo
 file "${newGenome}.chromInfo" into newInfo
script:
"""
 ##Equivalent command that can be run on FASTA files:
 seqkit fx2tab -nl ${oldGenome} | tr -s "\t" | sort -k2,2nr > ${oldGenome}.chromInfo
 seqkit fx2tab -nl ${newGenome} | tr -s "\t" | sort -k2,2nr > ${newGenome}.chromInfo

 ##Old way that used ucsc-twobitinfo from a 2bit file.
 ##twoBitInfo ${new_2bit} ${new_2bit}.chromInfo
 ##twoBitInfo ${old_2bit} ${old_2bit}.chromInfo
"""

}

process chainNet {
conda "ucsc-chainnet"
input:
 file allSortedChain from allSortedChain_1
 file oldInfo
 file newInfo
output:
 file "all.net" into netFile
script:
"""
chainNet ${allSortedChain} ${oldInfo} ${newInfo} all.net /dev/null
"""
}

process produceLiftOverFile {
conda "ucsc-netchainsubset"
input:
 file netFile
 file allSortedChain_2
output:
 file "final.liftOver" into liftOverFile_1, liftOverFile_2
script:
"""
netChainSubset ${netFile} ${allSortedChain_2} final.liftOver
"""
}

process crossmap_liftover {
conda "crossmap=0.3.6"
input:
 file gffFile from gffFile_1
 file liftOverFile from liftOverFile_1
output:
 file "crossmap-lifted_${gffFile}" into crossmap_lifted_gff
script:
"""
Crossmap.py -v
CrossMap.py gff ${liftOverFile} ${gffFile} lifted_unsorted_unfixed_${gffFile}
###Below line is to fix a bug in crossmap where it outputs coordinates as floats rather than integers
##Also bug where certain scores for features were set to null?
cat lifted_unsorted_unfixed_${gffFile} | sed \$'s/.0\t/\t/g' | sed \$'s/\t\t/\t0\t/g'> crossmap-lifted_${gffFile}
###
"""
}

process ucsc_liftover {
conda "ucsc-liftover"
input:
 file gffFile from gffFile_2
 file liftOverFile from liftOverFile_2
output:
 file "ucsc-lifted_${gffFile}" into ucsc_lifted_gff
 file "unmapped_${gffFile}"
script:
"""
liftOver -gff ${gffFile} ${liftOverFile} ucsc-lifted_${gffFile} unmapped_${gffFile}
"""
}

ucsc_lifted_gff.mix(crossmap_lifted_gff).set{liftedGffUnsorted}

process sort_gff {
publishDir './output/'
conda "genometools"
input:
 file liftedGffUnsorted
output:
 file "srt_${liftedGffUnsorted}"
tag "${liftedGffUnsorted}" 
script:
"""
 THENAME=${liftedGffUnsorted}
 NEWNAME=lifted_\${THENAME#lifted_unsorted_}
 gt gff3 -tidy -sort -retainids <(cat ${liftedGffUnsorted} | grep -v "#") > srt_${liftedGffUnsorted}
"""
}

