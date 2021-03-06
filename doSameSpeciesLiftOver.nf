//Inspired by: https://genome-source.gi.ucsc.edu/gitlist/kent.git/raw/master/src/hg/utils/automation/doSameSpeciesLiftOver.pl
//And this: https://iamphioxus.org/2013/06/25/using-liftover-to-convert-genome-assembly-coordinates/
//Also this: https://github.com/wurmlab/flo
//This is related: https://bioconductor.org/packages/release/bioc/vignettes/CNEr/inst/doc/PairwiseWholeGenomeAlignment.html

//Should preinstall the conda environment, as otherwise the script just spends all its time
//installing dependencies.
//conda activate doSameSpeciesLiftOver
//conda install ucsc-fatotwobit blat ucsc-fasplit ucsc-liftup ucsc-axtchain ucsc-chainmergesort ucsc-chainsplit ucsc-chainsort seqkit ucsc-chainnet ucsc-netchainsubset ucsc-liftover genometools gffutils

//-fastMap in blat: skip the time-consuming stage of merging alignment blocks that have gaps between them

params.oldFasta = "old.fasta"
params.newFasta = "new.fasta"
params.gff = "example.gff3"
params.recordSplit = 1 //Split MultiFasta into files with this many fasta records. Default = 1 .  Higher for less parallelization. 
params.splitSize = 1000000000 //length of a section, in base pairs. 5000 bp is the maximum size allowed for blat -fastmap. Interacts with param.extraBases. Larger to disable.
params.extraBases = 0 //Extra bases for splitFA for overlaps.  E.g. if splitSize is 2500, and extraBases is 2500, splitFA will make 5000 sections that overlap 2500 bp. 0 to disable.
params.splitDepth = 1000000000 //number of subrecord chunks (of length splitSize) per blat invocation. 100 for typical invocation (if chainMerge worked properly). Smaller for more parallelization. Larger to disable.

oldGenome = Channel.fromPath(params.oldFasta)
newGenome = Channel.fromPath(params.newFasta)
gffFiles = params.gff.tokenize(",")
gffFiles_ch = Channel.fromPath(gffFiles).flatten()
gffFiles_ch.into{ gffFile_1 ; gffFile_2 ; gffFile_3 }

oldGenome.into { oldGenome_1 ; oldGenome_2 ; oldGenome_3 ; oldGenome_4 }
newGenome.into { newGenome_1 ; newGenome_2 ; newGenome_3 }

//Split multi-FASTA file into muliple files with typically one FASTA record per file
newGenome_2.splitFasta(by:params.recordSplit,file:true).set{fastaChunks}

process convertFAto2bit_old {
    conda "ucsc-fatotwobit"
//Alternatively, can preinstall all the requirements to a particular location:    
//    conda params.totalCondaEnvPath
 
    tag "$fasta"
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
//    conda params.totalCondaEnvPath
    tag "$old_2bit"
    input:
      file old_2bit from old_2bit_2
    output:
      file "${old_2bit}.ooc" into ooc
    script:
    """
    ##TODO Should follow protocol for repMatch using the "Construct ooc file" instructions from http://genomewiki.ucsc.edu/index.php/DoSameSpeciesLiftOver.pl
    blat ${old_2bit} /dev/null /dev/null -stepSize=1 -tileSize=11 -makeOoc=${old_2bit}.ooc -repMatch=4096
    """
}

process subRecordChunks {
conda "ucsc-fasplit"
//conda params.totalCondaEnvPath
input:
 file fastaChunk from fastaChunks
output:
 set file("${fastaChunk}"), file("${fastaChunk}.lft"),file("${fastaChunk}.subsplit.fa") into subsplitFasta_liftUp
tag "${fastaChunk}"
script:
"""
faSplit size ${fastaChunk} ${params.splitSize} ${fastaChunk}.subsplit -lift=${fastaChunk}.lft -oneFile -extra=${params.extraBases}
"""
}

//Take the 1-per-scaffold FASTA files, which were subsplit internally into ~3000-4000 bp chunks, and separate them into 1-per-section chunks
//In memory splitting.

//Transpose is really the secret sauce below for getting the files setup properly.
subsplitFasta_liftUp.map{ values ->
 subChunks = values[2].splitFasta(by:params.splitDepth,file:true) //split a subrecord multiFASTA into individual FASTAs
 return tuple(values[0],values[1],subChunks)}.transpose().set{ subFastaChunks }

subFastaChunks.combine(old_2bit_1).combine(ooc).set{blatCmds}

process blat_align {
conda "blat ucsc-fasplit ucsc-liftup"
//conda params.totalCondaEnvPath
memory '4 GB'
input:
 set file(originalFasta),file(liftupFile),file(fastaSubChunk),file(old_2bit),file(ooc) from blatCmds
output:
 file "${fastaSubChunk}.lifted.psl" into axtChainCmds
tag "${fastaSubChunk}"
script:
"""
if [ "${params.splitSize}" -lt "4000" ]; then
  blat ${old_2bit} ${fastaSubChunk} -ooc=${ooc} -maxIntron=0 -stepSize=1 -tileSize=11 -minIdentity=98 -noHead -minScore=100 -fastMap -extendThroughN ${fastaSubChunk}.subsplit.psl
else
  blat ${old_2bit} ${fastaSubChunk} -ooc=${ooc} -maxIntron=0 -stepSize=1 -tileSize=11 -minIdentity=98 -noHead -minScore=100 -extendThroughN ${fastaSubChunk}.subsplit.psl
fi

liftUp -pslQ ${fastaSubChunk}.lifted.psl ${liftupFile} warn ${fastaSubChunk}.subsplit.psl
"""
}

process axtChain {
conda "ucsc-axtchain ucsc-fatotwobit"
//conda params.totalCondaEnvPath
input:
 file pslFile from axtChainCmds.collectFile(name:"merged.psl",keepHeader:true,skip:5)
 file oldFasta from oldGenome_3
 file newFasta from newGenome_3
output:
 file "${pslFile}.chain" into chains
tag "${pslFile}"
script:
"""
axtChain -linearGap=loose -faQ -faT -psl ${pslFile} ${oldFasta} ${newFasta} ${pslFile}.chain
"""

}

process chainSortFirst {
conda "ucsc-chainsort"
input:
 file chainFile from chains
output:
 file "sorted.${chainFile}" into sorted_chains
script:
"""
chainSort ${chainFile} sorted.${chainFile}
"""
}

process chainMergeSort_chainSplit {
conda "ucsc-chainmergesort ucsc-chainsplit"
//conda params.totalCondaEnvPath
tag "$chainFile"
input:
 file chainFile from sorted_chains.collect()

output:
 file "chainMerge/*.chain" into sortMergedChains

script:
"""
##chainMerge is the output directory
chainMergeSort ${chainFile} | chainSplit chainMerge stdin -lump=50
"""

}

process chainSortSecond {
conda "ucsc-chainsort"
//conda params.totalCondaEnvPath
tag "$chainFile"
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
//conda params.totalCondaEnvPath
tag "$oldGenome and $newGenome"
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
//conda params.totalCondaEnvPath
tag "$allSortedChain"
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
//conda params.totalCondaEnvPath
publishDir './liftover_output/',mode:'copy',overwrite:true
tag "$netFile & $allSortedChain_2"
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
conda "crossmap=0.3.7"
input:
 file gffFile from gffFile_3
 file liftOverFile from liftOverFile_1
output:
 file "crossmap-lifted_${gffFile}" into crossmap_lifted_gff
script:
"""
crossmap.py -v
crossmap.py gff ${liftOverFile} ${gffFile} crossmap-lifted_${gffFile}
###Below line is to fix a bug in crossmap where it outputs coordinates as floats rather than integers
##Also bug where certain scores for features were set to null?
##cat lifted_unsorted_unfixed_${gffFile} | sed \$'s/.0\t/\t/g' | sed \$'s/\t\t/\t0\t/g'> crossmap-lifted_${gffFile}
###
"""
}

gffFile_1.combine(oldGenome_4).set{normalizeCmds}

process normalizeGff {
publishDir './liftover_output/',mode:'copy',overwrite:true
tag "$gff by $fasta"
input:
 set file(gff),file(fasta) from normalizeCmds
output:
 file "target.${gff}.gff3" optional true into normalizedGff
 file "ignored.${gff}.gff3"
script:
"""
seqkit fx2tab --only-id -n ${fasta} | tr -s "\t" > target_scaffolds.txt
echo "##gff-version 3" >> target_scaffolds.txt
gt gff3 -tidy -sort -retainids -fixregionboundaries ${gff} > normalized.${gff}.gff3 
grep -f target_scaffolds.txt normalized.${gff}.gff3 > target.${gff}.gff3
grep -v -f target_scaffolds.txt normalized.${gff}.gff3 > ignored.${gff}.gff3

if [[ \$(wc -l <target.${gff}.gff3) -le 1 ]]
then
    echo "No targets were found. Deleting target gff3 file."
    rm -f target.${gff}.gff3
fi
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

normalizedGff.combine(liftOverFile_2).set{ liftoverCmds }

process ucsc_liftover {
conda "ucsc-liftover"
//conda params.totalCondaEnvPath
tag "$gffFile & liftOverFile"
input:
 set file(gffFile),file(liftOverFile) from liftoverCmds
output:
 file "original_${gffFile}" into gffOriginal
 file "ucsc-lifted_${gffFile}" into ucsc_lifted_gff, ucsc_lifted_gff_ch2
 file "unmapped_${gffFile}" into unmapped_gff, unmapped_gff_ch2
script:
"""
liftOver -gff ${gffFile} ${liftOverFile} ucsc-lifted_${gffFile} unmapped_${gffFile}
ln -s ${gffFile} original_${gffFile}
"""
}

process rescue_unlifted_features {
tag "$unmapped_gff_ch2"
input:
 file gffOriginal
 file ucsc_lifted_gff_ch2
 file unmapped_gff_ch2
output:
 file "rescued.gff" into rescuedGff
script:
"""
#!/usr/bin/env python
##Oftentimes gene and mRNA features are not lifted over, as they presumably span regions that changed
##Or also possible, they span "NNNN" regions that blat just doesn't want to deal with
##This node takes those unmapped parent features, checks if they originally simply spanned the extent
##of their child features, and if so simply updates the unmapped extent to the new extent of the child features
##1) Load unmapped features. Get the ID(s)
##2) Load the original file, check if the extent of the feature is == to extent of the children
##3) Load the lifted over file, find the aforementioned child ID.  Use the start/end of the left/right child feature respectively. 
import glob
import re
lifted_path = glob.glob("./ucsc-lifted_*")[0]
unmapped_path = glob.glob("./unmapped_*")[0]
original_path = glob.glob("./original_*")[0]
unmapped_ids = dict()

##Load unmapped features
handle = open(unmapped_path,"r")
for l in handle.readlines():
    if l[0] == "#":
        ##Skip comment lines
        continue
    splitline = l.split("\t")
    start = int(splitline[3]) ##Start of the feature
    end = int(splitline[4]) ##End of the feature
    re_result = re.search("ID=(.+)[;\$]",splitline[8])
    if re_result == None:
        continue

    feature_ID = re_result.group(1)
    unmapped_ids[feature_ID] = dict()
    unmapped_ids[feature_ID]["line"] = l 
    unmapped_ids[feature_ID]["type"] = splitline[2]
    unmapped_ids[feature_ID]["start"] = start
    unmapped_ids[feature_ID]["end"] = end
    unmapped_ids[feature_ID]["left_extent_defined_by_child"] = False
    unmapped_ids[feature_ID]["right_extent_defined_by_child"] = False
    unmapped_ids[feature_ID]["left_child"] = None
    unmapped_ids[feature_ID]["right_child"] = None
    unmapped_ids[feature_ID]["childs_scaffold"] = None
handle.close()

##Find the child features.
##Make a mapping of child IDs to parent IDs
child_ids = []
child_to_parent = dict()
handle = open(original_path,"r")
for l in handle.readlines():
    if l[0] == "#":
        continue
    splitline = l.split("\t")
    re_result = re.search("Parent=(.+)[;\$]",splitline[8])
    if re_result == None:
        continue
    parent_ID = re_result.group(1)
    re_result = re.search("ID=(.+)[;\$]",splitline[8])
    if re_result == None:
        continue
    child_ID = re_result.group(1)
    if parent_ID not in unmapped_ids.keys():
        continue
    start = int(splitline[3])
    end = int(splitline[4])
    scaffold = splitline[0]
    if start == unmapped_ids[parent_ID]["start"]:
        unmapped_ids[parent_ID]["left_extent_defined_by_child"] = True
        unmapped_ids[parent_ID]["left_child"] = child_ID 
        child_ids.append(child_ID)
        child_to_parent[child_ID] = parent_ID
    if end == unmapped_ids[parent_ID]["end"]:
        unmapped_ids[parent_ID]["right_extent_defined_by_child"] = True
        unmapped_ids[parent_ID]["right_child"] = child_ID 
        child_ids.append(child_ID)
        child_to_parent[child_ID] = parent_ID
handle.close()


handle = open(lifted_path,"r")
for l in handle.readlines():
    if l[0] == "#":
        continue
        print("skipping due to comment line...")
    splitline = l.split("\t")
    re_result = re.search("ID=(.+)[;\$]",splitline[8])
    if re_result == None:
        continue
        print("skipping due to no ID found...")
    feature_ID = re_result.group(1)
    if feature_ID not in child_ids:
        print("skipping due to feature_ID not in child_ids...")
        continue
    start = int(splitline[3])
    end = int(splitline[4])
    scaffold = splitline[0]
    print(scaffold)
    unmapped_ids[child_to_parent[feature_ID]]["childs_scaffold"] = scaffold
    if unmapped_ids[child_to_parent[feature_ID]]["left_child"] == feature_ID:
        unmapped_ids[child_to_parent[feature_ID]]["start"] = start
    if unmapped_ids[child_to_parent[feature_ID]]["right_child"] == feature_ID:
        unmapped_ids[child_to_parent[feature_ID]]["end"] = end
handle.close()
write_handle = open("rescued.gff","w")
for k in unmapped_ids.keys():
     unmapped_line = unmapped_ids[k]["line"]
     unmapped_splitline = unmapped_line.split("\t")
     if unmapped_ids[k]["childs_scaffold"] == None:
         theScaffold = unmapped_splitline[0]
     else:
         theScaffold = unmapped_ids[k]["childs_scaffold"]
     newline = "\t".join([theScaffold,unmapped_splitline[1],unmapped_splitline[2],str(unmapped_ids[k]["start"]),str(unmapped_ids[k]["end"]),unmapped_splitline[5],unmapped_splitline[6],unmapped_splitline[7],unmapped_splitline[8]])
     write_handle.write(newline)
"""
}

//ucsc_lifted_gff.mix(crossmap_lifted_gff).set{liftedGffUnsorted}

process sort_gff {
publishDir './liftover_output/',mode:'copy',overwrite:true
conda "genometools-genometools"
//conda params.totalCondaEnvPath
input:
 file gff from ucsc_lifted_gff
 file unmapped from unmapped_gff
 file rescued from rescuedGff
output:
 file "srt_${gff}" into final_gff
 file "${unmapped}"
tag "${gff}" 
script:
"""
 THENAME=${gff}
 NEWNAME=lifted_\${THENAME#lifted_unsorted_}
 cat ${gff} ${rescued} | grep -v "#" | gt gff3 -tidy -sort -retainids > srt_${gff}
"""
}

process compare_gffs {
echo true
tag "$ogff vs. $fgff"
input:
 file ogff from gffFile_2
 file fgff from final_gff
script:
"""
cat ${ogff} | grep -v "#" | cut -f 3 | sort | uniq > feature_types.txt
echo "Done with pipeline. Check the liftover_output folder"
"""
}
