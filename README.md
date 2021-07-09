# This file provides a general workflow for identifying TE junctions for a given Drosophila cell line

## This section describes the process of identifying TE junctions without prior information about the location of the TEs (slower)

Broadly speaking there here are the main steps in the process (each described in more detail below):
-	Mask reference genome for transposons of interest
-	Trim the reads for adapter/low quality bases
-	Map reads to transposon sequences
-	Demultiplex sample data based on transposon mapping
-	Map the demultiplexed reads to the genome and transposons
-	Identify informative reads that span TE junctions 
-	Aggregate informative reads to call TE junctions
-	Identify and describe putative transposable elements
-	Create a presence/absence tab-delimited file
-	Visualize and cluster the binary presence/absence file

### Mask reference genome for known transposons of interest
Masking of reference genome was performed with ncbi blast-2.2.26:<br>
```
blastall -i <transposon_consensus.fasta> -d <reference_genome.fasta> -p blastn -a 10 -e 1e-100 -F "m L" -U T -K 20000 -b 20000 -m 8 | perl maskRef.pl <reference_genome.fasta> - > <masked_reference_genome.fasta>
```

### Trim reads for adapter/low quality bases
There is no specific process for this, use your favorite trimming software. Examples that we use regularly are Trimmomatic or fastp

### Map reads to transposon sequences
Used your favorite mapping software to map R2 reads to the transposon consensus sequences. We used bowtie2 to do the mapping. R2 reads should be internal to the transposon.

### Demultiplex sample data based on transposon mapping
Based on the mapping of the R2 reads to the transposon sequences, we demultiplex the paired-end reads into the individual files tagged with the transposon that the reads are affiliated with.<br>
This is performed using the demultiplexPairedSamples.pl script, here is an example command line:<br>
perl demultiplexPairedSamples.pl <transposon_mapping_result.sam> <R1_fastq_file> <R2_fastq_file> <output_path_and_basename><br>
The final filename(s) will be the <output_path_and_basename>.<id_of_hit_transposon>.fastq.gz<br>
Here is an example commandline:<br>
```
perl demultiplexPairedSamples.pl transposonMapping/sample1.sam.gz trimmed/sample1_R1_001.paired.trimmed.fastq trimmed/ sample1_S1_R2_001.paired.trimmed.fastq demultiplexed/sample1
```
The output files generated would look like:<br>
```
demultiplexed/sample1.1731.R1.fastq.gz<br>
demultiplexed/sample1.1731.R2.fastq.gz<br>
demultiplexed/sample1.297.R1.fastq.gz<br>
demultiplexed/sample1.297.R2.fastq.gz<br>
demultiplexed/sample1.copia.R1.fastq.gz<br>
demultiplexed/sample1.copia.R2.fastq.gz<br>
demultiplexed/sample1.mdg1.R1.fastq.gz<br>
demultiplexed/sample1.mdg1.R2.fastq.gz<br>
demultiplexed/sample1.roo.R1.fastq.gz<br>
demultiplexed/sample1.roo.R2.fastq.gz<br>
```

### Map the demultiplexed reads to the genome and transposons
Every instance of the demultiplexed R1 reads need to be mapped to the transposon masked reference genome, in both the forward and reverse-complement orientation, and mapped to the transposon consensus sequences to identify reads that identify the junction between transposons and the genome. Reads must be mapped with soft clipping turned on to identify these junction reads.<br>
Mapping can be done with any mapping tool that meets these requirements and produces SAM formatted output. Here are the example commands for:<br>
```
bowtie2 -x masked_reference_forward –local –no-head -k 2 -U demultiplexed/sample1.copia.R1.fastq.gz | gzip -c > sam/sample1.copia.sam.dmel.gz<br>
bowtie2 -x masked_reference_reverse –local –no-head -k 2 -U demultiplexed/sample1.copia.R1.fastq.gz | gzip -c > sam/sample1.copia.sam.dmel.rc.gz<br>
bowtie2 -x transposon_consensus –local –no-head -k 2 -U demultiplexed/sample1.copia.R1.fastq.gz | gzip -c > sam/sample1.sam.TEs.gz<br>
```

### Identify informative reads that span TE junctions
This step takes the search results and looks for reads that straddle transposon boundary and outputs these for aggregation and transposon site calling. It works on existing and previously unidentified transposons insertions.
Example commands – they come in pairs (for the forward search and the reverse-complement search):
```
perl identifyGenomicTEJunctions.pl 108 5 5 xxx.log 1731 sam/sample1.copia.sam.dmel.gz sam/sample1.copia.sam.TEs.gz | sort -k 2,2 -k 3,3n > candidate/sample1.copia.for.cand.txt<br>
perl identifyGenomicTEJunctions.pl 108 5 5 xxx.log 1731 sam/sample1.copia.sam.dmel.rc.gz sam/sample1.copia.sam.TEs.gz | sort -k 2,2 -k 3,3n > candidate/sample1.copia.rev.cand.txt<br>
```

### Generate read count summary file
Given all the initially parsed mapping data is in a set of files, this command will take as input that list of files and generate a summary of the read count information.
```
grep TEs candidate/*txt | perl generateReadCountSummaryFile.pl - > readCountSummary.txt<br>
```

### Aggregate informative reads to call TE junctions
Taking the junction associated reads, data is aggregated to identify high confidence junctions. High confidence junctions have a number of independent reads starting at a number of different positions on the genome; this is done to filter out/exclude junctions-like cases resulting from amplification artifacts such as PCR-duplicates.
```
perl findJunctions3.pl readCountSummary.txt <transposon_type> <algorithm> candidate/sample1.copia.for.cand.txt candidate/sample1.copia.for.cand.txt candidate/sample1.copia.rev.cand.txt candidate/sample1.copia.rev.cand.txt masked_reference_forward.fasta | sort -k 1,1 -k 2,2n > sample1.copia.candidate.junctions.txt
```

### Predict transposable elements
Take the predicted junctions and consider whether they can be combined to identify intact transposable elements in the genome.
```
perl predictElement.pl <candidate_junction_file> <masked_reference_forward.fasta> > <results_file>
```

### Create a presence/absence tab-delimited file
Generating a table of the results files containing the full path to a particular result file, the sample id to use for this file, and the transposon name. This should be tab-delimited and is used as input.<br>
Run the following command to convert the results files into a binary table of results for clustering and visualization:
```
perl convert2Presence.pl <results_table> > combined.presence.tsv
```

### Cluster and visualize the results
Clustering and visualization was carried out in R with the following commands:
```
library("gplots")
data<-read.table("combined.presence.tsv",sep="\t",header=TRUE,row.names=1)
my_palette <- colorRampPalette(c("black", "green"))(n = 2)
heatmap.2(x=t(as.matrix(data)),distfun = function(x) dist(x,method='binary'),scale="column",dendrogram = c("row"),labCol=FALSE,trace="none",margins=c(5,20),key=FALSE, col=my_palette)
```
If clustering numerical data the following commands would be used:
```
library("gplots")
data<-read.table("combined.count.tsv",sep="\t",header=TRUE,row.names=1)
my_palette <- colorRampPalette(c("black", "black", "black, "black", "black", "green, "green", "green"))(n = 299)
heatmap.2(x=t(as.matrix(data)),distfun = function(x) dist(x,method='manhattan'),scale="column",dendrogram = c("row"),labCol=FALSE,labRow=colnames(data),trace="none",margins=c(5,20),key=TRUE, col=my_palette)
```

## This section describes the process of quantifying TEs given previously determined TE insertions (faster)

Broadly speaking there here are the main steps in the process (each described in more detail below):
-	Isolate sequences flanking established transposons
-	Trim the reads for adapter/low quality bases
-	Map reads to transposon sequences
-	Create a tab-delimited count table
-	Create a presence/absence tab-delimited file
-	Visualize and cluster the binary presence/absence file

### Isolate sequences flanking established transposons
The extractRegionsOfInterest.py takes existing transposon results determined by the process described above, and extracts fasta sequences flanking the identified transposons:
```
ls -1 <all_relevant_result_files_from_a_given_transposon> | python3 extractRegionsOfInterest.py -l - -n <transposon_name> -f <unmasked_reference_fasta> -w <window_size> > <transposon_specific_intervals>.fasta
```
**Note**: any overlapping intervals are merged.
Suggested window size is 300 bps given a typical Nextera library.

### Trim reads for adapter/low quality bases
There is no specific process for this, use your favorite trimming software. Examples that we use regularly are Trimmomatic or fastp

### Map reads to transposon specific intervals
Used your favorite mapping software to map R1 reads transposon specific interval sequences. There is no need to generate R2 reads that indicate which transposon is the source.
As an example:
```
bowtie2 -x copia_specific_intervals.fasta --local --no-head -k 2 -U sample_A.fastq.gz | gzip -c > results/sample_A.sam.copia.gz
```

###	Create a tab-delimited count table
Generating a table of the results files containing the full path to a particular result file, the sample id to use for this file, and the transposon name. This should be tab-delimited and is used as input.<br>
Run the following command to convert the results files into a counts table of results for clustering and visualization:
```
perl getIntervalCountsTable.pl <results_table> > combined.presence.tsv
```

### Create a presence/absence tab-delimited file

Here is an tiny example counts table:
	| Sample A | Sample B | Sample C | Sample D
1731.2L.10186650.10187243 | 89	| 68	| 84 | 76
copia.2L.10295155.10295748 | 88 | 73 | 68 | 82
293.2L.10338613.10339206 | 18 | 0 | 11 | 1

Here is an tiny example presence table:
  | Sample A | Sample B | Sample C | Sample D
--|----------|----------|----------|---------	
1731.2L.10186650.10187243 | 1 | 1 | 1 | 1
copia.2L.10295155.10295748 | 1 | 1 | 1 | 1
293.2L.10338613.10339206 | 1 | 0 | 1 | 0

### Visualize and cluster the binary presence/absence file
See the visualization approach described above
