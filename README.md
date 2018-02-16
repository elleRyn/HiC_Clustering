This is a workflow for using Hi-C paired end read linkages to cluster metagenomic contigs by 
species. 

Project website and data repository: [Using Hi-C to Track Plasmids and Antibiotic
Resistance in Microbial Communities](https://osf.io/gr2d7/)

By Elle J. Kohler, kohl5779@vandals.uidaho.edu

Project paper: link and citation coming soon

## Overview
Prior to performing contig clustering based on Hi-C read linkages, metagenomic shotgun 
reads must be cleaned and assembled into contigs and the Hi-C reads must be cleaned. See 
[project website](https://osf.io/7n8rx/wiki/home/) for how to complete these first two steps. 
Tutorials for these steps are given for microbial [Community 4a](https://osf.io/zacf7/) 
and a tutorial for the third step of the analysis pipeline, clustering metagenomic contigs 
based on Hi-C read linkages, is given here for the same sample, Community 4a. 
Step 2, the Hi-C read cleaning, results in an alignment file called: **hicup.sam**, that 
gives alignments specifying which contig each end of each Hi-C read aligns to. This is 
the data file used for the Hi-C cluster analysis that is done in R.

## Data prep
The file: [hicup.sam](https://osf.io/prnvy/), is too large to easily maneuver in R and 
only a small percentage of its information is actually needed as input. For this project 
the important data was extracted using the script available on this page, **Matrix.py**. 
**Matrix.py** should be downloaded, made an executable script, and implemented in the 
command line using the command:

```> python2.7 Matrix.py hicup.sam > Comm4a_interactions.csv```

This extracts column 3, which lists the contigs the Hi-C reads aligned to, of **hicup.sam** 
in a short format. Alternatively this third column could be extracted any way you choose.

## R: build Hi-C interaction matrix
Download the script **HiCclust.R**. The first step in running it will be to install the 3 
packages listed. Next there are 2 variables to set up. First one must specify the number 
of valid and unique Hi-C reads in the dataset. Using output from the proceeding portions 
of this analysis pipeline, the number of valid and unique Hi-C reads
can be found at the bottom of **AD004_S2_R1_2_001.HiCUP_summary_report.html**. For 
Community 4a, this number is: 5911667. The next piece of information that needs filled in is 
the path to the file **Comm4a_interactions.csv** (where it is stored on your computer). 

Now you can run the section of script called "build symmetrical matrix of Hi-C interactions". 
If you're working in Rstudio, you'll see a lot of objects added to your environment. Clicking 
on the one called: **imat**, shows the matrix of contigs with the numbers of Hi-C reads 
linking each contig pair.

## R: perform hierarchical clustering
Next run the 3 lines under the heading "perform clustering". This is done using the R 
[hclust](https://stat.ethz.ch/R-manual/R-devel/library/stats/html/hclust.html) command 
and will produce a dendrogram that looks like this: ![Community 4a dendrogram](dend.png)

For this example and all other samples run in our study, the dendrogram clearly shows the 
likely number of species clusters in the sample. Here there are 4 clear branches on the 
tree. This number must be manually input into the next section of code, "input k". There is 
also the option to save a file to your Desktop listing which cluster number each contig 
got assigned to. It is described later how to use this file, called **clusters.txt** to 
separate the contigs into a separate file for each cluster to be used in downstream analysis.

Alternatively, scroll down to the end of **HiCclust.R** (lines 97-102) and see the section 
on using silhouette plots to determine the optimal number of clusters. This can be tested 
on Community 4a to verify that the optimal number of clusters is indeed 4. This or similar 
methods may prove useful for determining the number of species in more complex samples 
where the optimal number of clusters is not visually obvious from the dendrogram.

## R: species identification
The next section of the script (lines 45-62 of **HiCclust.R**) computes some statistics on 
the clusters, such as which reference genomes the contigs in each cluster align to best and 
how many contigs in the dataset did not align to any of the provided reference genomes. It 
also creates a dataset that can be used to color the data points (contigs) in the dendrogram 
or following principle coordinate (PCo) plots by species as determined by alignment to 
reference genomes. This portion of the code requires another input file made available here: 
**outfile1**.

If you wish to reproduce the file yourself, **outfile1** was created by aligning all of 
the contigs to a compilation of reference genomes 
containing all the replicons known to be in the sample. The species and plasmids present in 
sample Community 4a are listed [here](https://osf.io/zacf7/wiki/home/). The reference 
genomes for these replicons can be collected [here](https://osf.io/xv3ud/) or pulled from 
the NCBI reference database. The required references should be concatenated into a single 
file called Comm4aref.fasta. The NCBI blast tool can then be used to perform the alignment. 
First turn the contigs into a database:

```> makeblastdb -dbtype nucl -in longcontigs.fasta -hash_index```

Then blast the file of combined reference sequences against the assembled contigs:

```> blastn -query Comm4aref.fasta -db longcontigs.fasta -outfmt 6 -evalue 0.0001 -out outfile -max_hsps 1```

The output from this command will be a file that again contains more information than is 
needed as input into R. The necessary information for script **HiCclust.r** being only 
columns 1, 2, and 4 of the **outfile**. These columns can be extracted using the command:

```> awk '{print$1,$2,$4}' outfile > outfile1```

**outfile1** can then be used as the input for this portion of the R script. Remember to 
specify the path to this file.

## R: visualization in PCo space
There is next the option to look at the clusters in 2D or 3D principle coordinate (PCo) space.
In this script the default is to color the contigs by species identity as was determined
by aligning the contigs to reference genomes. I usually jump right to viewing the clusters
in 3D PCo space as this visual is more informative (lines 71-77). If working on a Mac you may 
need to install [XQuartz](https://www.xquartz.org/) before being able to view 3D plots. 
XQuartz must already be open and running in the background any time this portion of the script is run (lines 71-77). 
After manually setting the plot to the desired size and direction there are the options 
to add a legend or take a snapshot of the 3D figure producing something like this: 
![Community 4a in 3D PCo space](Comm4a3D.png)

The following section of code (lines 89-94) show how to manually set the colors to something 
different. If you supply a vector of desired colors, you can then scroll back up to any of 
the plot commands and change the line "col = as.numeric(dataset$V1)+**some number**" to 
simply "col = as.numeric(dataset$V1)" and the plot will use the colors you specified.

## Downstream cluster analysis
If you previously saved the file, **clusters.txt**, that lists which cluster each contig 
was assigned to based on hierarchical cluster, you can use it and the script, **Clusters.py**, 
to separate the file of total contigs, **longcontigs.fasta** into multiple smaller files of 
contigs, one for each cluster. This can be done using the command:

```> python2.7 Clusters.py clusters.txt longcontigs.fasta```

This will automatically create and name 4 files that will be output into whatever directory 
you are currently working in. Note that this command will automatically create as many 
contig files as there were clusters (the value chosen for k). These files each represent a 
genome as determined by the hierarchical clustering based on Hi-C linkages. They can be used 
to align to reference sequences in order to determine how much of each species genome was 
recovered or used for blast type analysis methods.

An example of some blast searches that I used for tracking antibiotic resistance genes 
(ARG) and plasmid identifier genes within the clusters were the online search tools: 
[ResFinder](https://cge.cbs.dtu.dk/services/ResFinder/) and 
[PlasmidFinder](https://cge.cbs.dtu.dk/services/PlasmidFinder/). Either the file of total 
contigs, **longcontigs.fasta**, or the individual genome cluster files can be uploaded to 
these sites and compared to their databases of ARG and plasmid identifier genes to get a 
quick look at who may be carrying what. I recommend

