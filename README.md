# Workshop_pipeline_qiime2
This pipeline contains the information needed to run the fungal analysis pipeline using the Bilbo server and Qiime2 with the ITSxpress plugin. It is a combination of a script from Xinming, information from the Qiime2 webpage and a sprinkle of copilot

To start off, make sure you open the qiime2 environment on Bilbo

`micromamba activate /vol/local/conda_envs/qiime2`

Locate to the right folder where you want to do the analyses/store your samples (usually in your personal /vol/local folder)
As I always enter the server in my home folder, for me, that path is usually:
```
cd ..
cd ..
cd /vol/local/marieke
```
There, make a new directory for the new analysis, and within the directory, make two new directories to store the files in.
```
mkdir ITSxpress
cd ITSxpress
mkdir Fungal_raw
mkdir Fungal_ITS
```
Upload your raw files to the Fungal_raw directory. Normally I do this by dragging my files to the folder, but there is probably also a nerdy code way to do this. Depending on the amount of files, this can take a while.
Then, we will start with taking our samples through the [ITSxpress plugin]([GitHub - USDA-ARS-GBRU/itsxpress: Software to trim the ITS region of FASTQ sequences for amplicon sequencing analysis](https://github.com/USDA-ARS-GBRU/itsxpress)) by creating a loop that runs all files from the "Fungal_raw" data folder and outputs them into the "Fungal_ITS" folder. Make sure those names are the same in your directory, or change them in the code to fit your filenames. We will first make a script, save it and then we can run it. 

First, we make a new file with a name that finishes in .sh, telling the program that it has the format of an .sh file. Place it in the same folder as your data is located (so for example, if you have a folder called "ITSxpress" in which you have two subfolders "Fungal_raw" and "FungalITS", then place the file into the "ITSxpress" folder). You can choose to remove the logfile line from the code as it causes a lot of files to be created. If you do, remember to remove it from both parts! 

```
for r1 in Fungal_raw/*_1.fastq.gz; do
    base=$(basename "$r1" _1.fastq.gz)
    r2="Fungal_raw/${base}_2.fastq.gz"
    out1="Fungal_ITS/${base}_trimITS_1.fastq.gz"
    out2="Fungal_ITS/${base}_trimITS_2.fastq.gz"
    logfile="Fungal_ITS/${base}_logfile.txt" 
    
itsxpress --fastq "$r1" --fastq2 "$r2" \
              --region ITS2 --taxa Fungi \
              --log "$logfile" \
              --outfile "$out1" --outfile2 "$out2" \
              --threads 2
done
```

Making the code into a Linux script can be done by pasting the code into the new file (say it's called "run_itsxpress"), saving it and then typing the following into the shell:
```
chmod +x run_itsxpress.sh 
sed -i 's/\r$//' run_itsxpress.sh
bash run_itsxpress.sh
```
The second line tells it to remove any Windows-style formatting of the file, as I ran into an error with that, and the last line is to run the code. 
This step can take a while, about 30-60 seconds per sample with 2 cores assigned. You can assign more cores by increasing the number behind --threads but check through htop how much of the server is already in use and notify people beforehand through Slack if you plan on using more than 4.
I used 8 cores for 341 samples and it still took a couple of hours, but you can let it run in the background while doing other things. Do keep in mind that mobaXterm also stops if your computer goes to sleep. I would recommend the PowerToys application from Windows, which has an "awake" function that keeps your computer awake for an undefined amount of time.
## Prepare samples for dada2 pipeline

With the \_trimITS_1.fastq.gz files, you need to make a manifest that says which sample is which. The manifest should be written in notepad with tabs as spaces between the names as it will be read as a tab separated values (.tsv). It will look a bit like this (but then in .tsv format):

| sample-id   | forward-absolute-filepath                             | reverse-absolute-filepath                             |
| ----------- | ----------------------------------------------------- | ----------------------------------------------------- |
| samplename1 | /vol/local/NAME/FOLDER/samplename1_trimITS_1.fastq.gz | /vol/local/NAME/FOLDER/samplename1_trimITS_2.fastq.gz |
| samplename2 | /vol/local/NAME/FOLDER/samplename2_trimITS_1.fastq.gz | /vol/local/NAME/FOLDER/samplename2_trimITS_2.fastq.gz |
| samplename3 | /vol/local/NAME/FOLDER/samplename3_trimITS_1.fastq.gz | /vol/local/NAME/FOLDER/samplename3_trimITS_2.fastq.gz |

For a small dataset, you can make this manually, but for bigger datasets, we can make this manifest on the server. 
First, create your own python environment by running the command 
```
micromamba deactivate
micromamba create -n fastq_env python=3.10
micromamba activate fastq_env
```
And when in the environment, make a new code called run_manifest.py as a python script that contains (make sure to adjust the directories to how they fit your code and names):
```
import os
import csv

# Set the directory containing your FASTQ files

directory = "/vol/local/marieke/ITSxpress_total/Fungal_ITS"

# Get all files in the directory
all_files = os.listdir(directory)

# Filter forward and reverse FASTQ files
forward_files = sorted([f for f in all_files if f.endswith("_1.fastq.gz")])
reverse_files = sorted([f for f in all_files if f.endswith("_2.fastq.gz")])

# Check that the number of forward and reverse files match
assert len(forward_files) == len(reverse_files), "Mismatch in number of forward and reverse files"

# Write to TSV
with open("samples.tsv", "w", newline="") as tsvfile:
    writer = csv.writer(tsvfile, delimiter="\t")
    writer.writerow(["sample-ID", "forward-absolute-filepath", "reverse-absolute-filepath"])
    for i, (fwd, rev) in enumerate(zip(forward_files, reverse_files), start=1):
        sample_id = fwd.split("_")[3]
        fwd_path = os.path.abspath(os.path.join(directory, fwd))
        rev_path = os.path.abspath(os.path.join(directory, rev))
        writer.writerow([sample_id, fwd_path, rev_path])

print("samples.tsv file created successfully.")
```
Depending on the location of the samplename in your filename, sample_id might need to be adjusted. Now, it takes whatever is located after the third \_ in the name.

Run the code by typing 
```
python run_manifest.py
```

Then, you can import that manifest to qiime2 by typing:
```
micromamba deactivate
micromamba activate /vol/local/conda_envs/qiime2
qiime tools import --type 'SampleData[PairedEndSequencesWithQuality]' \
--input-path /vol/local/NAME/samples.tsv \
--output-path /vol/local/NAME/demultiplexed-seq.qza \
--input-format  PairedEndFastqManifestPhred33V2
```
You can visualize this using 
```
qiime demux summarize --i-data /vol/local/NAME/demultiplexed-seq.qza --o-visualization /vol/local/NAME/demultiplexed-seqs.qzv
```
This .qzv file can be downloaded and viewed in [Qiime2 View]([QIIME 2 View](https://view.qiime2.org/)) 

## Dada2 pipeline
Then we can run the dada2 pipeline, tweak with the parameters (of which an overview can be found on the [Qiime2 dada2 website]([denoise-paired: Denoise and dereplicate paired-end sequences — QIIME 2 2024.10.1 documentation](https://docs.qiime2.org/2024.10/plugins/available/dada2/denoise-paired/))) and visualize the track results.
I normally don't truncate but do set the max ee for both forward and reverse to 2, as is suggested in the dada2 pipeline for ITS in R.
```
qiime dada2 denoise-paired \
--i-demultiplexed-seqs /vol/local/NAME/demultiplexed-seq.qza \
--p-trunc-len-f  0 \
--p-trunc-len-r  0 \
--p-min-overlap  12 \
--p-max-ee-f 2 \
--p-max-ee-r 2 \
--p-trunc-q 2 \
--p-n-threads 2 \
--o-table /vol/local/NAME/table.qza \
--o-representative-sequences /vol/local/NAME/rep-seqs.qza \
--o-denoising-stats /vol/local/NAME/denoising-stats.qza
```
In case of very big datasets, consider using more threads to increase the speed, but always check what is available with htop AND let people know through Slack that you're planning on using many cores!! (many = more than four)

**Output artifacts:**

-   `denoising-stats.qza`
-   `table.qza`
-   `rep_seqs.qza`

The denoising-stats is similar to the track table, showing how many reads you have maintained per step, including merging and removing non-chimeric sequences. The other two files will be used later.
```
qiime metadata tabulate --m-input-file /vol/local/NAME/denoising-stats.qza --o-visualization /vol/local/NAME/denoising-stats.qzv
```
This .qzv file can be downloaded and viewed in [Qiime2 View]([QIIME 2 View](https://view.qiime2.org/)) 

## Vsearch cluster
This step can be used to cluster ASVs together into OTUs yet there has been some debate as to whether ASVs or OTUs are better to use. For bacteria, the consensus seems to be ASVs but for fungi, it could be good to use OTUs depending on your question. You can skip this step if you want to work with ASVs instead of OTUs.

```
    qiime vsearch cluster-features-de-novo \
    --i-table /vol/local/NAME/table.qza \
    --i-sequences /vol/local/NAME/rep-seqs.qza \
    --p-perc-identity 0.99 \
    --o-clustered-table /vol/local/NAME/table-dn-99.qza \
    --o-clustered-sequences /vol/local/NAME/rep-seqs-dn-99.qza \
    --p-threads 0
```

**Output artifacts:**

-   `table-dn-99.qza`
-   `rep-seqs-dn-99.qza`

The first file is used as the basis for the OTU table and the second for the taxonomy table. However, qiime2 makes interesting ASV names such as 091cea7f43712472a172f4aa7c69a093 which could also simply be called ASV1. To do this, you have to take a bit of a scenic route through R.

## Scenic route through R to make the tables legible
First, we make a new directory called phyloseq, in which we will put our table. We will need to transform it from .qza to .biom to .tsv which can be downloaded and imported into R.

Step 1: table-dn-99.qza --> phyloseq/feature-table.biom
```
mkdir phyloseq
    qiime tools export \
    --input-path /vol/local/NAME/table-dn-99.qza  \
    --output-path phyloseq
```
Step 2: phyloseq/feature-table.biom --> phyloseq/otu_table.tsv
```
    # Convert biom format to tsv format
    biom convert \
    -i phyloseq/feature-table.biom \
    -o phyloseq/otu_table.tsv \
    --to-tsv
```
Step 3: Modify the otu_table.tsv to make it easier to read into R
```
    # Use sed to delete the first line and "#OTU ID" from the
    # second line.
    cd phyloseq
    sed -i '1d' otu_table.tsv
    sed -i 's/#OTU ID//' otu_table.tsv
    cd ../
```
Upload this file into R
``` r
library(pacman)
pacman::p_load(tidyverse,magrittr,stringr)
setwd("C:\\Users\\NAME\\OneDrive - Universiteit Leiden\\Documenten\\FOLDER_WITH_OTU_TABLE")
otu <- "otu_table.tsv" %>%
  read.delim(check.names = FALSE,header = T,sep="\t")

rown <- paste0("ASV",seq_len(nrow(otu)))
otu[,1] <- rown
colnames(otu)[1] <- paste0("asv",colnames(data)[1])
filtered_table <- otu[, colSums(otu != 0) > 0] #this step is to guarantee that all the samples contain ASVs
write.table (filtered_table,file ="otu-table-2.tsv", sep ="\t", row.names = F)   
```
Now, we have changed the content of `otu_table.tsv` and cleaned it from any empty samples, and then we will restore it back as follows:

    otu_table.tsv --> feature-table.biom ---> otu_table.qza

Notice, we should store the otu_table_2.tsv in the folder phyloseq.

```
    cd phyloseq
```
```
    biom convert -i otu_table_2.tsv -o feature-table.biom --to-hdf5 --table-type="OTU table"
```
```
    qiime tools import \
      --input-path feature-table.biom \
      --type 'FeatureTable[Frequency]' \
      --input-format BIOMV210Format \
      --output-path otu-table-rename.qza
```

**Output artifacts:**

-   `otu_table-rename.qza`

### The same has to be done for the rep-seqs-dn-99.qza file (or the rep-seqs when you're working with ASV's instead of OTU's)

```
    # Export representative sequences
    qiime tools export \
    --input-path /vol/local/marieke/rep-seqs-dn-99.qza \
    --output-path phyloseq
```

**Output artifacts:**

-   `dna-sequences.fasta`

```{=html}
<!-- -->
```
    less dna-sequences.fasta |paste - -|sed '1i ASVID,seq' > rep.fa

-   `rep.fa`
Download the rep.fa file and upload it into R
``` r
library(pacman)
pacman::p_load(tidyverse,magrittr,stringr)
rep <- "rep.fa" %>%
read.delim(check.names = FALSE, row.names = 1) %>%
set_rownames(paste0(">ASV", seq_len(nrow(.))))
write.table (rep,file ="rep.xls", sep ="\t", row.names = T,quote = F)
```

-   `rep.xls`
Upload this file into the general directory
```{=html}
<!-- -->
```
    less rep.xls|sed '1d'|sed 's/"//g'|sed 's/\r//g'|tr "\t" "\n" > rep-seqs.fasta
    cd ..

-   `rep-seqs.fasta`

``` go
qiime tools import \
--type 'FeatureData[Sequence]' \
--input-path phyloseq/rep-seqs.fasta \
--output-path phyloseq/rep-seqs-rename.qza
```

-   `rep-seqs-rename.qza`

## Assigning taxonomy
Usually for assigning taxonomy in Qiime2, you need to train your dataset. For the Unite dataset from 2025, somebody already did that so you can download that version from [this website]([Release UNITE v10.0 2025-02-19 for qiime2-2024.10 · colinbrislawn/unite-train · GitHub](https://github.com/colinbrislawn/unite-train/releases/tag/v10.0-2025-02-19-qiime2-2024.10)) and choose which dataset you would like to use. Make sure to cite the UNITE database in your paper!
For this code, it is important that database and rep-seqs-rename.qza file are in the same directory.
```
qiime feature-classifier classify-sklearn \
  --i-classifier unite_ver10_dynamic_all_19.02.2025-Q2-2024.10.qza \
  --i-reads rep-seqs-rename.qza \
  --o-classification taxonomy.qza
```


Export it to the phyloseq folder
```
qiime tools export \
--input-path taxonomy.qza \
--output-path phyloseq
```
Sort annotation results
```
cd phyloseq
less taxonomy.tsv|sort -k1.4n|\
awk 'BEGIN{OFS=FS="\t"}{print $1,$2 }' > otu_taxa.xls
```
This taxa table will put all taxonomy levels in one cell, eg, 'D_0\_\_Bacteria;D_1\_\_Firmicutes;D_2\_\_Bacilli;D_3\_\_Bacillales;D_4\_\_Bacillaceae;D_5\_\_Bacillus'. This is hard to read for R so let's split them into different cells

    cat otu_taxa.xls| sed 's/D_[0-9]__//g'| sed 's/\t/;/g' | sed 's/;/,/g' | sed '1s/Feature ID,Taxon/OTU,domain,phylum,class,order,family,genus,species/'> otu_taxa1.csv

And you're done! Further visualization and analysis can all be done in R or another software of your preference. 
