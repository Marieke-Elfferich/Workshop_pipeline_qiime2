# Workshop_pipeline_qiime2
This pipeline contains the information needed to run the fungal analysis pipeline using the Bilbo server and Qiime2 with the ITSxpress plugin. It is a combination of a script from Xinming, information from the Qiime2 webpage and a sprinkle of copilot

To start off, make sure you open the qiime2 environment on Bilbo

```
micromamba activate /vol/local/conda_envs/qiime2
```

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

    echo "Processing: $r1 and $r2"
    echo "Output: $out1 and $out2"

    itsxpress --fastq "$r1" --fastq2 "$r2" \
        --region ITS2 --taxa Fungi \
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

For a small dataset, you can make this manually, but for bigger datasets, we can make this manifest on the server. Let's practice this!

In the qiime2 environment, make a new code called run_manifest.py as a python script that contains:
```python
import os
import csv
import argparse
import sys

def main():
    parser = argparse.ArgumentParser(
        description="Generate a samples.tsv file listing paired-end FASTQ files."
    )
    parser.add_argument(
        "input_dir",
        help="Path to the directory containing FASTQ files."
    )
    parser.add_argument(
        "-o", "--output",
        default="samples.tsv",
        help="Output TSV filename (default: samples.tsv)"
    )
    args = parser.parse_args()

    directory = args.input_dir

    # Verify directory exists
    if not os.path.isdir(directory):
        sys.exit(f"Error: '{directory}' is not a valid directory.")

    # Get all files in the directory
    all_files = os.listdir(directory)

    # Filter forward and reverse FASTQ files
    forward_files = sorted([f for f in all_files if f.endswith("_1.fastq.gz")])
    reverse_files = sorted([f for f in all_files if f.endswith("_2.fastq.gz")])

    # Check that the number of forward and reverse files match
    if len(forward_files) != len(reverse_files):
        sys.exit("Error: Mismatch in number of forward and reverse FASTQ files.")

    # Write to TSV
    with open(args.output, "w", newline="") as tsvfile:
        writer = csv.writer(tsvfile, delimiter="\t")
        writer.writerow(["sample-ID", "forward-absolute-filepath", "reverse-absolute-filepath"])
        for fwd, rev in zip(forward_files, reverse_files):
            # Adjust the sample ID extraction rule as needed
            parts = fwd.split("_")
            sample_id = parts[3] if len(parts) > 3 else os.path.splitext(fwd)[0]
            fwd_path = os.path.abspath(os.path.join(directory, fwd))
            rev_path = os.path.abspath(os.path.join(directory, rev))
            writer.writerow([sample_id, fwd_path, rev_path])

    print(f"{args.output} file created successfully for directory: {directory}")

if __name__ == "__main__":
    main()
```
The sample name is now decided based on the underscores in the file name. In this case, the name after the third underscore is taken as the sample ID. If your file name is different with either different deliminators or the sample ID at a different location, you will have to change the sample_id row in the code. 

Run the code by typing 
```
python run_manifest.py input_dir
```
input_dir is the directory where the files were stored after running itsxpress.
Optional: add a new name for an output file by adding -o output_file_name.tsv

Then, you can import that manifest to qiime2 by typing:
```
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
The rest of the pipeline will continue with the original output from DADA2 so don't forget to change the names of the files to table-dn-99.qza or rep-seqs-dn-99.qza when you decide to include this step in your analysis.

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

The first file is used as the basis for the OTU table and the second for the taxonomy table. However, qiime2 makes interesting ASV names such as 091cea7f43712472a172f4aa7c69a093 which could also simply be called ASV1. This is called a hash and is used by programs such as dada2 to make unique identifiers for each species, yet it makes tables very hard to read. Therefore, we change them as follows.

## Make the tables legible (example with the non-clustered tables)
First, we make a new directory called phyloseq, in which we will put our table. We will need to transform it from .qza to .biom to .tsv which can be downloaded and imported into R.
We need to transform both the table.qza and the rep-seqs.qza files. We start with the table.qza file.

Step 1: table.qza --> phyloseq/feature-table.biom
```
mkdir phyloseq
    qiime tools export \
    --input-path /vol/local/NAME/table.qza  \
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
Step 3: Modify the otu_table.tsv to clean it up a little
```
    # Use sed to delete the first line and "#OTU ID" from the
    # second line.
    cd phyloseq
    sed -i '1d' otu_table.tsv
    sed -i 's/#OTU ID//' otu_table.tsv
```
Save the following code in a python script (make a new file called table_conversion.py)
```python
import pandas as pd
from pathlib import Path

INFILE = "otu_table.tsv"
OUTFILE = "otu_table_converted.tsv"

# Load TSV, keep header row; handle odd spacing
otu = pd.read_csv(INFILE, sep="\t")

# If the first column header is blank/NaN, give it a name (it's your feature/OTU IDs)
first_col = otu.columns[0]
if pd.isna(first_col) or str(first_col).strip() == "":
    otu.rename(columns={first_col: "feature_id"}, inplace=True)
    first_col = "feature_id"

# Create ASV row labels
otu.index = [f"ASV{i+1}" for i in range(len(otu))]

# Ensure numeric sample columns are numeric
for col in otu.columns[1:]:
    otu[col] = pd.to_numeric(otu[col], errors="coerce")

# Replace the first column values with the new ASV IDs
otu.iloc[:, 0] = otu.index

# Optional: rename that first column to include 'asv' (as you intended)
otu.rename(columns={first_col: f"asv_{first_col}"}, inplace=True)

# Drop columns where all values are zero (but keep the first column)
nonzero_mask = (otu.iloc[:, 1:] != 0).any(axis=0)
filtered = pd.concat([otu.iloc[:, :1], otu.iloc[:, 1:].loc[:, nonzero_mask]], axis=1)

# Write output
filtered.to_csv(OUTFILE, sep="\t", index=False)
print(f"Wrote {OUTFILE} with shape {filtered.shape}")
```
Go back to the server and type:
```
python table_conversion.py
```

**Output artifacts:**

-   `otu_table_converted.tsv`

Download this table to use in further visualization steps in other programs (like R, you can read .tsv files with the read.delim function). 
This is your OTU-table and can be used for calculating alpha and beta diversity, or other analyses.

### The same has to be done for the rep-seqs.qza file
First locate to the folder including the rep-seqs.qza. If you've followed all the steps correctly, you should only have to type 
```
cd ..
```
Then, we follow the same step as we did for the OTU-file
```
    # Export representative sequences
    qiime tools export \
    --input-path /vol/local/NAME/rep-seqs.qza \
    --output-path phyloseq
```

**Output artifacts:**

-   `dna-sequences.fasta`

As you can see, the same step on a different dataset leads to a different data format. So we need to follow some different steps to work with this data.
```
cd phyloseq
```
```
awk '/^>/ {if (seq) print id"\t"desc"\t"seq; split(substr($0,2), a, " "); id=a[1]; desc=substr($0, index($0,$2)); seq=""} !/^>/ {seq=seq $0} END {print id"\t"desc"\t"seq}' dna-sequences.fasta > fasta_tax_table.tsv
```

Make the following python script for the tax conversion as tax_conversion.py
```python
import pandas as pd

# Input and output file names
INFILE = "fasta_tax_table.tsv"
OUTFILE = "tax_table_converted.tsv"

# Load the TSV file (assumes it has a header)
df = pd.read_csv(INFILE, sep="\t", header=None)

# Remove the second column (index 1)
df.drop(df.columns[1], axis=1, inplace=True)

# Rename the first column values to ASV1, ASV2, ...
df.iloc[:, 0] = [f"ASV{i+1}" for i in range(len(df))]

# Save the modified DataFrame to a new TSV file
df.to_csv(OUTFILE, sep="\t", index=False, header=False)

print(f"Saved converted file to {OUTFILE} with shape {df.shape}")
```
Save this and run.
```
python tax_conversion.py
```

**Output artifacts:**

- `tax_table_converted.tsv`

You can check if this has the right format. Then, we need to make this file into a qza file again. For this, we again need python. Save it as a script again and run with python script.py 
```python
# Convert a two-column TSV file (ASV ID and sequence) to FASTA format
import pandas as pd

# Load the TSV file (no header)
df = pd.read_csv("tax_table_converted.tsv", sep="\t", header=None)

# Write to FASTA format
with open("tax_table_converted.fasta", "w") as fasta:
    for i, row in df.iterrows():
        fasta.write(f">{row[0]}\n{row[1]}\n")
```
Then type:
``` 
qiime tools import \
  --input-path tax_table_converted.fasta \
  --output-path tax_table_converted.qza \
  --type 'FeatureData[Sequence]' \
  --input-format DNAFASTAFormat
```

**Output artifacts:**
-   `tax_table_converted.qza`

## Assigning taxonomy
Usually for assigning taxonomy in Qiime2, you need to train your dataset. For the Unite dataset from 2025, somebody already did that so you can download that version from [this website]([Release UNITE v10.0 2025-02-19 for qiime2-2024.10 · colinbrislawn/unite-train · GitHub](https://github.com/colinbrislawn/unite-train/releases/tag/v10.0-2025-02-19-qiime2-2024.10)) and choose which dataset you would like to use. Make sure to cite the UNITE database in your paper!
For this code, it is important that database and tax_table_converted.qza file are in the same directory.
```
qiime feature-classifier classify-sklearn \
  --i-classifier unite_ver10_dynamic_all_19.02.2025-Q2-2024.10.qza \
  --i-reads tax_table_converted.qza \
  --o-classification taxonomy.qza
```
Move the taxonomy.qza file one folder up
 ```
 mv taxonomy.qza ..
 cd ..
 ```
Then, export it to the phyloseq folder to make it the right format
```
qiime tools export \
--input-path taxonomy.qza \
--output-path phyloseq
```
Sort annotation results
```
cd phyloseq
less taxonomy.tsv|sort -k1.4n|\
awk 'BEGIN{OFS=FS="\t"}{print $1,$2 }' > tax_table.xls
```
This taxa table will put all taxonomy levels in one cell, eg, 'D_0\_\_Bacteria;D_1\_\_Firmicutes;D_2\_\_Bacilli;D_3\_\_Bacillales;D_4\_\_Bacillaceae;D_5\_\_Bacillus'. This is hard to read for R so let's split them into different cells

    cat tax_table.xls| sed 's/D_[0-9]__//g'| sed 's/\t/;/g' | sed 's/;/,/g' | sed '1s/Feature ID,Taxon/OTU,domain,phylum,class,order,family,genus,species/'> tax_table1.csv

**Output artifacts:**

- `tax_table1.csv`

And you're done! You can download the tax_table1.csv to use in further visualization or processing like further filtering of the data. 
