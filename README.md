# Log of the pipeline I used to filter lncRNAs from the Tropaeolum tuberosum transcriptome

## Dependencies, in order of usage:

- seqkit
- CPC2
- PlantLncBoost
- CPAT
- LncFinder
- diamond
- TransDecoder
- InterProScan
- minimap2
- samtools
- kent (UCSC bedtogenepred)
- kent (UCSC genpredtogtf)
- FEELnc

>Remember to define the file placeholders with your own paths!

## Transcriptome ncRNA filtering

The transcriptome, sequenced with poly-A capture, contains 297,918 transcripts, with isoforms (taken with `seqkit stats $file`)

```
file            format   type    num_seqs   sum_len   min_len  avg_len  max_len
*Transcriptome*  FASTA   DNA     297.918    295,386,517   184    991.5   99,287
stats: SUCCESS
```

### Filter by length

lncRNAs are generally defined as transcripts longer than 200nt, so length filtering is the first step. In our case, the initial transcripts didn't diverge much from 200

```bash
seqkit seq -m 200 "$IN_transcriptome"  > "$min200" && echo "200nt filter: SUCCESS" || echo "200nt filter: FAILED"

seqkit stats "$min200"
```

In fact, only 6 transcripts were between 184 and 200 nt in length.

### Filter by coding potential

Here, we used CPC2, PlantLncBoost, CPAT-plant and LncFinder-plant. These last two carry "plant", as they are models re-trained by the Plant-LncRNA-pipeline-v2 with the highest benchmarking metrics.

Since my `cpat` and `LncFinder` implementations are the retrained versions, I downloaded the retraining data from the `Plant-LncRNA-Pipeline` github repository clone, and set the retraining folder locally where the tools' scripts could reach it. `CPC2` and `PlantLncBoost` were more straightforward.

First, we create dedicated conda environments:

```bash
# cpc2
conda create -n cpc2_env -c cbp44 -c bioconda -c conda-forge cpc2

# CPAT
conda create -n cpat_env python=2.7 bioconda::cpat==1.2.4 -y
```
`PlantLncBoost` is not on bioconda and runs on python from its repository, no env needed
`LncFinder` runs on R locally from its repository, no env needed

I used different environment formulae in learning how to navigate the command line. `-n`names the env, `-c`sets the channel through which a program is to be found and installed by conda.


#### CPC2

```bash
#!/bin/bash

echo "===START==="

# add the SOURCE path if needed
conda activate cpc2_env

# define your files, these are generic examples

INPATH=$PATH1
OUTPATH=$PATH2
IN_transcriptome= "$INPATH/filename"
CPC2_RESULT= "$OUTPATH/filename"

set -euo pipefail
# this kills the whole job if any part fails


echo "===== Running CPC2 ====="

CPC2.py -i "$IN_transcriptome" -o "$CPC2_RESULT" && echo "cpc2: SUCCESS" || echo "cpc2: FAILED"

echo "===== DONE ====="
```

#### PlantLncBoost

```bash
#!/bin/bash

set -euo pipefail

echo "=== Inicio ==="

PLANTLNCBOOST="$PATH1"
TRANSCRIPTOME="$PATH2"
OUTDIR="$PATH3"

echo "PlantLncBoost: $PLANTLNCBOOST"
echo "Transcriptome: $TRANSCRIPTOMA"
echo "Output: $OUTDIR"

echo "=== Feature extraction ==="

python "${PLANTLNCBOOST}/Script/Feature_extraction.py" \
    -i "${TRANSCRIPTOME}" \
    -o "${OUTDIR}/PlantLncBoost_feature.csv"

echo "=== LncRNA prediction ==="

python "${PLANTLNCBOOST}/Script/PlantLncBoost_prediction.py" \
    -i "${OUTDIR}/PlantLncBoost_feature.csv" \
    -m "${PLANTLNCBOOST}/Model/PlantLncBoost_model.cb" \
    -t 0.5 \
    -o "${OUTDIR}/PlantLncBoost_prediction.csv"

echo "===count each call==="
awk -F'\t' '{print $2}' "${OUTDIR}/PlantLncBoost_prediction.csv" | sort | uniq -c

echo "=== fin  ==="

```

#### CPAT

- Default coding probability (CP) cutoff: 0.46 (CP >=0.46 indicates coding sequence, CP < 0.46 indicates noncoding sequence).

```bash
#!/bin/bash

echo "=== Inicio ==="

RETRAINING="PATH/CPAT-plant/LncRNA-retraining"
TRANSCRIPTOMA="/PATH"
OUTDIR="/Users/eduardocarriazobarajas/Desktop/TESIS/cod_pot_tools/CPAT-plant"

echo "Retraining cloned git: $RETRAINING"
echo "Transcriptoma: $TRANSCRIPTOMA"
echo "CPAT Output: $OUTDIR"

set -euo pipefail

echo "=== CPAT run ==="

cpat.py -x "$RETRAINING"/Model/Plant_Hexamer.tsv \
	-d "$RETRAINING"/Model/Plant.logit.RData \
	-g "$TRANSCRIPTOMA" \
	-o "$OUTDIR"/CPAT-plant.out \

if [ -f "$OUTDIR"/CPAT-plant.out ];
then echo "success";
else echo "CPAT output not found" >&2
fi

echo "=== fin  ==="

```


#### LncFinder

This is an `R`-based tool, so we used an Rmd

```r
---
title: "lncfinder-plant_script"
---


0. If packages were not previously installed, run this box

```R
install.packages("LncFinder")
install.packages("seqinr")


1. Load the packages:
   
```{R}
library(LncFinder)
library(seqinr)

setwd("./LncRNA-retraining")
getwd()


2. Load training data: this tells *LncFinder* how to differentiate mRNAs from noncoding transcripts

```{R}
mRNA <- seqinr::read.fasta(file ="./LncRNA-retraining/data/training/mRNA.fasta")
lncRNA <- seqinr::read.fasta(file ="./LncRNA-retraining/data/training/lncRNA.fasta")


3. features
```{R}
frequencies <- make_frequencies(cds.seq = mRNA, lncRNA.seq = lncRNA, SS.features = FALSE, cds.format = "DNA", lnc.format = "DNA", check.cds = TRUE, ignore.illegal = TRUE)


4. load plant-retrained model
```{R}
setwd("./LncRNA-retraining")
plant = readRDS("./Model/Plant_model.rda")


5. run the model on own data
```{R}
Seqs <- seqinr::read.fasta(file ="/Your_PATH_to_transcriptome")
Plant_results <- LncFinder::lnc_finder(Seqs, SS.features = FALSE, format = "DNA", frequencies.file = frequencies, svm.model = plant, parallel.cores = 2)

6. write the output table
```{R}
write.table(Plant_results, file ="cubio_plant-lncFinder.txt", sep ="\t",row.names =TRUE, col.names =TRUE,quote =FALSE)

7. inspect results:

```{R}
# df= data frame, could be named any other way
# sep= tab delimited
# header = TRUE   -> first row has column names (default, shown for clarity)
df <- read.table("cubio_plant-lncFinder.txt", sep = "\t", header = TRUE)

# names() - lists column names, so you can confirm which column holds the label
names(df)

# table() - tabulates counts of each unique value in the given vector/column
table(df$Pred)

```

#### nc ID extraction

This extracted noncoding IDs from the specific format of each tool's output. therer's a gene-only count, just got curious.

```sh
#!/bin/bash

INDIR="/path/to/cod_pot_tools"
OUTDIR="/path/to/results"

coding="$OUTDIR/coding"
noncoding="$OUTDIR/noncoding"

mkdir -p "$coding" "$noncoding"

cpat="$INDIR/CPAT-plant-copy.txt"
lncfinder="$INDIR/cubio_plant-lncFinder-copy.txt"
plncboost="$INDIR/cubio_PlantLncBoost_prediction-copy.txt"
cpc2="$INDIR/cpc2_output-copy.txt"

set -euo pipefail

echo "===Definición de archivos==="
echo "Input directory: $INDIR"
echo "Output directory: $OUTDIR"

echo "===Extracción de nc-IDs==="

echo "===CPAT-plant==="
awk -F'\t' 'NR>1 && $6 < 0.46' "$cpat" | wc -l
awk -F'\t' 'NR>1 && $6 > 0.46' "$cpat" | wc -l

awk -F'\t' 'NR>1 && $6 < 0.46 {print $1}' "$cpat" \
    > "$noncoding/cpat_noncoding_IDs.txt"
awk -F'\t' 'NR>1 && $6 > 0.46 {print $1}' "$cpat" \
    > "$coding/cpat_coding_IDs.txt"

echo "CPAT-plant coding IDs: $(wc -l < "$coding/cpat_coding_IDs.txt")"
echo "CPAT-plant noncoding IDs: $(wc -l < "$noncoding/cpat_noncoding_IDs.txt")"

echo "===PlantLncBoost==="
awk -F'\t' 'NR>1 && $2==1' "$plncboost" | wc -l
awk -F'\t' 'NR>1 && $2==0' "$plncboost" | wc -l

awk -F'\t' 'NR>1 && $2==1 {print $1}' "$plncboost" \
    > "$noncoding/plncboost_noncoding_IDs.txt"
awk -F'\t' 'NR>1 && $2==0 {print $1}' "$plncboost" \
    > "$coding/plncboost_coding_IDs.txt"

echo "PlantLncBoost coding IDs: $(wc -l < "$coding/plncboost_coding_IDs.txt")"
echo "PlantLncBoost noncoding IDs: $(wc -l < "$noncoding/plncboost_noncoding_IDs.txt")"

echo "===LncFinder-plant==="
awk -F'\t' 'NR>1 && $2=="NonCoding"' "$lncfinder" | wc -l
awk -F'\t' 'NR>1 && $2=="Coding"' "$lncfinder" | wc -l

awk -F'\t' 'NR>1 && $2=="NonCoding" {print $1}' "$lncfinder" \
    > "$noncoding/lncFinder_noncoding_IDs.txt"
awk -F'\t' 'NR>1 && $2=="Coding" {print $1}' "$lncfinder" \
    > "$coding/lncFinder_coding_IDs.txt"

echo "LncFinder-plant coding IDs: $(wc -l < "$coding/lncFinder_coding_IDs.txt")"
echo "LncFinder-plant noncoding IDs: $(wc -l < "$noncoding/lncFinder_noncoding_IDs.txt")"

echo "===CPC2==="
awk -F'\t' 'NR>1 && $8=="noncoding"' "$cpc2" | wc -l
awk -F'\t' 'NR>1 && $8=="coding"' "$cpc2" | wc -l

awk -F'\t' 'NR>1 && $8=="noncoding" {print $1}' "$cpc2" \
    > "$noncoding/cpc2_noncoding_IDs.txt"
awk -F'\t' 'NR>1 && $8=="coding" {print $1}' "$cpc2" \
    > "$coding/cpc2_coding_IDs.txt"

echo "CPC2 coding IDs: $(wc -l < "$coding/cpc2_coding_IDs.txt")"
echo "CPC2 noncoding IDs: $(wc -l < "$noncoding/cpc2_noncoding_IDs.txt")"

echo "=== c/nc gene count (no isoforms) ==="
    echo "cpat noncoding genes = $(sed -E 's/_i[0-9]+$//' "$noncoding/cpat_noncoding_IDs.txt" | sort -u | wc -l)"
    echo "cpat coding genes = $(sed -E 's/_i[0-9]+$//' "$coding/cpat_coding_IDs.txt" | sort -u | wc -l)"
    echo "====="

    echo "lncFinder noncoding genes = $(sed -E 's/_i[0-9]+$//' "$noncoding/lncFinder_noncoding_IDs.txt" | sort -u | wc -l)"
    echo "lncFinder coding genes = $(sed -E 's/_i[0-9]+$//' "$coding/lncFinder_coding_IDs.txt" | sort -u | wc -l)"
    echo "====="

    echo "plncboost noncoding genes = $(sed -E 's/_i[0-9]+$//' "$noncoding/plncboost_noncoding_IDs.txt" | sort -u | wc -l)"
    echo "plncboost coding genes = $(sed -E 's/_i[0-9]+$//' "$coding/plncboost_coding_IDs.txt" | sort -u | wc -l)"
    echo "====="

    echo "cpc2 noncoding genes = $(sed -E 's/_i[0-9]+$//' "$noncoding/cpc2_noncoding_IDs.txt" | sort -u | wc -l)"
    echo "cpc2 coding genes = $(sed -E 's/_i[0-9]+$//' "$coding/cpc2_coding_IDs.txt" | sort -u | wc -l)"

echo "=== Fin ==="
```

### Filter by protein homology



## Structural annotation with FEELnc_classifier.pl