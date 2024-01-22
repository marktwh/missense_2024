## approach_1_inefficient_RNG

This seems like the most straightforward intuitive approach. It's the one I'm going to try. I should be able to get it to work. It's probably the least computationally efficient, but with enough computing power that might not matter so much.

### source data preparation

A list of canonical transcripts can be constructed with [this](https://www.ensembl.org/biomart/martview). The reference human transcriptome is available [here](https://ftp.ensembl.org/pub/release-111/fasta/homo_sapiens/cds/Homo_sapiens.GRCh38.cds.all.fa.gz). The AlphaMissense datasets are available [here](https://console.cloud.google.com/storage/browser/dm_alphamissense;tab=objects?prefix=&forceOnObjectsSortingFiltering=false).

**canonical transcript data**

The minimal data for processing are a .tsv with transcript identifiers and sequences. I'm using canonical transcripts so that genes with lots of transcripts aren't overrepresented in the downstream analysis. Biomart doesn't like you to download lots of sequence data at once so I used it to download a list of identifiers for canonical transcripts (as of GRCh38.p14). For the identifiers, I used ENSEMBL 'Transcript stable ID version's so that I can check that they match the transcriptome, and the AlphaMissense data later.

I downloaded the current sequences for the whole transcriptome as a FASTA. To extract the sequence data, match the transcript IDs to the canonical list, and generate the .tsv I want, I used the following R code:

```R
library(Biostrings)
library(tidyverse)
fasta <- readDNAStringSet("Homo_sapiens.GRCh38.cds.all.fa.gz")
df <- tibble(ID=names(fasta), sequence=as.character(fasta))
df2 <- as_tibble(df$ID %>% str_split(" ", simplify = TRUE))
df$ID <- df2$V1
mart <- read_csv("mart_export.csv")
mart <- mart %>% subset(select=`Transcript stable ID version`) %>% rename_at('Transcript stable ID version', ~'ID')
df <- mart %>% inner_join(df, by = "ID")
write_tsv(df, "canonical_full.tsv")
df <- df %>% filter(str_detect(sequence, "N", negate=TRUE))
write_tsv(df, "canonical_partial.tsv")
rm(list=ls())
```
I noticed that even among canonical transcripts, there are some that contain unspecified "N" nucleotides. This presents a problem because codons with Ns in can't be reliably interpreted. As a stopgap, I removed these in the 'canonical_partial.tsv'. However, a downstream workaround would be to record but not interpret these codons if they're substituted. That way, I wouldn't drop any source data unnecessarily.

I wondered if Ns have been put at the start of transcripts to maintain them in frame, so I checked whether the lengths of the transcripts were divisible by 3 without a remainder with:

```R
df <- read_tsv("canonical_full.tsv")
df$sequence_lengths <- str_length(df$sequence)
df$remainder <- df$sequence_lengths %% 3
```
They weren't. I realised I should have been working with protein-coding transcripts only, so I downloaded a filtered-version of the canonical transcript list and re-ran the script. They still weren't. Rats!

This is bad because if the number of nucleotides in a transcript is divisible by 3 without a remainder then not everything in that transcript is a codon. Worse, there may be other transcripts in the set containing non-codons whose lengths are perfectly divisible by three by chance.

I can further filter the canonical transcript list to include only transcripts with [CCDS](https://en.wikipedia.org/wiki/Consensus_CDS_Project)s. Then the canonical_full.tsv and canonical_partial.tsv files have the same number of transcripts and all of the transcript lengths are perfectly divisible by 3. These *probably* don't contain non-codons, but a lot of transcripts are dropped.

A better solution to this was to grab the list of transcripts from the AlphaMissense dataset rather than Biomart, then double-check it. Happily, it turns out the ["AlphaMissense_gene_hg38.tsv.gz"](https://console.cloud.google.com/storage/browser/_details/dm_alphamissense/AlphaMissense_gene_hg38.tsv.gz) dataset is actually a list of the transcripts they analysed rather than a list of genes.

I re-ran the script. Not all of their IDs match the IDs from the current transcriptome. They've used a previous release. I did what I should've done in the first place and read their materials and methods. They've used Gencode [hg38v32](https://www.gencodegenes.org/human/release_32.html). All of their IDs match with this release, but the CDSs in the Gencode FASTA aren't proper CDSs. The FASTA does contain coordinates for the CDSs, however. The CDSs can be extracted with:

```R
library(Biostrings)
library(tidyverse)
fasta <- readDNAStringSet("gencode.v32.pc_transcripts.fa.gz")
df <- tibble(ID=names(fasta), sequence=as.character(fasta))
df2 <- as_tibble(df$ID %>% str_split("\\|", simplify = TRUE))
df3 <- df$ID %>% str_extract("\\|CDS:\\d+\\-\\d+\\|")
df$start <- as.numeric(df3 %>% str_extract("\\|CDS:\\d+\\-") %>% str_extract("\\d+"))
df$end <- as.numeric(df3 %>% str_extract("\\-\\d+\\|") %>% str_extract("\\d+"))
df$ID <- df2$V1
df$CDS <- df$sequence %>% str_sub(df$start, df$end)
```
Among the matches to AlphaMissense, there are no transcripts containing "N"s, but there are transcripts containing non-codons/partial-codons. Balls.

The closest ENSMBL release to Gencode v32 is release 98. There's only 18 missing sequences and these would be easy to get. However, among the matches to AlphaMissense are transcripts containing "N"s and transcripts containing non-codons.

In their paper, the AlphaMissense authors state they "generated predictions for every possible single amino acid substitution within each UniProt canonical isoform" and that this approach means "every genetic missense variant can be mapped to a UniProt protein variant, regardless of the genome build". Therefore, it would be good to re-map the UniProt canonical isoform IDs to the current transcriptome. For this I need the list of 'UniProt canonical isoform IDs' from Uniref release 90, I think.

* Note: it seems like there are fewer 'UniProt canonical isoforms' than current ENSEMBL 'canonical protein-coding transcripts' (19233 vs. 23056). I'll update with a table at some point. Some of the missing ones may be in the ["AlphaMissense_isoforms_hg38.tsv.gz"](https://console.cloud.google.com/storage/browser/_details/dm_alphamissense/AlphaMissense_isoforms_hg38.tsv.gz) and ["AlphaMissense_isoforms_aa_substitutions.tsv.gz"](https://console.cloud.google.com/storage/browser/_details/dm_alphamissense/AlphaMissense_isoforms_aa_substitutions.tsv.gz) files. If they are, then it would be good to include these in the transcript list so that the source data is as complete as possible (and assemble a combined version of the prospective AlphaMissense lookup table accordingly). If they aren't, then it seems like it might be good to run translations of them through the AlphaMissense pipeline as an addition to the dataset (if one had an obscene amount of computing power at one's disposal).

To do:
1. Get Gencode hg38v32 CDSs matching gene list for 'good enough' transcript set.
2. Remap Uniprot canonical release 90 to current transcriptome for a potentially 'better' transcript set.
3. Compare and contrast.
Also:
4. See if any Uniprot noncanonical proteins in the AlphaMissense dataset are encoded by ENSEMBL canonical protein-coding transcripts. If they are, incorporate these into a 'more complete' transcript set.

Welp, this is more involved than I thought it would be. 

**AlphaMissense data**

My laptop does not have enough RAM to open either the ["AlphaMissense_hg38.tsv.gz"](https://console.cloud.google.com/storage/browser/_details/dm_alphamissense/AlphaMissense_hg38.tsv.gz) or the ["AlphaMissense_aa_substitutions.tsv.gz"](https://console.cloud.google.com/storage/browser/_details/dm_alphamissense/AlphaMissense_aa_substitutions.tsv.gz) files in R.

To do:
1. Try to open these data in chunks.
2. Spin-up an R environment in the cloud with more RAM and open them there.

### data generation

1. Approximate the probability of inosine misincorporation in place of each other nucleotide with the data from the literature above.

2. Read the sequences of transcripts and store as an appropriate type. Each transcript is a vector with a sequence of nucleotides, each with a defined position.

3. Iterate through each nucleotide of each transcript. Use a random number generator (RNG) to determine whether each nucleotide is substituted or not. For example, if there's a 1% chance of substitution and a RNG generating a number between 0 and 1, then a roll of between 0 and 0.01 means a substitution, and a roll of between 0.01 and 1 means no substitution.

4. In case of a substitution(s) in a codon, capture the codons-sequence (original and substituted), its position in the transcript, and the transcript identifier. Save these data for each run-through of the transcriptome.

### substitution effect

1. To transform the original and substituted codons into original and substituted amino acids I will initially assume that 'I' is interpreted as 'G'. This has been an assumption in the literature in the past. It means that determining original and substituted amino acids just requires [the genetic code](https://en.wikipedia.org/wiki/Genetic_code#/media/File:GeneticCode21-version-2.svg) as a lookup table.

2. The changes that replace an amino acid with a different amino acid (missense changes) will be cross-referenced with AlphaMissense. The changes that replace an amino acid with the same amino acid (synonymous changes) will be noted but not cross-referenced with AlphaMissense. The changes that introduce premature stop codons (nonsense changes) will also be noted but not cross-referenced with AlphaMissense, though they may be used in later analysis. Note that nonsense changes will often completely nullify a transcript because of the phenomenon of [nonsense-mediated decay](https://en.wikipedia.org/wiki/Nonsense-mediated_decay), but that in terminal exons they may instead cause production of truncated protein.

* In reality, 'I' may not always be interpreted as 'G'. [This](https://doi.org/10.1093%2Fnar%2Fgky1163) paper shows that. If I can get everything to work where 'I' is always interpreted as 'G', then I'll refine the code to include the possibility of alternative substitutions. The possibility of 'stalling changes' will also be considered. These possibilities will require another layer of iteration with a RNG.

### substitution functional effect

1. With the defined amino acid changes per-run, the corresponding scores in AlphaMissense can be looked up (since the dataset contains every possible missense change).

2. A first-pass metric for the 'per-transcript' effect will just be adding up the AlphaMissense scores for each transcript.

3. After enough runs, the mean and variance for the predicted functional effect on each transcript should stabilise.

