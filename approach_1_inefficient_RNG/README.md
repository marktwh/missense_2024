## approach_1_inefficient_RNG

This seems like the most straightforward intuitive approach. It's the one I'm going to try. I should be able to get it to work. It's probably the least computationally efficient, but with enough computing power that might not matter so much.

### source data preparation

A list of canonical transcripts can be constructed with [this](https://www.ensembl.org/biomart/martview). The reference human transcriptome is available [here](https://ftp.ensembl.org/pub/release-111/fasta/homo_sapiens/cds/Homo_sapiens.GRCh38.cds.all.fa.gz). The AlphaMissense datasets are available [here](https://console.cloud.google.com/storage/browser/dm_alphamissense;tab=objects?prefix=&forceOnObjectsSortingFiltering=false).

**canonical transcript data**

The minimal data for processing are a .tsv with transcript identifiers and sequences. I'm using canonical transcripts so that genes with lots of transcripts aren't overrepresented in the downstream analysis. Biomart doesn't like you to download lots of sequence data at once so I used it to download a list of identifiers for canonical transcripts (as of GRCh38.p14). For the identifiers, I used Ensembl 'Transcript stable ID version's so that I can check that they match the transcriptome, and the AlphaMissense data.

I didn't realise this initially, but the Ensembl canonical set includes transcripts from genes on scaffolds as well as transcripts from genes on chromosomes. I excluded the scaffold-transcripts as these are probably duplicates.

Note: I am excluding transcripts from the mitochondrial genome (though there are only 13). These are a special case because mitochondria have their own transcription and translation machinery. Some mitochondrial transcript codons are translated differently as compared to the same codon in a typical mRNA. Therefore, mitochondrial transcript substitutions would have to be analysed differently. I have saved them for later in case they'd be useful.

```R
library(tidyverse)
mart <- read_tsv("mart_export_with_chr.tsv")
mart <- mart %>% select(`Transcript stable ID version`, `Chromosome/scaffold name`) %>% rename_at('Transcript stable ID version', ~'ID')
mart$scaf <- mart$`Chromosome/scaffold name` %>% str_detect("^[X|Y]|^[0-9]{1,2}", negate = TRUE)
mart <- mart %>% filter(mart$scaf == FALSE)
mart <- mart %>% subset(select=`ID`)
# 19,689 transcripts
```

I downloaded the current sequences for the whole transcriptome as a FASTA. To extract the sequence data, match the transcript IDs to the canonical list, and generate the .tsv I want, I used the following R code:

```R
library(Biostrings)
library(tidyverse)
fasta <- readDNAStringSet("Homo_sapiens.GRCh38.cds.all.fa.gz")
df <- tibble(ID=names(fasta), sequence=as.character(fasta))
df2 <- as_tibble(df$ID %>% str_split(" ", simplify = TRUE))
df$ID <- df2$V1
df <- mart %>% inner_join(df, by = "ID")
df <- df %>% distinct(.keep_all = TRUE)
write_tsv(df, "canonical_full.tsv")
# 19,689 transcripts
```

I noticed that even among canonical transcripts, there are some that contain unspecified "N" nucleotides. This presents a problem because codons with Ns in can't be reliably interpreted.

I wondered if Ns have been put at the start of transcripts to maintain them in frame, so I checked whether the lengths of the transcripts were divisible by 3 without a remainder with:

```R
df <- read_tsv("canonical_full.tsv")
df$sequence_lengths <- str_length(df$sequence)
df$remainder <- df$sequence_lengths %% 3
```
They weren't. This is bad because if the number of nucleotides in a transcript is divisible by 3 without a remainder then not everything in that transcript is a codon. Worse, there may be other transcripts in the set containing non-codons whose lengths are perfectly divisible by three by chance.

I can further filter the canonical transcript list to include only transcripts with [CCDS](https://en.wikipedia.org/wiki/Consensus_CDS_Project)s. Then the canonical_full.tsv and canonical_partial.tsv files have the same number of transcripts and all of the transcript lengths are perfectly divisible by 3. These *probably* don't contain non-codons, but a lot of transcripts are dropped.

A better solution to this was to grab the list of transcripts from the AlphaMissense dataset rather than Biomart, then double-check it. Happily, it turns out the ["AlphaMissense_gene_hg38.tsv.gz"](https://console.cloud.google.com/storage/browser/_details/dm_alphamissense/AlphaMissense_gene_hg38.tsv.gz) dataset is actually a list of the transcripts they analysed rather than a list of genes.

I re-ran the script. Not all of their IDs match the IDs from the current transcriptome. They've used a previous release, Gencode [hg38v32](https://www.gencodegenes.org/human/release_32.html). All of their IDs match with this release, but the CDSs in the Gencode FASTA aren't proper CDSs. The FASTA does contain coordinates for the CDSs, however. The CDSs can be extracted with:

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
In their paper, the AlphaMissense authors state they "generated predictions for every possible single amino acid substitution within each UniProt canonical isoform" and that this approach means "every genetic missense variant can be mapped to a UniProt protein variant, regardless of the genome build".

It would be better to use the up-to-date ENSEMBL canonical transcriptome rather than the older Gencode hg38v32 transcriptome. Therefore, I remapped the Uniprot IDs to 'Ensembl canonical' transcript IDs  with the tool available [here](https://www.uniprot.org/id-mapping). To do this, I first had to get a list of Uniprot IDs to map. These weren't in the ["AlphaMissense_gene_hg38.tsv.gz"](https://console.cloud.google.com/storage/browser/_details/dm_alphamissense/AlphaMissense_gene_hg38.tsv.gz) file but were in the larger ["AlphaMissense_hg38.tsv.gz"](https://console.cloud.google.com/storage/browser/_details/dm_alphamissense/AlphaMissense_hg38.tsv.gz) file. I initially had some difficulties opening this, but was able to do so at the cost of maxxing-out my pagefile and depositing a series of giant tempfiles on my HDD over the course of several crashes. I made a list of Uniprot IDs with:

```R
library(tidyverse)
am_canonical <- read_tsv("AlphaMissense_hg38.tsv.gz", comment = "#", col_names = FALSE, col_select = c(6,7,8,9))
am_canonical <- am_canonical %>% rename("uniprot_ID" = "X6", "ID" = "X7", "change" = "X8", "score" = "X9")
id_list <- am_canonical %>% group_by(uniprot_ID, ID) %>% group_keys()
# 19,233 transcripts
uniprot_id_list <- id_list$uniprot_ID
uniprot_id_list <- unique(uniprot_id_list)
write_tsv(as_tibble(uniprot_id_list), "uniprot_id_list.tsv")
# 19,117 Uniprot IDs
```

There are a few more transcript IDs than Uniprot IDs indicating that some distinct transcripts encode the same protein. These are probably transcripts from duplicated genes.

Remapping with the [Uniprot tool](https://www.uniprot.org/id-mapping) resulted in 49,328 hits. These probably include instances where both canonical and noncanonical transcripts from the same gene encode the same protein. 18,902 of the 19,117 Uniprot IDs mapped. 215 Uniprot IDs did not map.

I matched the mapped IDs with the Ensembl canonical IDs with:

```R
library(tidyverse)
canonical_mapped <- read_tsv("idmapping.tsv.gz")
# 18,902 of 19,117 Uniprot IDs mapped
canonical_mapped <- dplyr::rename(canonical_mapped, "uniprot_ID" = "From", "ID" = "To")
am_mapped_to_ensembl_canonical <- df %>% inner_join(canonical_mapped, by = "ID")
# 17,996 of 18,902 mapped Uniprot IDs matched to 18,120 of 19,689 Ensembl canonical IDs
```

The multiple matches for some Uniprot IDs are again probably because distinct transcripts from duplicated genes encode the same protein. The Uniprot IDs that mapped to Ensembl IDs but don't match to Ensembl canonical IDs are probably mapped to proteins encoded by noncanonical transcripts, but I should check this. If there are proteins encoded by noncanonical transcripts in the AlphaMissense "canonical" set, it follows that there may be proteins encoded by canonical transcripts in the AlphaMissense "isoforms" set.

I checked for duplicate Ensembl canonical IDs in the matched list. It turns out that a few of the Uniprot IDs have now been merged and so I removed the few entries with obsolete Uniprot IDs:

```R
library(tidyverse)
obsolete_ids <- c("Q32Q52", "Q8IXS6", "Q5VZT2", "Q9UPP5", "Q6ZW33")
am_mapped_to_ensembl_canonical <- am_mapped_to_ensembl_canonical %>% filter(!uniprot_ID %in% obsolete_ids)
# 18,120 rows for 18,120 Ensembl canonical IDs
```

I checked the 215 Uniprot IDs from AlphaMissense that did not map at all manually. The list can be generated with:

```R
library(tidyverse)
am_uniprot_unmatched <- id_list %>% anti_join(canonical_mapped, by = "uniprot_ID")
# 215 Uniprot IDs
```

It turns out that these are in three groups:

1. **121 'non-proteins'**. These are Uniprot identifiers for proteins that were thought to be real proteins but are not any more. They include artifacts, products of pseudogenes, lncRNAs, dubious gene predictions, dubious CDSs, etc. Corresponding entries will be removed from the AlphaMissense reference set to be constructed.
2. **69 'changed-proteins'**. These are Uniprot identifiers for proteins that are now thought to contain a different number of amino acids. Because the AlphaMissense scoring was done on a protein with the "wrong" number of amino acids, corresponding entries will also be removed from the AlphaMissense reference set to be constructed. However, the new Uniprot identifiers for the new proteins of the "correct" length can be used to produce a list of "real" proteins that could be re-run through the AlphaMissense pipeline to add to that dataset in future.
3. **25 'proteins that seem perfectly fine'**. These didn't map with the [Uniprot tool](https://www.uniprot.org/id-mapping) but are of the expected lengths and do seem to map to Ensembl canonical transcripts. These can be added to the `am_mapped_to_ensembl_canonical` set.

To do the addition, I used:

```R
library(tidyverse)
perfectly_fine_proteins <- read_csv("perfectly_fine_proteins.csv")
perfectly_fine_proteins <- perfectly_fine_proteins %>% rename("ID" = "new_ID", "uniprot_ID" = "am_uniprot_ID")
perfectly_fine_proteins <- df %>% inner_join(perfectly_fine_proteins, by = "ID") %>% select('ID', 'sequence', 'uniprot_ID')
am_mapped_to_ensembl_canonical <- rows_insert(am_mapped_to_ensembl_canonical, perfectly_fine_proteins)
# 18,144 transcripts; 24 added (one is only on a scaffold)
```

Checking the `am_mapped_to_ensembl_canonical` set for sequences containing "N"s showed that none did. Checking the set for transcripts with lengths not divisible by three without a remainder identified only 5 transcripts. 4 were transcripts with incomplete CDSs that are still the canonical transcripts for genes.I chose to truncate the non-coding 3' nucleotides with:

```R
library(tidyverse)
am_mapped_to_ensembl_canonical$sequence[am_mapped_to_ensembl_canonical$ID == "ENST00000506619.2"] <- str_sub(am_mapped_to_ensembl_canonical$sequence[am_mapped_to_ensembl_canonical$ID == "ENST00000506619.2"], 1, end = -3)
am_mapped_to_ensembl_canonical$sequence[am_mapped_to_ensembl_canonical$ID == "ENST00000434431.2"] <- str_sub(am_mapped_to_ensembl_canonical$sequence[am_mapped_to_ensembl_canonical$ID == "ENST00000434431.2"], 1, end = -3)
am_mapped_to_ensembl_canonical$sequence[am_mapped_to_ensembl_canonical$ID == "ENST00000621907.1"] <- str_sub(am_mapped_to_ensembl_canonical$sequence[am_mapped_to_ensembl_canonical$ID == "ENST00000621907.1"], 1, end = -2)
am_mapped_to_ensembl_canonical$sequence[am_mapped_to_ensembl_canonical$ID == "ENST00000636699.1"] <- str_sub(am_mapped_to_ensembl_canonical$sequence[am_mapped_to_ensembl_canonical$ID == "ENST00000636699.1"], 1, end = -2)
# points for a nicer way to do this
```

I also removed one transcript for a non-functional immunoglobulin:

```R
library(tidyverse)
am_mapped_to_ensembl_canonical <- am_mapped_to_ensembl_canonical %>% filter(!ID == "ENST00000569103.2")
#18,143 transcripts
```

I checked the 'start' and 'stop' codons of the remaining transcripts. Start codons are usually 'ATG' but sometimes 'CTG' or 'TTG'. Whatever the codon, they always encode methionine, but that won't matter downstream because substituted start codons will be considered differently (as possible start-loss variants). They won't be checked against AlphaMissense. 'Stop' codons are 'TGA', 'TAG' and 'TAA'.

```R
am_mapped_to_ensembl_canonical$start_codon <- am_mapped_to_ensembl_canonical$sequence %>% str_detect("^ATG|^CTG|^TTG")
am_mapped_to_ensembl_canonical$stop_codon <- am_mapped_to_ensembl_canonical$sequence %>% str_detect("TGA$|TAG$|TAA$")
no_start <- am_mapped_to_ensembl_canonical %>% filter(start_codon == FALSE)
write_tsv(no_start, "no_start.tsv")
no_stop <- am_mapped_to_ensembl_canonical %>% filter(stop_codon == FALSE)
write_tsv(no_stop, "no_stop.tsv")
```

7 transcripts had unusual start codons and 4 had no stop codon, but after checking them I'm keeping all of these in the reference set.




**AlphaMissense data**

The AlphaMissense data are of two types each split into two sets. The two types are:

1. Scores for every amino acid change that it is possible to produce with a single nucleotide change.
2. Scores for every possible amino acid change at every position.

The latter, much larger, data are preferable as they can interpret the effects of rare edge-cases such as the effects more than one simulated nucleotide substitution affecting the same codon. However, it may take much longer to search through these. Therefore, it would be good to prepare both types for flexibility.

The two sets are:

1. Scores for 'canonical' proteins. However, these contain both proteins encoded by Ensembl canonical transcripts and proteins not encoded by Ensembl canonical transcripts. The existing IDs for these are Uniprot IDs.
2. Scores for 'isoforms'. However, these contain both proteins encoded by Ensembl canonical transcripts and proteins not encoded by Ensembl canonical transcripts. The existing IDs for these are old Ensembl IDs.

To make the AlphaMissense reference set I will:

* Filter type 1&2 and set 1&2 to contain only proteins encoded by Ensembl canonical transcripts.
* Replace the previous IDs with the most recent Ensembl canonical IDs.
* Combine set 1&2 into a single set.


**Checking**

* Make sure there is a 1:1 correspondence between Ensembl IDs in the prepared transcript set and the prepared AlphaMissense set.
* Make sure that the transcripts included in the prepared transcript set encode proteins of the lengths of the proteins in the prepared AlphaMissense set.


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