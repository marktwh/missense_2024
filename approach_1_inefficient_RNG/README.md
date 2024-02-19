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
library(tidyverse)
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
library(tidyverse)
am_mapped_to_ensembl_canonical$start_codon <- am_mapped_to_ensembl_canonical$sequence %>% str_detect("^ATG|^CTG|^TTG")
am_mapped_to_ensembl_canonical$stop_codon <- am_mapped_to_ensembl_canonical$sequence %>% str_detect("TGA$|TAG$|TAA$")
no_start <- am_mapped_to_ensembl_canonical %>% filter(start_codon == FALSE)
write_tsv(no_start, "no_start.tsv")
no_stop <- am_mapped_to_ensembl_canonical %>% filter(stop_codon == FALSE)
write_tsv(no_stop, "no_stop.tsv")
#18,143 transcripts
```

7 transcripts had unusual start codons and 4 had no stop codon, but after checking them I'm keeping all of these in the reference set.

There are 1,546 Ensembl canonical transcripts that the Uniprot IDs from the AlphaMissense 'canonical' set did not match to. It was possible some of these were in the AlphaMissense 'isoforms' sets. Unfortunately, the 'isoforms' sets have Ensembl IDs rather than Uniprot IDs. To fully open the files and generate lists of these IDs, it was necessary to use a VM in the cloud with 64GB of RAM. Strangely, the "AlphaMissense_isoforms_hg38.tsv.gz" file contained data for 63,099 while the "AlphaMissense_isoforms_aa_substitutions.tsv.gz" file contained data for an additional 13.

```R
library(tidyverse)
am_isoforms <- read_tsv("AlphaMissense_isoforms_hg38.tsv.gz", comment = "#", col_select = c(6), col_names = FALSE)
am_isoforms <- unique(am_isoforms)
am_isoforms <- rename(am_isoforms, "ID" = "X6")
write_tsv(am_isoforms, "am_isoforms_list.tsv")
am_isoforms_full <- read_tsv("AlphaMissense_isoforms_aa_substitutions.tsv.gz", comment = "#", col_select = c(1))
am_isoforms_full <- unique(am_isoforms_full)
am_isoforms <- rename(am_isoforms_full, "ID" = "transcript_id")
write_tsv(am_isoforms_full, "am_isoforms_full_list.tsv")
#VM with 64GB RAM
#63,099 transcripts in "AlphaMissense_isoforms_hg38.tsv.gz"
#63,112 transcripts in "AlphaMissense_isoforms_aa_substitutions.tsv.gz"
```

To match the Ensembl IDs from the Ensembl 'canonical' set to the IDs from the 'isoforms' sets, I converted the 'transcript stable IDs with version' to just 'transcript stable IDs'. I matched with these.

This yielded 1,224 matches. The matches included some very long transcripts encoding very large proteins. These may be particularly relevant in the planned analyses.

```R
library(tidyverse)
canonical_missing <- df %>% filter(!ID %in% am_mapped_to_ensembl_canonical$ID)
canonical_missing$short_ID <- canonical_missing$ID %>% str_extract("ENST\\d+")
am_isoforms_full_list$short_ID <- am_isoforms_full_list$ID %>% str_extract("ENST\\d+")
canonical_missing_matched <- canonical_missing %>% inner_join(am_isoforms_full_list, by = "short_ID")
canonical_missing_matched <- canonical_missing_matched %>% rename("ID" = "ID.x", "am_ID" = "ID.y")
# 1,224 transcripts
```

None of these transcripts contained "N"s. However, 86 contained non-codons and/or lacked a conventional start codon and/or lacked a conventional start codon. I checked these manually. A high proportion were incomplete 'novel transcript's encoding 'novel protein'. Nevertheless, I'll keep them in the reference set. I truncated non-coding 3' nucleotides where appropriate with:

```R
library(tidyverse)
clip_1 <- read_csv("clip_1.csv")
clip_2 <- read_csv("clip_2.csv")
canonical_missing_matched <- canonical_missing_matched %>%
  mutate(sequence = ifelse(ID %in% clip_1$ID, str_sub(sequence, 1, end = -2), sequence))
canonical_missing_matched <- canonical_missing_matched %>%
  mutate(sequence = ifelse(ID %in% clip_2$ID, str_sub(sequence, 1, end = -3), sequence))
# thanks to Pierfrancesco Butti
```

Unexpectedly, whilst checking, I found the data from the "AlphaMissense_isoforms_aa_substitutions.tsv.gz" file contained data for multiple versions of the same transcript in some cases. I removed the duplicates manually leaving the most recent version numbers.

```R
dupes <- read_csv("dupes.csv")
canonical_missing_matched <- canonical_missing_matched %>% filter(!am_ID %in% dupes$ID)
# 1,166 transcripts
```

The "AlphaMissense_hg38.tsv.gz" and "AlphaMissense_aa_substitutions" datasets were purported to contain data on the exact-same proteins. Just with more amino acid changes scored in the "AlphaMissense_aa_substitutions" set. However, as I had found that there were fewer proteins in the "AlphaMissense_isoforms_hg38.tsv.gz" set than the "AlphaMissense_isoforms_aa_substitutions.tsv.gz" set, I checked this:

```R
library(tidyverse)
am_canonical_full <- read_tsv("AlphaMissense_aa_substitutions.tsv.gz", comment = "#", col_select = c(1))
am_canonical_full <- unique(am_canonical_full)
am_canonical <- read_tsv("AlphaMissense_hg38.tsv.gz", comment = "#", col_names = FALSE, col_select = c(6,7,8,9))
am_canonical <- am_canonical %>% rename("uniprot_ID" = "X6", "ID" = "X7", "change" = "X8", "score" = "X9")
id_list <- am_canonical %>% group_by(uniprot_ID, ID) %>% group_keys()
am_canonical_full <- am_canonical_full %>% rename("uniprot_ID" = "ID")
more_missing <- am_canonical_full %>% anti_join(id_list, by = "uniprot_ID")
write_tsv (more_missing, "more_missing.tsv")
# 1,399 additional Uniprot IDs
```

There were 1,399 additional Uniprot IDs in the "AlphaMissense_aa_substitutions" set, 1,222 of which mapped to Ensembl transcripts. However, only 282 of these mapped to Ensembl canonical transcripts. Of that 282, 278 were among those canonical transcripts that had not been matched to anything yet. The small overlap was due to obsolete IDs.

```R
library(tidyverse)
more_missing <- read_tsv("more_missing.tsv")
more_missing_mapped <- read_tsv("extra_idmapping_2024_02_14.tsv")
more_missing_mapped <- dplyr::rename(more_missing_mapped, "uniprot_ID" = "From", "ID" = "To")
more_missing_mapped_to_ensembl_canonical <- df %>% inner_join(more_missing_mapped, by = "ID")
ensembl_canonical_unmatched <- df %>% anti_join(am_mapped_to_ensembl_canonical, by = "ID")
ensembl_canonical_unmatched <- unique(ensembl_canonical_unmatched)
more_missing_mapped_to_ensembl_unmatched <- ensembl_canonical_unmatched %>% inner_join(more_missing_mapped, by = "ID")
discrepancy <- more_missing_mapped_to_ensembl_canonical %>% anti_join(more_missing_mapped_to_ensembl_unmatched, by = "ID")
rm1 <- c("ENST00000410005.2", "ENST00000291552.9")
am_mapped_to_ensembl_canonical <- am_mapped_to_ensembl_canonical %>% filter(!ID %in% rm1)
rm2 <- c("ENST00000450895.2", "ENST00000291554.6")
more_missing_mapped_to_ensembl_canonical <- more_missing_mapped_to_ensembl_canonical %>% filter(!ID %in% rm2)
# 280 transcripts
```

Of these, one transcript contained "N"s and was removed. One had no start codon but was kept.

```R
library(tidyverse)
more_missing_mapped_to_ensembl_canonical$sequence_lengths <- str_length(more_missing_mapped_to_ensembl_canonical$sequence)
more_missing_mapped_to_ensembl_canonical$remainder <- more_missing_mapped_to_ensembl_canonical$sequence_lengths %% 3
more_missing_mapped_to_ensembl_canonical$start_codon <- more_missing_mapped_to_ensembl_canonical$sequence %>% str_detect("^ATG|^CTG|^TTG")
more_missing_mapped_to_ensembl_canonical$stop_codon <- more_missing_mapped_to_ensembl_canonical$sequence %>% str_detect("TGA$|TAG$|TAA$")
no_start <- more_missing_mapped_to_ensembl_canonical %>% filter(start_codon == FALSE)
no_stop <- more_missing_mapped_to_ensembl_canonical %>% filter(stop_codon == FALSE)
more_missing_mapped_to_ensembl_canonical <- more_missing_mapped_to_ensembl_canonical %>% filter(!ID == "ENST00000477973.5")
# 279 transcripts
```

There was substantial overlap between the 'more_missing_mapped_to_ensembl_canonical' set and the 'canonical_missing_matched' sets indicating that there's some overlap in the content of the original AlphaMissense 'canonical' and 'isoforms' sets. However, I won't remove the overlapping transcripts from either 'more_missing_mapped_to_ensembl_canonical' or 'canonical_missing_matched'. These will be removed when the data are merged, and similarly removed when the AlphaMissense reference datasets are created.

```R
library(tidyverse)
no_overlap <- more_missing_mapped_to_ensembl_canonical %>% filter(!ID %in% canonical_missing_matched$ID)
# 102 of 279 transcripts
```

I combined the extra transcripts with 'am_mapped_to_ensembl_canonical':

```R
library(tidyverse)
am_mapped_to_ensembl_canonical <- am_mapped_to_ensembl_canonical %>% select(`ID`, `sequence`, `uniprot_ID`)
more_missing_mapped_to_ensembl_canonical <- more_missing_mapped_to_ensembl_canonical %>% select(`ID`, `sequence`, `uniprot_ID`)
am_mapped_to_ensembl_canonical <- rows_insert(am_mapped_to_ensembl_canonical, more_missing_mapped_to_ensembl_canonical)
# 18,420 transcripts
```

So that's 18,420 transcripts to match with the AlphaMissense 'canonical' sets and 1,166 to match with the 'isoforms' sets, though there's some overlap. I made the reference set with:

```R
library(tidyverse)
am_mapped_to_ensembl_canonical <- am_mapped_to_ensembl_canonical %>% select(`ID`, `sequence`)
canonical_missing_matched <- canonical_missing_matched %>% select(`ID`, `sequence`)
reference_transcripts_full <- rows_insert(am_mapped_to_ensembl_canonical, canonical_missing_matched, conflict = "ignore")
reference_transcripts_full <- unique(reference_transcripts_full)
write_tsv(reference_transcripts_full, "reference_transcripts_full.tsv")
# 19,409 transcripts
```

**AlphaMissense data**

The AlphaMissense data are of two types each split into two sets. The two types are:

1. Scores for every amino acid change that it is possible to produce with a single nucleotide change.
2. Scores for every possible amino acid change at every position.

The latter, much larger, data are preferable as they can interpret the effects of rare edge-cases such as the effects more than one simulated nucleotide substitution affecting the same codon. However, it may take much longer to search through these. Therefore, it would be good to prepare both types for flexibility. Note: from preparing the 

The two sets are:

1. Scores for 'canonical' proteins. However, these contain both proteins encoded by Ensembl canonical transcripts and proteins not encoded by Ensembl canonical transcripts. The existing IDs for these are Uniprot IDs.
2. Scores for 'isoforms'. However, these contain both proteins encoded by Ensembl canonical transcripts and proteins not encoded by Ensembl canonical transcripts. The existing IDs for these are old Ensembl IDs.

To make the AlphaMissense reference set I:

* Filter type 1&2 and set 1&2 to contain only proteins encoded by Ensembl canonical transcripts.

```R
library(tidyverse)
am_canonical_full <- read_tsv("AlphaMissense_aa_substitutions.tsv.gz", comment = "#")
am_canonical_full <- unique(am_canonical_full)
am_canonical_full <- am_canonical_full %>% rename("uniprot_ID" = "uniprot_id", "change" = "protein_variant", "score" = "am_pathogenicity")
am_mapped_to_ensembl_canonical <- read_tsv("am_mapped_to_ensembl_canonical.tsv")
am_mapped_to_ensembl_canonical <- am_mapped_to_ensembl_canonical %>% select(`ID`, `uniprot_ID`)
am_canonical_full <- am_mapped_to_ensembl_canonical %>% left_join(am_canonical_full, by = "uniprot_ID")
# "AlphaMissense_aa_substitutions.tsv.gz"
am_canonical <- read_tsv("AlphaMissense_hg38.tsv.gz", comment = "#", col_names = FALSE, col_select = c(6,8,9))
am_canonical <- unique(am_canonical)
am_canonical <- am_canonical %>% rename("uniprot_ID" = "X6", "change" = "X8", "score" = "X9")
am_mapped_to_ensembl_canonical <- am_mapped_to_ensembl_canonical %>% select(`ID`, `uniprot_ID`)
am_canonical <- am_mapped_to_ensembl_canonical %>% left_join(am_canonical, by = "uniprot_ID")
#"AlphaMissense_hg38.tsv.gz"
am_isoforms_full <- read_tsv("AlphaMissense_isoforms_aa_substitutions.tsv.gz", comment = "#")
canonical_missing_matched <- canonical_missing_matched %>% select(`ID`, `am_ID`)
am_isoforms_full <- am_isoforms_full %>% rename("am_ID" = "transcript_id", "change" = "protein_variant", "score" = "am_pathogenicity")
am_isoforms_full <- canonical_missing_matched %>% left_join(am_isoforms_full, by = "am_ID")
# "AlphaMissense_isoforms_aa_substitutions.tsv.gz"
am_isoforms <- read_tsv("AlphaMissense_isoforms_hg38.tsv.gz", comment = "#", col_names = FALSE, col_select = c(6, 7, 8))
am_isoforms <- am_isoforms %>% rename("am_ID" = "X6", "change" = "X7", "score" = "X8")
am_isoforms <- canonical_missing_matched %>% left_join(am_isoforms, by = "am_ID")
# "AlphaMissense_isoforms_hg38.tsv.gz"
# VM with 128GB RAM
```

* Replace the previous IDs with the most recent Ensembl canonical IDs.

```R
library(tidyverse)
am_canonical_full <- am_canonical_full %>% select(`ID`, `change`, `score`)
am_canonical_full <- unique(am_canonical_full)
am_canonical_full <- am_canonical_full %>% drop_na(score)
write_tsv(am_canonical_full, "am_canonical_full.tsv.gz")
# "AlphaMissense_aa_substitutions.tsv.gz"
am_canonical <- am_canonical %>% select(`ID`, `change`, `score`)
am_canonical <- unique(am_canonical)
am_canonical <- am_canonical %>% drop_na(score)
write_tsv(am_canonical, "am_canonical.tsv.gz")
#"AlphaMissense_hg38.tsv.gz"
am_isoforms_full <- am_isoforms_full %>% select(`ID`, `change`, `score`)
am_isoforms_full <- unique(am_isoforms_full)
am_isoforms_full <- am_isoforms_full %>% drop_na(score)
write_tsv(am_isoforms_full, "am_isoforms_full.tsv.gz")
# "AlphaMissense_isoforms_aa_substitutions.tsv.gz"
am_isoforms <- am_isoforms %>% select(`ID`, `change`, `score`)
am_isoforms <- unique(am_isoforms)
am_isoforms <- am_isoforms %>% drop_na(score)
write_tsv(am_isoforms, "am_isoforms.tsv.gz")
# "AlphaMissense_isoforms_hg38.tsv.gz"
# VM with 128GB RAM
```

* Combine set 1&2 into a single set.

```R
library(tidyverse)
am_isoforms_full <- read_tsv("am_isoforms_full.tsv.gz")
am_full <- rows_insert(am_canonical_full, am_isoforms_full, conflict = "ignore")
write_tsv(am_full, "am_full.tsv.gz")
# full set (2)
am_canonical <- read_tsv("am_canonical.tsv.gz")
am_isoforms <- read_tsv("am_isoforms.tsv.gz")
am <- rows_insert(am_canonical, am_isoforms, conflict = "ignore")
write_tsv(am, "am.tsv.gz")
# partial set (1)
# VM with 128GB RAM
```

* am contains 19,302 transcripts and 65,052,376 changes
* am_full contains 19,409 transcripts and 212,067,740 changes


**Checking**

* Make sure there is a 1:1 correspondence between Ensembl IDs in the prepared transcript set and the prepared AlphaMissense set.
* Make sure that the transcripts included in the prepared transcript set encode proteins of the lengths of the proteins in the prepared AlphaMissense set.

In order to extract amino acid lengths of proteins and their sequences from the AlphaMissense data I used the following. Note here that amino acid lengths are given by 'max' amino acid numbers rather than by counting amino acids per protein. This is because come of the proteins start at amino acid position '1' and some at position '2'.

```R
library(tidyverse)
am_full <- read_tsv("am_full.tsv.gz")
am_full$aa_no <- am_full$change %>% str_extract("\\d+")
am_full$ref_aa <- am_full$change %>% str_extract("[^0-9]+")
am_full <- am_full %>% select(`ID`, `aa_no`, `ref_aa`)
am_full <- unique(am_full)
am_full <- am_full %>% transform(aa_no = as.numeric(am_full$aa_no))
am_full_protein_length <- am_full %>%
  group_by(ID) %>%
  filter(aa_no == max(aa_no)) %>%
  select(`ID`, `aa_no`) %>%
  rename("amino_acids" = "aa_no")
write_tsv(am_full_protein_length, "am_full_protein_length.tsv")
am_full <- arrange_all(am_full)
am_full_protein_sequence <- am_full %>%
  group_by(ID) %>%
  mutate(ref_aa = paste0(ref_aa, collapse = "")) %>%
  rename("aa_sequence" = "ref_aa") %>%
  select(`ID`, `aa_sequence`)
am_full_protein_sequence <- unique(am_full_protein_sequence)
am_full_protein_details <- am_full_protein_sequence %>% inner_join(am_full_protein_length, by = "ID")
am_full_protein_details <- unique(am_full_protein_details)
am_full_protein_details$seq_length <- str_length(am_full_protein_details$aa_sequence)
am_full_protein_details <- am_full_protein_details %>%
  mutate(comparison = ifelse(amino_acids == seq_length, TRUE, FALSE))
write_tsv(am_full_protein_details, "am_full_protein_details.tsv")
#am_full_protein_details; 19,409 proteins
am <- read_tsv("am.tsv.gz")
am$aa_no <- am$change %>% str_extract("\\d+")
am$ref_aa <- am$change %>% str_extract("[^0-9]+")
am <- am %>% select(`ID`, `aa_no`, `ref_aa`)
am <- unique(am)
am <- am %>% transform(aa_no = as.numeric(am$aa_no))
am_protein_length <- am %>%
  group_by(ID) %>%
  filter(aa_no == max(aa_no)) %>%
  select(`ID`, `aa_no`) %>%
  rename("amino_acids" = "aa_no")
write_tsv(am_protein_length, "am_protein_length.tsv")
am <- arrange_all(am)
am_protein_sequence <- am %>%
  group_by(ID) %>%
  mutate(ref_aa = paste0(ref_aa, collapse = "")) %>%
  rename("aa_sequence" = "ref_aa") %>%
  select(`ID`, `aa_sequence`)
am_protein_sequence <- unique(am_protein_sequence)
am_protein_details <- am_protein_sequence %>% inner_join(am_protein_length, by = "ID")
am_protein_details <- unique(am_protein_details)
am_protein_details$seq_length <- str_length(am_protein_details$aa_sequence)
am_protein_details <- am_protein_details %>%
  mutate(comparison = ifelse(amino_acids == seq_length, TRUE, FALSE))
write_tsv(am_protein_details, "am_protein_details.tsv")
# am_protein_details; 19,302 proteins
# VM with 128GB RAM
```

To get the number of codons in the Ensembl reference sequences and compare these with the protein data I used the following. Note that I checked for stop codons and removed these, as they are not translated into amino acids.

```R
library(tidyverse)
coding_length <- reference_transcripts_full
coding_length$stop_codon <- coding_length$sequence %>% str_detect("TGA$|TAG$|TAA$")
coding_length <- coding_length %>%
  mutate(sequence = ifelse(stop_codon == TRUE, str_sub(sequence, 1, end = -4), sequence))
coding_length$codons <- str_length(coding_length$sequence) / 3
am_full_protein_details <- read_tsv("am_full_protein_details.tsv")
comparison <- coding_length %>% inner_join(am_full_protein_details, by = "ID")
comparison <- comparison %>%
  mutate(comparison = ifelse(codons == amino_acids, TRUE, FALSE))
mismatch <- comparison %>% filter(comparison == FALSE)
# In the full set, 18,223 transcripts matched to 19,409 proteins. 1,186 did not.
am_protein_details <- read_tsv("am_protein_details.tsv")
comparison2 <- coding_length %>% inner_join(am_protein_details, by = "ID")
comparison2 <- comparison2 %>%
  mutate(comparison = ifelse(codons == amino_acids, TRUE, FALSE))
mismatch2 <- comparison2 %>% filter(comparison == FALSE)
# In the partial set, 18,115 transcripts matched to 19,302 proteins. 1,147 did not.
```

1,186 of the transcript lengths in the full set did not match to corresponding protein lengths. Because I had previously found missing transcripts in the 'isoforms' sets, I searched them for the *mismatched* transcripts. For the full set:

```R
library(tidyverse)
am_isoforms_full_list <- read_tsv("am_isoforms_full_list.tsv")
am_isoforms_full_list<- am_isoforms_full_list %>% rename("ID" = "transcript_id")
am_isoforms_full_list$short_ID <- am_isoforms_full_list$ID %>% str_extract("ENST\\d+")
mismatch$short_ID <- mismatch$ID %>% str_extract("ENST\\d+")
mismatch_isoforms <- mismatch %>% inner_join(am_isoforms_full_list, by = "short_ID")
mismatch_isoforms <- mismatch_isoforms %>% rename("ID" = "ID.x", "am_ID" = "ID.y")
write_tsv(mismatch_isoforms, "mismatch_isoforms.tsv")
# 1,036 of the 1,186 mismatched transcript IDs match to IDs in the full 'isoforms' set.
am_isoforms_full <- read_tsv("AlphaMissense_isoforms_aa_substitutions.tsv.gz", comment = "#")
mismatch_isoforms <- read_tsv("mismatch_isoforms.tsv")
mismatch_isoforms <- mismatch_isoforms %>% select(`ID`, `am_ID`)
am_isoforms_full <- am_isoforms_full %>% rename("am_ID" = "transcript_id", "change" = "protein_variant", "score" = "am_pathogenicity")
am_isoforms_full <- mismatch_isoforms %>% left_join(am_isoforms_full, by = "am_ID")
am_isoforms_full <- am_isoforms_full %>% select(`ID`, `change`, `score`)
am_isoforms_full <- unique(am_isoforms_full)
am_isoforms_full <- am_isoforms_full %>% drop_na(score)
am_isoforms_full_rematched <- am_isoforms_full
write_tsv(am_isoforms_full_rematched, "am_isoforms_full_rematched.tsv")
rm(am_isoforms_full_rematched)
am_isoforms_full$aa_no <- am_isoforms_full$change %>% str_extract("\\d+")
am_isoforms_full$ref_aa <- am_isoforms_full$change %>% str_extract("[^0-9]+")
am_isoforms_full <- am_isoforms_full %>% select(`ID`, `aa_no`, `ref_aa`)
am_isoforms_full <- unique(am_isoforms_full)
am_isoforms_full <- am_isoforms_full %>% transform(aa_no = as.numeric(am_isoforms_full$aa_no))
mismatch_protein_length <- am_isoforms_full %>%
  group_by(ID) %>%
  filter(aa_no == max(aa_no)) %>%
  select(`ID`, `aa_no`) %>%
  rename("amino_acids" = "aa_no")
am_isoforms_full <- arrange_all(am_isoforms_full)
mismatch_protein_sequence <- am_isoforms_full %>%
  group_by(ID) %>%
  mutate(ref_aa = paste0(ref_aa, collapse = "")) %>%
  rename("aa_sequence" = "ref_aa") %>%
  select(`ID`, `aa_sequence`)
mismatch_protein_sequence <- unique(mismatch_protein_sequence)
mismatch_protein_details <- mismatch_protein_sequence %>% inner_join(mismatch_protein_length, by = "ID")
mismatch_protein_details <- unique(mismatch_protein_details)
# 1,015 protein sequences generated from the 1,036 matched IDs
reference_transcripts_full <- read_tsv("reference_transcripts_full.tsv")
coding_length <- reference_transcripts_full
coding_length$stop_codon <- coding_length$sequence %>% str_detect("TGA$|TAG$|TAA$")
coding_length <- coding_length %>%
  mutate(sequence = ifelse(stop_codon == TRUE, str_sub(sequence, 1, end = -4), sequence))
coding_length$codons <- str_length(coding_length$sequence) / 3
comparison <- coding_length %>% inner_join(mismatch_protein_details, by = "ID")
comparison <- comparison %>%
  mutate(comparison = ifelse(codons == amino_acids, TRUE, FALSE))
mismatch <- comparison %>% filter(comparison == FALSE)
am_isoforms_full_rematched_to_add <- comparison %>% filter(comparison == TRUE)
write_tsv(am_isoforms_full_rematched_to_add, "am_isoforms_full_rematched_to_add.tsv")
# clear workspace
am_isoforms_full_rematched <- read_tsv("am_isoforms_full_rematched.tsv")
am_isoforms_full_rematched_to_add <- read_tsv("am_isoforms_full_rematched_to_add")
# 962 of the 1,015 protein sequences are of the appropriate lengths
# VM with 128GB RAM
```

For the partial set:

```R
library(tidyverse)
am_protein_details <- read_tsv("am_protein_details.tsv")
comparison2 <- coding_length %>% inner_join(am_protein_details, by = "ID")
comparison2 <- comparison2 %>%
  mutate(comparison = ifelse(codons == amino_acids, TRUE, FALSE))
mismatch2 <- comparison2 %>% filter(comparison == FALSE)
am_isoforms_list <- read_tsv("am_isoforms_list.tsv")
am_isoforms_ist<- am_isoforms_list %>% rename("ID" = "transcript_id")
am_isoforms_list$short_ID <- am_isoforms_list$ID %>% str_extract("ENST\\d+")
mismatch2$short_ID <- mismatch2$ID %>% str_extract("ENST\\d+")
mismatch_isoforms2 <- mismatch2 %>% inner_join(am_isoforms_list, by = "short_ID")
write_tsv(mismatch_isoforms2, "mismatch_isoforms2.tsv")
# 998 of the 1,147 mismatched transcript IDs match to IDs in the partial 'isoforms' set.
am_isoforms <- read_tsv("AlphaMissense_isoforms_hg38.tsv.gz", comment = "#", col_names = FALSE, col_select = c(6, 7, 8))
am_isoforms <- am_isoforms %>% rename("am_ID" = "X6", "change" = "X7", "score" = "X8")
mismatch_isoforms2 <- read_tsv("mismatch_isoforms2.tsv")
mismatch_isoforms2 <- mismatch_isoforms2 %>% rename("ID" = "ID.x", "am_ID" = "ID.y")
mismatch_isoforms2 <- mismatch_isoforms2 %>% select(`ID`, `am_ID`)
am_isoforms <- am_isoforms %>% rename("am_ID" = "transcript_id", "change" = "protein_variant", "score" = "am_pathogenicity")
am_isoforms <- mismatch_isoforms2 %>% left_join(am_isoforms, by = "am_ID")
am_isoforms <- am_isoforms %>% select(`ID`, `change`, `score`)
am_isoforms <- unique(am_isoforms)
am_isoforms <- am_isoforms %>% drop_na(score)
am_isoforms_rematched <- am_isoforms
write_tsv(am_isoforms_rematched, "am_isoforms_rematched.tsv")
rm(am_isoforms_rematched)
am_isoforms$aa_no <- am_isoforms$change %>% str_extract("\\d+")
am_isoforms$ref_aa <- am_isoforms$change %>% str_extract("[^0-9]+")
am_isoforms <- am_isoforms %>% select(`ID`, `aa_no`, `ref_aa`)
am_isoforms <- unique(am_isoforms)
am_isoforms <- am_isoforms %>% transform(aa_no = as.numeric(am_isoforms$aa_no))
mismatch_protein_length <- am_isoforms %>%
  group_by(ID) %>%
  filter(aa_no == max(aa_no)) %>%
  select(`ID`, `aa_no`) %>%
  rename("amino_acids" = "aa_no")
am_isoforms <- arrange_all(am_isoforms)
mismatch_protein_sequence <- am_isoforms %>%
  group_by(ID) %>%
  mutate(ref_aa = paste0(ref_aa, collapse = "")) %>%
  rename("aa_sequence" = "ref_aa") %>%
  select(`ID`, `aa_sequence`)
mismatch_protein_sequence <- unique(mismatch_protein_sequence)
mismatch_protein_details <- mismatch_protein_sequence %>% inner_join(mismatch_protein_length, by = "ID")
mismatch_protein_details <- unique(mismatch_protein_details)
# 980 protein sequences generated from the 998 matched IDs
reference_transcripts_full <- read_tsv("reference_transcripts_full.tsv")
coding_length <- reference_transcripts_full
coding_length$stop_codon <- coding_length$sequence %>% str_detect("TGA$|TAG$|TAA$")
coding_length <- coding_length %>%
  mutate(sequence = ifelse(stop_codon == TRUE, str_sub(sequence, 1, end = -4), sequence))
coding_length$codons <- str_length(coding_length$sequence) / 3
comparison <- coding_length %>% inner_join(mismatch_protein_details, by = "ID")
comparison <- comparison %>%
  mutate(comparison = ifelse(codons == amino_acids, TRUE, FALSE))
mismatch <- comparison %>% filter(comparison == FALSE)
am_isoforms_rematched_to_add <- comparison %>% filter(comparison == TRUE)
write_tsv(am_isoforms_rematched_to_add, "am_isoforms_rematched_to_add.tsv")
# clear workspace
am_isoforms_rematched <- read_tsv("am_isoforms_rematched.tsv")
am_isoforms_rematched_to_add <- read_tsv("am_isoforms_rematched_to_add.tsv")
am_isoforms_rematched_to_add <- am_isoforms_rematched_to_add %>% select(ID)
am_isoforms_rematched_to_add <- am_isoforms_rematched_to_add %>% inner_join(am_isoforms_rematched, by = "ID")
write_tsv(am_isoforms_rematched_to_add, "am_isoforms_rematched_to_add.tsv")
# 844 of the 980 protein sequences are of the appropriate lengths
# VM with 128GB RAM
```

A substantial fraction of the proteins with the wrong numbers of amino acids in the initial versions of the combined AlphaMissense data were found with the "correct" numbers of amino acids in the 'isoforms' sets. To remove the mismatched protein data and replace it with the data that matched from the 'isoforms' sets, I used:

```R
library(tidyverse)
mismatch <- read_tsv("mismatch.tsv")
am_full <- read_tsv("am_full.tsv.gz")
am_isoforms_full_rematched_to_add <- read_tsv("am_isoforms_full_rematched_to_add.tsv")
am_full <- am_full %>% filter(!ID %in% mismatch$ID)
am_full <- rows_insert(am_full, am_isoforms_full_rematched_to_add)
write_tsv(am_full, "am_full.tsv.gz")
# 210,270,379 variants
mismatch2 <- read_tsv("mismatch2.tsv")
am <- read_tsv("am.tsv.gz")
am_isoforms_rematched_to_add <- read_tsv("am_isoforms_rematched_to_add.tsv")
am <- am %>% filter(!ID %in% mismatch2$ID)
am <- rows_insert(am, am_isoforms_rematched_to_add)
write_tsv(am, "am.tsv.gz")
# 64,274,881 variants
# VM with 128GB RAM
```

I rechecked the updated data to make sure it was okay:

```R
library(tidyverse)
am_full <- read_tsv("am_full.tsv.gz")
am_full$aa_no <- am_full$change %>% str_extract("\\d+")
am_full$ref_aa <- am_full$change %>% str_extract("[^0-9]+")
am_full <- am_full %>% select(`ID`, `aa_no`, `ref_aa`)
am_full <- unique(am_full)
am_full <- am_full %>% transform(aa_no = as.numeric(am_full$aa_no))
am_full_protein_length <- am_full %>%
  group_by(ID) %>%
  filter(aa_no == max(aa_no)) %>%
  select(`ID`, `aa_no`) %>%
  rename("amino_acids" = "aa_no")
am_full <- arrange_all(am_full)
am_full_protein_sequence <- am_full %>%
  group_by(ID) %>%
  mutate(ref_aa = paste0(ref_aa, collapse = "")) %>%
  rename("aa_sequence" = "ref_aa") %>%
  select(`ID`, `aa_sequence`)
am_full_protein_sequence <- unique(am_full_protein_sequence)
am_full_protein_details <- am_full_protein_sequence %>% inner_join(am_full_protein_length, by = "ID")
am_full_protein_details <- unique(am_full_protein_details)
am_full_protein_details$seq_length <- str_length(am_full_protein_details$aa_sequence)
am_full_protein_details <- am_full_protein_details %>%
  mutate(comparison = ifelse(amino_acids == seq_length, TRUE, FALSE))
reference_transcripts_full <- read_tsv("reference_transcripts_full.tsv")
coding_length <- reference_transcripts_full
coding_length$stop_codon <- coding_length$sequence %>% str_detect("TGA$|TAG$|TAA$")
coding_length <- coding_length %>%
  mutate(sequence = ifelse(stop_codon == TRUE, str_sub(sequence, 1, end = -4), sequence))
coding_length$codons <- str_length(coding_length$sequence) / 3
comparison <- coding_length %>% inner_join(am_full_protein_details, by = "ID")
comparison <- comparison %>%
  mutate(comparison = ifelse(codons == amino_acids, TRUE, FALSE))
mismatch <- comparison %>% filter(comparison == FALSE)
# full data; 19,185 transcripts matching proteins of appropriate lengths
am <- read_tsv("am.tsv.gz")
am$aa_no <- am$change %>% str_extract("\\d+")
am$ref_aa <- am$change %>% str_extract("[^0-9]+")
am <- am %>% select(`ID`, `aa_no`, `ref_aa`)
am <- unique(am)
am <- am %>% transform(aa_no = as.numeric(am$aa_no))
am_protein_length <- am %>%
  group_by(ID) %>%
  filter(aa_no == max(aa_no)) %>%
  select(`ID`, `aa_no`) %>%
  rename("amino_acids" = "aa_no")
am <- arrange_all(am)
am_protein_sequence <- am %>%
  group_by(ID) %>%
  mutate(ref_aa = paste0(ref_aa, collapse = "")) %>%
  rename("aa_sequence" = "ref_aa") %>%
  select(`ID`, `aa_sequence`)
am_protein_sequence <- unique(am_protein_sequence)
am_protein_details <- am_protein_sequence %>% inner_join(am_protein_length, by = "ID")
am_protein_details <- unique(am_protein_details)
am_protein_details$seq_length <- str_length(am_protein_details$aa_sequence)
am_protein_details <- am_protein_details %>%
  mutate(comparison = ifelse(amino_acids == seq_length, TRUE, FALSE))
reference_transcripts_full <- read_tsv("reference_transcripts_full.tsv")
coding_length <- reference_transcripts_full
coding_length$stop_codon <- coding_length$sequence %>% str_detect("TGA$|TAG$|TAA$")
coding_length <- coding_length %>%
  mutate(sequence = ifelse(stop_codon == TRUE, str_sub(sequence, 1, end = -4), sequence))
coding_length$codons <- str_length(coding_length$sequence) / 3
comparison <- coding_length %>% inner_join(am_protein_details, by = "ID")
comparison <- comparison %>%
  mutate(comparison = ifelse(codons == amino_acids, TRUE, FALSE))
mismatch <- comparison %>% filter(comparison == FALSE)
# Partial data; 18,999 transcripts matching proteins of appropriate lengths
# VM with 128GB RAM
```

At this point, 19,689  Ensembl canonical transcripts had been matched to 19,185 and 18,999 proteins in the AlphaMissense complete and partial datasets, respectively. The 'canonical' AlphaMissense proteins were matched by Ensembl IDs mapped from Uniprot IDs. The 'isoforms' AlphaMissense proteins were matched by Ensembl IDs from which version numbers were removed. The only IDs that had not been cross-checked to the 19,689  Ensembl canonical transcripts were the Ensembl IDs given in the "AlphaMissense_hg38.tsv.gz" set. Therefore, I cross-checked these against the few missing ones:

```R
library(tidyverse)
df <- read_tsv("df.tsv")
df$short_ID <- df$ID %>% str_extract("ENST\\d+")
am_full <- read_tsv("am_full.tsv.gz")
am_full <- am_full %>% select(ID)
am_full <- unique(am_full)
am_full$short_ID <- am_full$ID %>% str_extract("ENST\\d+")
hg38_canonical <- read_tsv("AlphaMissense_hg38.tsv.gz", comment = "#", col_names = FALSE, col_select = c(6,7,8,9))
hg38_canonical <- hg38_canonical %>% rename("uniprot_ID" = "X6", "am_ID" = "X7", "change" = "X8", "score" = "X9")
hg38_ensembl_IDs <- hg38_canonical %>% select(am_ID)
hg38_ensembl_IDs <- unique(hg38_ensembl_IDs)
hg38_ensembl_IDs$short_ID <- hg38_ensembl_IDs$am_ID %>% str_extract("ENST\\d+")
hg38_ensembl_IDs_not_in_full <- hg38_ensembl_IDs %>% anti_join(am_full, by = "short_ID")
matches1 <- hg38_ensembl_IDs_not_in_full %>% inner_join(df, by = "short_ID")
matches1 <- matches1 %>% filter(!ID == "ENST00000569103.2")
matches1 <- matches1 %>% select(`ID`, `am_ID`, `sequence`)
coding_length <- matches1
coding_length$stop_codon <- coding_length$sequence %>% str_detect("TGA$|TAG$|TAA$")
coding_length <- coding_length %>%
  mutate(sequence = ifelse(stop_codon == TRUE, str_sub(sequence, 1, end = -4), sequence))
coding_length$codons <- str_length(coding_length$sequence) / 3
am <- hg38_canonical
am$aa_no <- am$change %>% str_extract("\\d+")
am$ref_aa <- am$change %>% str_extract("[^0-9]+")
am <- am %>% select(`am_ID`, `aa_no`, `ref_aa`)
am <- unique(am)
am <- am %>% transform(aa_no = as.numeric(am$aa_no))
am_protein_length <- am %>%
  group_by(am_ID) %>%
  filter(aa_no == max(aa_no)) %>%
  select(`am_ID`, `aa_no`) %>%
  rename("amino_acids" = "aa_no")
am <- arrange_all(am)
am_protein_sequence <- am %>%
  group_by(am_ID) %>%
  mutate(ref_aa = paste0(ref_aa, collapse = "")) %>%
  rename("aa_sequence" = "ref_aa") %>%
  select(`am_ID`, `aa_sequence`)
am_protein_sequence <- unique(am_protein_sequence)
am_protein_details <- am_protein_sequence %>% inner_join(am_protein_length, by = "am_ID")
am_protein_details <- unique(am_protein_details)
am_protein_details$seq_length <- str_length(am_protein_details$aa_sequence)
comparison <- coding_length %>% inner_join(am_protein_details, by = "am_ID")
comparison <- comparison %>%
  mutate(comparison = ifelse(codons == amino_acids, TRUE, FALSE))
mismatch <- comparison %>% filter(comparison == FALSE)
hg38_rematched_to_add <- comparison %>% filter(comparison == TRUE)
write_tsv(hg38_rematched_to_add, "hg38_rematched_to_add.tsv")
# 24 additional transcripts and proteins to add to 'full' set
# VM with 128GB RAM
```

To add the 24 additional transcripts to the full set, I used:

```R
library(tidyverse)
am_full <- read_tsv("am_full.tsv.gz")
hg38_rematched_to_add <- read_tsv("hg38_rematched_to_add.tsv")
hg38_canonical <- read_tsv("AlphaMissense_hg38.tsv.gz", comment = "#", col_names = FALSE, col_select = c(6,7,8,9))
hg38_canonical <- hg38_canonical %>% rename("uniprot_ID" = "X6", "am_ID" = "X7", "change" = "X8", "score" = "X9")
hg38_rematched_to_add <- hg38_rematched_to_add %>% inner_join(hg38_canonical, by = "am_ID")
hg38_rematched_to_add <- hg38_rematched_to_add %>% select(`ID`, `change`, `score`)
am_full <- rows_insert(am_full, hg38_rematched_to_add)
# 210,329,860 variants; 19,209 transcripts
# VM with 128GB RAM
```

I then compared the full and partial sets to see which transcripts weren't present in both:

```R
am_full <- read_tsv("am_full.tsv.gz")
am_full <- am_full %>% select(ID)
am_full <- unique(am_full)
am <- read_tsv("am.tsv.gz")
am <- am %>% select(ID)
am <- unique(am)
am_full_not_am <- am_full %>% anti_join(am, by = "ID")
write_tsv(am_full_not_am, "am_full_not_am.tsv")
am_not_am_full <- am %>% anti_join(am_full, by = "ID")
write_tsv(am_not_am_full, "am_not_am_full.tsv")
# 1 transcript in am not am_full; 211 transcripts in am_full not am
# VM with 128GB RAM
```

Although the 'partial' set contains 'scores for every amino acid change that it is possible to produce with a single nucleotide change' and the 'full' set contains 'scores for every possible amino acid change at every position', it's acceptable for the completeness that the 'partial' set contains a few proteins with full sets of scores and the 'full' set contains a few proteins with partial sets of scores.

```R
am_full <- read_tsv("am_full.tsv.gz")
am <- read_tsv("am.tsv.gz")
am_full_not_am <- am_full_not_am %>% inner_join(am_full, by = "ID")
write_tsv(am_full_not_am, "am_full_not_am.tsv")
am_not_am_full <- am_not_am_full %>% inner_join(am, by = "ID")
write_tsv("am_not_am_full.tsv")
am <- rows_insert(am, am_full_not_am)
am_full <- rows_insert(am_full, am_not_am_full)
write_tsv(am, "am1.tsv.gz")
write_tsv(am_full, "am_full1.tsv.gz")
# am1.tsv.gz (376.4MB) contains 19,210 transcripts and 66,426,471 variants
# am_full1.tsv.gz (1.1GB) contains 19,210 transcripts and 210,335,478 variants
# VM with 128GB RAM
```

Finally, I added the additional transcript nucleotide sequences to the reference set:

```R
am_full <- read_tsv("am_full1.tsv.gz")
am_full <- am_full %>% select(`ID`)
am_full <- unique(am_full)
reference <- read_tsv("reference_transcripts_full.tsv")
reference <- am_full %>% inner_join(reference, by = "ID")
df <- read_tsv("df.tsv")
extra <- am_full %>% filter(!`ID` %in% reference$ID)
extra <- extra %>% inner_join(df, by = "ID")
reference <- rows_insert(reference, extra)
write_tsv(reference, "reference1.tsv")
# 19,210 reference sequences; 33,257,151 nucleotides
# VM with 128GB RAM
```

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