### approach_1_inefficient_RNG

This seems like the most straightforward intuitive approach. It's the one I'm going to try. I should be able to get it to work. It's probably the least computationally efficient, but with enough computing power that might not matter so much.

**data generation**

1. Approximate the probability of inosine misincorporation in place of each other nucleotide with the data from the literature above.

2. Read the sequences of transcripts and store as an appropriate type. Each transcript is a vector with a sequence of nucleotides, each with a defined position.

3. Iterate through each nucleotide of each transcript. Use a random number generator (RNG) to determine whether each nucleotide is substituted or not. For example, if there's a 1% chance of substitution and a RNG generating a number between 0 and 1, then a roll of between 0 and 0.01 means a substitution, and a roll of between 0.01 and 1 means no substitution.

4. In case of a substitution(s) in a codon, capture the codons-sequence (original and substituted), its position in the transcript, and the transcript identifier. Save these data for each run-through of the transcriptome.

**substitution effect**

1. To transform the original and substituted codons into original and substituted amino acids I will initially assume that 'I' is interpreted as 'G'. This has been an assumption in the literature in the past. It means that determining original and substituted amino acids just requires [the genetic code](https://en.wikipedia.org/wiki/Genetic_code#/media/File:GeneticCode21-version-2.svg) as a lookup table.

2. The changes that replace an amino acid with a different amino acid (missense changes) will be cross-referenced with AlphaMissense. The changes that replace an amino acid with the same amino acid (synonymous changes) will be noted but not cross-referenced with AlphaMissense. The changes that introduce premature stop codons (nonsense changes) will also be noted but not cross-referenced with AlphaMissense, though they may be used in later analysis. Note that nonsense changes will often completely nullify a transcript because of the phenomenon of [nonsense-mediated decay](https://en.wikipedia.org/wiki/Nonsense-mediated_decay), but that in terminal exons they may instead cause production of truncated protein.

* In reality, 'I' may not always be interpreted as 'G'. [This](https://doi.org/10.1093%2Fnar%2Fgky1163) paper shows that. If I can get everything to work where 'I' is always interpreted as 'G', then I'll refine the code to include the possibility of alternative substitutions. The possibility of 'stalling changes' will also be considered. These possibilities will require another layer of iteration with a RNG.

**substitution functional effect**

1. With the defined amino acid changes per-run, the corresponding scores in AlphaMissense can be looked up (since the dataset contains every possible missense change).

2. A first-pass metric for the 'per-transcript' effect will just be adding up the AlphaMissense scores for each transcript.

3. After enough runs, the mean and variance for the predicted functional effect on each transcript should stabilise.