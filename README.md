## leveraging AlphaMissense to predict the functional effects of transcription errors in coding RNA

Coding RNA is transcribed from template DNA with high - but not perfect - fidelity. Stochastic (random) transcription errors affecting individual RNA molecules can introduce sequence changes that cause amino acid substitutions within the proteins they encode.

AlphaMissense (1) is a dataset published by scientists at Google DeepMind that scores the functional effects of all possible amino acid substitutions in humans. It is intended primarily to inform the interpretation of potentially pathogenic genomic variants.

Our project will repeatedly simulate probabilistic transcription errors across the 39 million nucleotides of the canonical human coding transcriptome in silico. We will then leverage the 216 million substitutions scored in AlphaMissense as a lookup table to predict the functional effects of transcriptome-wide errors on a per-transcript basis. Our hypothesis is that these data will distinguish error-sensitive transcripts that contribute to disease pathogenesis.

Our particular research focus is severe ITPase deficiency, in which pathological inosine misincorporation into RNA occurs at high frequency. This rare disease has been the subject of three of our recent peer-reviewed publications (2-4). The endpoint of the project will be functioning open-source software, the data it has produced, and a publication describing these. We anticipate publication either as a standalone work or contributary to a larger study for which microarray, RNAseq, long-read sequencing and proteomics' experiments are already complete.

The wider significance of our work includes known or readily-inferred relevance to heart disease, chemotherapy, and ageing. Moreover, the 'central dogma' of molecular biology is concerned with the flow of information from DNA to RNA to protein. The effects of 'noise' on this flow are of general interest.

Together with funding, we hope that ResearchHub will provide a forum through which open development can proceed. We hope to benefit from the critical assessment and suggestions of the community as we progress.
 

1.	Cheng J, Novati G, Pan J, Bycroft C, Žemgulytė A, Applebaum T, et al. Accurate proteome-wide missense variant effect prediction with AlphaMissense. Science. 2023 Sep 22;381(6664):eadg7492. 
2.	Handley MT, Reddy K, Wills J, Rosser E, Kamath A, Halachev M, et al. ITPase deficiency causes a Martsolf-like syndrome with a lethal infantile dilated cardiomyopathy. Plagnol V, editor. PLoS Genet. 2019 Mar 11;15(3):e1007605. 
3.	Schroader JH, Jones LA, Meng R, Shorrock HK, Richardson JI, Shaughnessy SM, et al. Disease-associated inosine misincorporation into RNA hinders translation. Nucleic Acids Res. 2022 Sep 9;50(16):9306–18. 
4.	Schroader JH, Handley MT, Reddy K. Inosine triphosphate pyrophosphatase: A guardian of the cellular nucleotide pool and potential mediator of RNA function. WIREs RNA. 2023 Sep;14(5):e1790. 


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


### approach_2_efficient_RNG

This seems like it would be a more efficient approach but I don't know how to do it. It may be trivial to a statistician or computer scientist, but it's not trivial to me. Once I get through 'approach_1', I may have learned enough to attempt it.

* In approach_2, instead of 'rolling the dice' for each nucleotide of each transcript, the RNG would be used to select the nucleotides at which substitutions happen. That is, the probabilities of inosine misincorporation in place of each other nucleotide would be used to calculate the range of the number of times the dice is rolled. Then, each roll of the dice would return a substituted codon, its position, and the transcript identifier.

Because this approach is more opaque, a formal mathematical proof that it the data it generates is equivalent to that generated by approach_1 would be necessary for it to be reportable.


### approach_3_no_RNG

This seems like it would be the most efficient approach. Again I don't know how to do it. Again it may be trivial to a statistician or computer scientist. Again, it's not trivial to me.

Approach_3 works backwards from AlphaMissense rather than forwards from the transcriptome. Instead of approximate mean and variance for the functional effects on each transcript, it would give absolute ones.

1. Calculate the probability of every possible 'positionless' codon change from the approximated probabilities of each nucleotide substitution.

2. Use these data to calculate the probabilities of every possible 'positionless' amino acid change.

3. Iterate through AlphaMissense and assign a probability to every amino acid change.

4. Use that probability and the AlphaMissense scores to determine exact mean and variance for the functional effects on each transcript.


## A note about Open Science

I'm interested in Open Science. As an experimentalist, I'm going to try it like this:

1. Do a science project in the open.

2. Consider what worked and what didn't.

3. Make adjustments.

4. Repeat.

This is my test project. For me, science is about finding things out. With this project, I'd like to find out whether some transcripts are more sensitive to transcription errors than others. With publicly-available data and a computer I think I can. I'm not sitting back and expecting anyone else to do it. However, I hope to find out what I want to know faster with your ongoing critical assessment and suggestions. If you like, you could even try to do this project more quickly and better than me for the kudos. That way, I'll find out what I want to know *even faster*.


