use csv::ReaderBuilder;
use std::fs::OpenOptions;
use std::io::prelude::*;
use crate::r_n_g;

// download reference transcriptome from https://storage.googleapis.com/am_data_transfer_test1/reference1.tsv
// download full AlphaMissense dataset from https://storage.googleapis.com/am_data_transfer_test1/am_full1.tsv.gz
// download partial AlphaMissense dataset from https://storage.googleapis.com/am_data_transfer_test1/am1.tsv.gz

struct RowData {
    run_number: u32,
    id: String,
    sequence: String,
}

struct Substitutions {
    run_number: u32,
    id: String,
    sequence: String,
    substituted_positions: Vec<u32>
}

struct SubstitutedSequences {
    run_number: u32,
    id: String,
    original_sequence: String,
    substituted_sequence: String,
    substituted_positions: Vec<u32>
}

struct Codons {
    run_number: u32,
    id: String,
    original_codon: String,
    codon_number: u32,
    substituted_codon: String
}


pub fn parse_run(num_of_times: u32) {

    // call parse_file a specified number of times
    // count the run number
    // pass parse_file the run number
    
    let mut run_number: u32 = 0;
    for _ in 0..num_of_times {
        run_number += 1;
        println!("Running parse_file for run number: {}", run_number);
        parse_file(run_number);
    }
    
}


fn parse_file(run_number: u32) {

    // parse the reference transcript file line-by-line and pass to the parse_sequence function

    let file_path = "reference_toy.tsv";
    let mut reader = ReaderBuilder::new().delimiter(b'\t').has_headers(true).from_path(file_path).unwrap();

    for result in reader.records() {
        let record = result.unwrap();
        if record.len() < 2 {
            continue;
        }
        let id = record[0].to_string();
        let sequence = record[1].to_string();
        let row_data = RowData {
            run_number,
            id,
            sequence,
        };
        parse_sequence(row_data); // This is a blocking call. It will wait until parse_sequence has completed.
    }
}


fn parse_sequence (row_data: RowData) {

    // iterate character-by-character through the sequence, calling r_n_g to determine if a substitution is made at each position
    // store run number, ID, corresponding sequence, and substituted positions in a data structure called "Substitutions"
    // pass the data to the substitute_sequence function if a substitution is made at any position.
    
    let mut substitutions = Substitutions {
        run_number: row_data.run_number,
        id: row_data.id.clone(),
        sequence: row_data.sequence.clone(),
        substituted_positions: Vec::new(),
    };

    // substitution thresholds can be changed here
    // heart (Itpa-null IMP:AMP - control IMP:AMP) / 4 = 0.002432
    // kidney (Itpa-null IMP:AMP - control IMP:AMP) / 4 = 0.000635
    // brain (Itpa-null IMP:AMP - control IMP:AMP) / 4 = 0.00033493
    //
    // proportionate "substitution rates" with T7 are @A=0.109295, @T=0.099429, @G=0.548914, @C=0.242362

    for (i, character) in row_data.sequence.chars().enumerate() {
        match character {
            'A' => {
                if r_n_g() < 0.00106322 {
                    substitutions.substituted_positions.push(i as u32);
                }
            },
            'T' => {
                if r_n_g() < 0.00096724 {
                    substitutions.substituted_positions.push(i as u32);
                }
            },
            'G' => {
                if r_n_g() < 0.00533984 {
                    substitutions.substituted_positions.push(i as u32);
                }
            },
            'C' => {
                if r_n_g() < 0.0023577 {
                    substitutions.substituted_positions.push(i as u32);
                }
            },
            _ => println!("Warning: Unexpected character {} at position {}", character, i),
        }
    }

    // not passing empty data

    if !substitutions.substituted_positions.is_empty() {
        substitute_sequence(substitutions);
    }
    
}


fn substitute_sequence (substitutions: Substitutions) {

    // copy the sequence from the substitutions data structure
    // substitute the character at each specified position with an 'I'
    // store the original and substituted sequences in the data structure called "SubstitutedSequences"
    // pass the SubstitutedSequences data to the extract_codons function

    let mut substituted_sequence = substitutions.sequence.clone();

    let substituted_positions = substitutions.substituted_positions.clone(); // Clone the substituted_positions vector

    for position in substitutions.substituted_positions {
        let index = position as usize;
        substituted_sequence.replace_range(index..index+1, "I");
    }

    let substituted_sequences = SubstitutedSequences {
        run_number: substitutions.run_number,
        id: substitutions.id.clone(),
        original_sequence: substitutions.sequence,
        substituted_sequence,
        substituted_positions,
    };

    extract_codons(substituted_sequences);
}


fn extract_codons (substituted_sequences: SubstitutedSequences) {

    // calculates codon number of substitutions
    // extracts the original and substituted codons from the original and substituted sequences
    // note that in Rust, indexes are zero-based
    // the dedup is used because a twice-substituted codon will still only be translated once
    // probably only codon number needs to be compared for dedup to work 

    let mut codons_vec: Vec<Codons> = Vec::new();

    for &position in &substituted_sequences.substituted_positions {
        let codon_number = ((position as f64 + 1f64) / 3f64).ceil() as u32;

        let original_codon_start = if ((position as f64 + 1f64) / 3f64) % 3f64 == 0.0 { position - 2 } else { position - position % 3 };
        let original_codon_end = original_codon_start + 3;
        let original_codon = substituted_sequences.original_sequence[original_codon_start as usize..original_codon_end as usize].to_string();

        let substituted_codon_start = if ((position as f64 + 1f64) / 3f64) % 3f64 == 0.0 { position - 2 } else { position - position % 3 };
        let substituted_codon_end = substituted_codon_start + 3;
        let substituted_codon = substituted_sequences.substituted_sequence[substituted_codon_start as usize..substituted_codon_end as usize].to_string();

        let codon = Codons {
            run_number: substituted_sequences.run_number,
            id: substituted_sequences.id.clone(),
            original_codon,
            codon_number,
            substituted_codon,
        };

        codons_vec.push(codon);
    }

    impl PartialEq for Codons {
        fn eq(&self, other: &Self) -> bool {
            // Implement the equality comparison logic for Codons struct
            // Return true if the two Codons are equal, false otherwise
            // You need to compare all the fields of Codons struct for equality
            self.run_number == other.run_number &&
            self.id == other.id &&
            self.original_codon == other.original_codon &&
            self.codon_number == other.codon_number &&
            self.substituted_codon == other.substituted_codon
        }
    }

    codons_vec.dedup();

    for codon in codons_vec {
        write_codon_output(codon);
    }

}


fn write_codon_output (codon: Codons) {

    // write the "codons" data to a file called "codon_output.tsv"
    // the file is tab separated with columns for run number, ID, original codon, codon number and substituted codon.

    let mut file = OpenOptions::new()
        .write(true)
        .create(true) // Create the file if it doesn't exist
        .append(true)
        .open("codon_output_toy.tsv")
        .unwrap();

    if let Err(e) = writeln!(file, "{}\t{}\t{}\t{}\t{}", codon.run_number, codon.id, codon.original_codon, codon.codon_number, codon.substituted_codon) {
        eprintln!("Couldn't write to file: {}", e);
    }
}