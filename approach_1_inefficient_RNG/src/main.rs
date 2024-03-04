use csv::ReaderBuilder;
use rand::Rng;
use std::fs::OpenOptions;
use std::io::prelude::*;

// download reference transcriptome from https://storage.googleapis.com/am_data_transfer_test1/reference1.tsv
// download full AlphaMissense dataset from https://storage.googleapis.com/am_data_transfer_test1/am1.tsv.gz
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

fn main() {
    parse_run(1);
}

fn parse_run(num_of_times: u32) {

    // call parse_file a specified number of times
    // count the run number
    // pass parse_file the run number
    
    let mut run_number: u32 = 0;
    for _ in 0..num_of_times {
        run_number += 1;
        parse_file(run_number);
    }

}

fn parse_file(run_number: u32) {

    // parse the reference transcript file line-by-line and pass to the parse_sequence function

    let file_path = "reference1.tsv";
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
    
    let mut substitutions = Substitutions {
        run_number: row_data.run_number,
        id: row_data.id.clone(),
        sequence: row_data.sequence.clone(),
        substituted_positions: Vec::new(),
    };

    for (i, character) in row_data.sequence.chars().enumerate() {
        match character {
            'A' => {
                if r_n_g() < 0.0001 {
                    substitutions.substituted_positions.push(i as u32);
                }
            },
            'T' => {
                if r_n_g() < 0.0002 {
                    substitutions.substituted_positions.push(i as u32);
                }
            },
            'G' => {
                if r_n_g() < 0.0003 {
                    substitutions.substituted_positions.push(i as u32);
                }
            },
            'C' => {
                if r_n_g() < 0.0004 {
                    substitutions.substituted_positions.push(i as u32);
                }
            },
            _ => println!("Warning: Unexpected character {} at position {}", character, i),
        }
    }

    if !substitutions.substituted_positions.is_empty() {
        substitute_sequence(substitutions);
    }
    // iterate character-by-character through the sequence, calling the RNG to determine if a substitution is made at each position
    // store run number, ID, corresponding sequence, and substituted positions in a data structure called "substitutions"
    // call parseFile when the end of the sequence is reached
    // pass the data structure to the substituteSequence function if a substitution is made at any position.

}

fn r_n_g () -> f64 {
    
    // generate a random number between 0 and 1
    // return the number
    // this is a separate function to allow for the RNG to be swapped-out if desired
    
    let mut rng = rand::thread_rng();
    return rng.gen::<f64>();
}

fn substitute_sequence (substitutions: Substitutions) {

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
    
    
    // copy the sequence from the substitutions data structure
    // substitute the character at each specified position with an 'I'
    // store the original and substituted sequences in the data structure called "SubstitutedSequences"
    // pass the SubstitutedSequences data to the extract_codons function

}

fn extract_codons (substituted_sequences: SubstitutedSequences) {

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

    // iterate through the original sequence, extracting the codon at each substituted position
    // the codon is 3 characters long and starts two characters before a substituted position perfectly divisible by 3
    // without a remainder. It starts at one character before a substituted position with a remainder of 2, and at a 
    // the same position as a substituted position with a remainder of 1.
    // the 'codon number' is the floored integer division of the substituted position by 3
    // iterate through the substituted sequence, extracting the codon at each substituted position
    // store the original and substituted codons in a data structure called "codons"
    // the "codons" data structure should contain the run numbers, IDs, original codons, codon numbers and substituted codons.
    // 
    // note that where two characters in the same codon are substituted, duplicate entries will be made in the "codons" data structure.
    // duplicate entries should be removed before the data structure is passed to the writeCodonOutput function.

}



fn write_codon_output (codon: Codons) {

    // write the "codons" data structure to a file called "codon_output.txt"
    // the file should be tab separated with columns for run number, ID, original codon, codon number and substituted codon.

    let mut file = OpenOptions::new()
        .write(true)
        .create(true) // Create the file if it doesn't exist
        .append(true)
        .open("codon_output.tsv")
        .unwrap();

    if let Err(e) = writeln!(file, "{}\t{}\t{}\t{}\t{}", codon.run_number, codon.id, codon.original_codon, codon.codon_number, codon.substituted_codon) {
        eprintln!("Couldn't write to file: {}", e);
    }
}