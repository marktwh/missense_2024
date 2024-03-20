use csv::ReaderBuilder;
use rand::Rng;
use std::fs::OpenOptions;
use std::collections::HashMap;
use std::io::prelude::*;
use crate::determine_translation::determine_translation;
use crate::Codons;


struct OriginalTranslations {
    run_number: u32,
    id: String,
    original_codon: String,
    original_translation: String,
    codon_number: u32,
    substituted_codon: String,
}

struct SubstitutedTranslations {
    run_number: u32,
    id: String,
    original_codon: String,
    original_translation: String,
    codon_number: u32,
    substituted_codon: String,
    substituted_translation: String,
}

#[derive(Clone)]
struct Flags {
    run_number: u32,
    id: String,
    original_codon: String,
    original_translation: String,
    codon_number: u32,
    substituted_codon: String,
    substituted_translation: String,
    flag: String
}


pub fn process_codons_file () {

    let original_codon_lookup_table = load_original_codon_lookup();
    let substituted_codon_lookup_table = load_substituted_codon_lookup();
    parse_codon_output(&original_codon_lookup_table, &substituted_codon_lookup_table);
    
}


fn parse_codon_output (original_codon_lookup_table: &HashMap<String, String>, substituted_codon_lookup_table: &HashMap<String, String>) {

    // parse the codon_output file and print the contents to the console

    let file_path = "codon_output_toy.tsv";
    let mut reader = ReaderBuilder::new().delimiter(b'\t').has_headers(false).from_path(file_path).unwrap();

    for result in reader.records() {
        let record = result.unwrap();
        if record.len() < 4 {
            continue;
        }
        
        let run_number = record[0].parse::<u32>().unwrap();
        let id = record[1].to_string();
        let original_codon = record[2].to_string();
        let codon_number = record[3].parse::<u32>().unwrap();
        let substituted_codon = record[4].to_string();

        let codon = Codons {
            run_number,
            id,
            original_codon,
            codon_number,
            substituted_codon,
        };

        let original_translation: OriginalTranslations = original_codon_lookup(codon, &original_codon_lookup_table);
        let substituted_translation: SubstitutedTranslations = substituted_codon_lookup(original_translation, &substituted_codon_lookup_table);
        if let Some(non_synonymous_translation) = remove_synonymous(substituted_translation) {
            select_missense(non_synonymous_translation);
        }

    }

    // call original_codon_lookup function with change
}


fn original_codon_lookup (codon: Codons, original_codon_lookup: &HashMap<String, String>) -> OriginalTranslations {

    
    let original_translation = match original_codon_lookup.get(&codon.original_codon) {
        Some(translation) => translation.clone(),
        None => {
            eprintln!("Error: Original codon not found in lookup table");
            "".to_string()
        }
    };

    OriginalTranslations {
        run_number: codon.run_number,
        id: codon.id.clone(),
        original_codon: codon.original_codon.clone(),
        original_translation,
        codon_number: codon.codon_number,
        substituted_codon: codon.substituted_codon.clone(),
    }
    
    // return original_translation

}


fn substituted_codon_lookup (original_translation: OriginalTranslations, substituted_codon_lookup: &HashMap<String, String>) -> SubstitutedTranslations {

    let substituted_translation = match substituted_codon_lookup.get(&original_translation.substituted_codon) {
        Some(translation) => translation.clone(),
        None => {
            // Call the determine_substituted_codon function when the translation is not found
            determine_translation(&original_translation.substituted_codon)

        }
    };


    SubstitutedTranslations {
        run_number: original_translation.run_number,
        id: original_translation.id.clone(),
        original_codon: original_translation.original_codon.clone(),
        original_translation: original_translation.original_translation.clone(),
        codon_number: original_translation.codon_number,
        substituted_codon: original_translation.substituted_codon.clone(),
        substituted_translation: substituted_translation,
    }

}


fn remove_synonymous(substituted_translation: SubstitutedTranslations) -> Option<SubstitutedTranslations> {
    // Check if the original and substituted translations are different
    if substituted_translation.original_translation != substituted_translation.substituted_translation {
        // If they are different, return the struct
        Some(substituted_translation)
    } else {
        // If they are the same, return None
        None
    }
}


fn select_missense (non_synonymous_translation: SubstitutedTranslations) {
    
    let mut flag = String::new();

    if non_synonymous_translation.codon_number == 1 {
        flag = "STARTLOSS".to_string();
    } else if non_synonymous_translation.original_translation == "*" {
        flag = "STOPLOSS".to_string();
    } else if non_synonymous_translation.substituted_translation == "*" {
        flag = "NONSENSE".to_string();
    } else if non_synonymous_translation.substituted_translation == "#" {
        flag = "STALL".to_string();
    } else {
        flag = "MISSENSE".to_string();
    }

    let flag = Flags {
        run_number: non_synonymous_translation.run_number,
        id: non_synonymous_translation.id.clone(),
        original_codon: non_synonymous_translation.original_codon.clone(),
        original_translation: non_synonymous_translation.original_translation.clone(),
        codon_number: non_synonymous_translation.codon_number,
        substituted_codon: non_synonymous_translation.substituted_codon.clone(),
        substituted_translation: non_synonymous_translation.substituted_translation.clone(),
        flag,
    };

    if flag.flag == "MISSENSE" {
        write_missense_changes(flag.clone());
    } else {
        write_other_changes(flag.clone());
    }

    write_archive(flag);

}


fn write_missense_changes (flag: Flags) {

    let mut file = OpenOptions::new()
        .write(true)
        .create(true) // Create the file if it doesn't exist
        .append(true)
        .open("missense_changes_toy.tsv")
        .unwrap();

    let change = format!("{}{}{}", flag.original_translation, flag.codon_number, flag.substituted_translation);

    if let Err(e) = writeln!(file, "{}\t{}\t{}",
    flag.run_number,
    flag.id,
    change,
    ) {
        eprintln!("Couldn't write to file: {}", e);
    }

}


fn write_other_changes (flag: Flags) {

    // write the "other_changes" data to a file called "other_changes.tsv"
    // the file is tab separated with columns for run number, ID, original translation, substituted translation, and flag.
    
    let mut file = OpenOptions::new()
        .write(true)
        .create(true) // Create the file if it doesn't exist
        .append(true)
        .open("other_changes_toy.tsv")
        .unwrap();

    if let Err(e) = writeln!(file, "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
    flag.run_number,
    flag.id,
    flag.original_codon,
    flag.original_translation,
    flag.codon_number,
    flag.substituted_codon,
    flag.substituted_translation,
    flag.flag,
    ) {
        eprintln!("Couldn't write to file: {}", e);
    }

    
}


fn write_archive (flag: Flags) {

    // write the "other_changes" data to a file called "other_changes.tsv"
    // the file is tab separated with columns for run number, ID, original translation, substituted translation, and flag.
    
    let mut file = OpenOptions::new()
        .write(true)
        .create(true) // Create the file if it doesn't exist
        .append(true)
        .open("archive_toy.tsv")
        .unwrap();

    if let Err(e) = writeln!(file, "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
    flag.run_number,
    flag.id,
    flag.original_codon,
    flag.original_translation,
    flag.codon_number,
    flag.substituted_codon,
    flag.substituted_translation,
    flag.flag,
    ) {
        eprintln!("Couldn't write to file: {}", e);
    }

}


fn load_original_codon_lookup () -> HashMap<String, String> {

    // load the codon lookup table from a file called "codon_table1.tsv"
    // the file is tab separated with columns for codon and translation.
    // store the data in a HashMap called "original_codon_lookup"

    let mut original_codon_lookup_table: HashMap<String, String> = HashMap::new();

    let file_path = "codon_table1.txt";
    let mut reader = ReaderBuilder::new().delimiter(b'\t').has_headers(false).from_path(file_path).unwrap();

    for result in reader.records() {
        let record = result.unwrap();
        if record.len() < 1 {
            continue;
        }
        let codon = record[0].to_string();
        let translation = record[1].to_string();
        original_codon_lookup_table.insert(codon, translation);
    }

    return original_codon_lookup_table;

}


fn load_substituted_codon_lookup () -> HashMap<String, String> {

    // load the codon lookup table from a file called "codon_table1.tsv"
    // the file is tab separated with columns for codon and translation.
    // store the data in a HashMap called "original_codon_lookup"

    let mut substituted_codon_lookup_table: HashMap<String, String> = HashMap::new();

    // let file_path = "inosine_codon_table_IisG.txt";
    let file_path = "inosine_codon_table1.txt";
    let mut reader = ReaderBuilder::new().delimiter(b'\t').has_headers(false).from_path(file_path).unwrap();

    for result in reader.records() {
        let record = result.unwrap();
        if record.len() < 1 {
            continue;
        }
        let codon = record[0].to_string();
        let translation = record[1].to_string();
        substituted_codon_lookup_table.insert(codon, translation);
    }

    return substituted_codon_lookup_table;

}
