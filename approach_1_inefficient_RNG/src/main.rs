use csv::ReaderBuilder;
use rand::Rng;
use std::fs::OpenOptions;
use std::collections::HashMap;
use std::io::prelude::*;
use datagen::*;
use annotate::*;
use determine_translation::*;

mod datagen;
mod annotate;
mod determine_translation;

pub struct Codons {
    run_number: u32,
    id: String,
    original_codon: String,
    codon_number: u32,
    substituted_codon: String
}

fn main() {
    
    // parse_run(2);
    process_codons_file();

}


pub fn r_n_g () -> f64 {
    
    // generate a random number between 0 and 1
    // return the number
    // this is a separate function to allow for the RNG to be swapped-out if desired
    
    let mut rng = rand::thread_rng();
    return rng.gen::<f64>();
}
