use rand::Rng;
use datagen::*;
use annotate::*;
use am_lookup::*;
use std::error::Error;

mod datagen;
mod annotate;
mod determine_translation;
mod am_lookup;

pub struct Codons {
    run_number: u32,
    id: String,
    original_codon: String,
    codon_number: u32,
    substituted_codon: String
}

fn main() {
    
    // parse_run(10);
    // process_codons_file();
    am_lookup();

}


pub fn r_n_g () -> f64 {
    
    // generate a random number between 0 and 1
    // return the number
    // this is a separate function to allow for the RNG to be swapped-out if desired
    
    let mut rng = rand::thread_rng();
    return rng.gen::<f64>();
}



