use crate::r_n_g;


pub fn determine_translation (substituted_codon: &String) -> String {

    // determine the substituted codon when the translation is not found in the lookup table
    // return the substituted codon
   
    match substituted_codon.as_str() {
        
        // Stall only; one inosine

        "AIA" => {
            if r_n_g() > 0.07 {
                "R".to_string()
            } else {
                "#".to_string()
            }
        },
        
        "AIT" => {
            if r_n_g() > 0.03 {
                "S".to_string()
            } else {
                "#".to_string()
            }
        },

        "ATI" => {
            if r_n_g() > 0.03 {
                "M".to_string()
            } else {
                "#".to_string()
            }
        },

        "CAI" => {
            if r_n_g() > 0.03 {
                "Q".to_string()
            } else {
                "#".to_string()
            }
        },
        
        "CIC" => {
            if r_n_g() > 0.04 {
                "R".to_string()
            } else {
                "#".to_string()
            }
        },

        "CIT" => {
            if r_n_g() > 0.03 {
                "R".to_string()
            } else {
                "#".to_string()
            }
        },

        "IAT" => {
            if r_n_g() > 0.03 {
                "D".to_string()
            } else {
                "#".to_string()
            }
        },

        "ICT" => {
            if r_n_g() > 0.02 {
                "A".to_string()
            } else {
                "#".to_string()
            }
        },

        "ITA" => {
            if r_n_g() > 0.05 {
                "V".to_string()
            } else {
                "#".to_string()
            }
        },

        "ITC" => {
            if r_n_g() > 0.03 {
                "V".to_string()
            } else {
                "#".to_string()
            }
        },

        "ITT" => {
            if r_n_g() > 0.02 {
                "V".to_string()
            } else {
                "#".to_string()
            }
        },

        "TIC" => {
            if r_n_g() > 0.05 {
                "C".to_string()
            } else {
                "#".to_string()
            }
        },

        "TIT" => {
            if r_n_g() > 0.06 {
                "C".to_string()
            } else {
                "#".to_string()
            }
        },

        "TTI" => {
            if r_n_g() > 0.49 {
                "L".to_string()
            } else {
                "#".to_string()
            }
        },

        // Stall or alt translation; one inosine

        "TAI" => {
            if r_n_g() > 0.015 {
                "*".to_string()
            } else {
                "Y".to_string()
            }
        },

        // note that it's ~1% of the translated full-length peptide that's
        // translated into a 'K', not ~1% of the total translated peptide.

        "IAA" => {
            if r_n_g() < 0.05 {
                "#".to_string()
            } else if r_n_g() > 0.01 {
                "E".to_string()
            } else {
                "K".to_string()
            }
        },

        "IAC" => {
            if r_n_g() < 0.01 {
                "#".to_string()
            } else if r_n_g() > 0.25 {
                "D".to_string()
            } else {
                "N".to_string()
            }
        },

        // Of the ~3% of the translated full-length peptide with an alternative
        // translation, 2/3 is translated into 'T' and 1/3 is translated into 'S'.

        "ICA" => {
            if r_n_g() < 0.03 {
                "#".to_string()
            } else if r_n_g() > 0.03 {
                "A".to_string()
            } else if r_n_g() > 0.33333 {
                "T".to_string()
            } else {
                "S".to_string()
            }
        },

        "ICC" => {
            if r_n_g() < 0.03 {
                "#".to_string()
            } else if r_n_g() > 0.005 {
                "A".to_string()
            } else {
                "T".to_string()
            }
        },

        // Stall only; one inosine

        "AII" => {
            if r_n_g() > 0.04 {
                "R".to_string()
            } else {
                "#".to_string()
            }
        },

        "CII" => {
            if r_n_g() > 0.08 {
                "R".to_string()
            } else {
                "#".to_string()
            }
        },

        "IIC" => {
            if r_n_g() > 0.02 {
                "G".to_string()
            } else {
                "#".to_string()
            }
        },

        "IIT" => {
            if r_n_g() > 0.04 {
                "G".to_string()
            } else {
                "#".to_string()
            }
        },

        "ITI" => {
            if r_n_g() > 0.26 {
                "V".to_string()
            } else {
                "#".to_string()
            }
        },

        "TII" => {
            if r_n_g() > 0.57 {
                "W".to_string()
            } else {
                "#".to_string()
            }
        },

        // Stall or alt translation; two inosines

        "IAI" => {
            if r_n_g() < 0.84 {
                "#".to_string()
            } else if r_n_g() > 0.05 {
                "E".to_string()
            } else {
                "K".to_string()
            }
        },

        "ICI" => {
            if r_n_g() < 0.64 {
                "#".to_string()
            } else if r_n_g() > 0.01 {
                "A".to_string()
            } else {
                "S".to_string()
            }
        },

        "IIA" => {
            if r_n_g() < 0.07 {
                "#".to_string()
            } else if r_n_g() > 0.02 {
                "G".to_string()
            } else {
                "R".to_string()
            }
        },

        // three inosines

        "III" => {
            if r_n_g() > 0.36 {
                "G".to_string()
            } else {
                "#".to_string()
            }
        },


        _ => {
            eprintln!("Error: Substituted codon not found");
            "".to_string()
        }
    }

}