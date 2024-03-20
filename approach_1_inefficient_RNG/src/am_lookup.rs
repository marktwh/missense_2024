use polars::prelude::*;
use std::fs::File;

pub fn am_lookup() {

    let mut changes_df = CsvReader::from_path("missense_changes_toy.tsv")
        .unwrap()
        .has_header(false)
        .with_separator(b'\t')
        .finish()
        .unwrap();

    changes_df.rename("column_1", "run").unwrap();
    changes_df.rename("column_2", "ID").unwrap();
    changes_df.rename("column_3", "change").unwrap();

    let lookup_df = CsvReader::from_path("am_toy.tsv")
        .unwrap()
        .has_header(true)
        .with_separator(b'\t')
        .finish()
        .unwrap();
    
    println!("{:?}", changes_df.head(Some(10)));
    println!("{:?}", lookup_df.head(Some(10)));

    let mut df_left_join = changes_df
        .clone()
        .lazy()
        .join(
            lookup_df.clone().lazy(),
            [col("ID"), col("change")],
            [col("ID"), col("change")],
            JoinArgs::new(JoinType::Left),
        )
        .collect().unwrap();
    println!("{:?}", df_left_join.head(Some(10)));
    let file = File::create("scored_missense.tsv").unwrap();
    let mut writer = CsvWriter::new(file);
    writer.finish(&mut df_left_join).unwrap();
}