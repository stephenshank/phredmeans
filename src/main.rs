use std::fs;

use bio::io::fastq;
use clap::{arg, Command};
use serde_json::json;

fn main() {
    let matches = Command::new("phreadmeans")
        .version("0.1.0")
        .author("Stephen D. Shank, Ph. D. <sshank@temple.edu>")
        .about("High performance average phred score calculator")
        .arg(arg!(--fastq <VALUE>).required(true))
        .arg(arg!(--json <VALUE>).required(true))
        .get_matches();

    let fastq_file_path = matches.get_one::<String>("fastq").expect("required");
    let json_file_path = matches.get_one::<String>("json").expect("required");

    let fastq_file = std::fs::File::open(fastq_file_path);
    let reader = fastq::Reader::new(fastq_file.unwrap());
    let all_phred_data: Vec<(f32, f32, f32, f32, f32, f32, f32, f32)> = reader
        .records()
        .map(|record| {
            let result = record.unwrap();
            let base_score_pairs = result.seq().iter().zip(result.qual().iter());
            let mut a_total_phred = 0.0;
            let mut a_counts = 0.0;
            let mut c_total_phred = 0.0;
            let mut c_counts = 0.0;
            let mut g_total_phred = 0.0;
            let mut g_counts = 0.0;
            let mut t_total_phred = 0.0;
            let mut t_counts = 0.0;
            let mut adjusted_qual;
            for (base, qual) in base_score_pairs {
                adjusted_qual = *qual as f32 - 33.0;
                if *base == 65 {
                    a_total_phred += adjusted_qual;
                    a_counts += 1.0;
                } else if *base == 67 {
                    c_total_phred += adjusted_qual;
                    c_counts += 1.0;
                } else if *base == 71 {
                    g_total_phred += adjusted_qual;
                    g_counts += 1.0;
                } else if *base == 84 {
                    t_total_phred += adjusted_qual;
                    t_counts += 1.0;
                }
            }
            (
                a_total_phred,
                a_counts,
                c_total_phred,
                c_counts,
                g_total_phred,
                g_counts,
                t_total_phred,
                t_counts,
            )
        })
        .collect();
    let mut a_average_phred = 0.0;
    let mut c_average_phred = 0.0;
    let mut g_average_phred = 0.0;
    let mut t_average_phred = 0.0;
    let mut a_count = 0.0;
    let mut c_count = 0.0;
    let mut g_count = 0.0;
    let mut t_count = 0.0;
    for phred_data in all_phred_data {
        let (
            a_phred,
            a_partial_count,
            c_phred,
            c_partial_count,
            g_phred,
            g_partial_count,
            t_phred,
            t_partial_count,
        ) = phred_data;
        if a_partial_count > 0.0 {
            a_average_phred = a_count * a_average_phred + a_phred;
            a_count += a_partial_count;
            a_average_phred /= a_count;
        }
        if c_partial_count > 0.0 {
            c_average_phred = c_count * c_average_phred + c_phred;
            c_count += c_partial_count;
            c_average_phred /= c_count;
        }
        if g_partial_count > 0.0 {
            g_average_phred = g_count * g_average_phred + g_phred;
            g_count += g_partial_count;
            g_average_phred /= g_count;
        }
        if t_partial_count > 0.0 {
            t_average_phred = t_count * t_average_phred + t_phred;
            t_count += t_partial_count;
            t_average_phred /= t_count;
        }
    }

    let data = json!({
        "a_average_phred": a_average_phred,
        "c_average_phred": c_average_phred,
        "t_average_phred": g_average_phred,
        "g_average_phred": t_average_phred
    });

    let json_str = serde_json::to_string(&data).unwrap();

    fs::write(json_file_path, json_str).expect("Unable to write to file.");
}
