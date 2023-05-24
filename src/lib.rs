use bio::alignment::pairwise::*;
use bio::alignment::Alignment;
use bio::alignment::AlignmentOperation::*;
use bio::alphabets;
use bio::scores::blosum62;
use lazy_static::lazy_static;
use pyo3::prelude::*;
use regex::Regex;
use std::collections::HashMap;

type EnvMap = HashMap<String, String>;

// Tools to quantify methylation in reduced representation bisulfite sequencing reads.

// Adapted from QUMA CLI: http://quma.cdb.riken.jp/
// Licensed under GPLv3
// initial Perl conversion by: https://freelancer.com/u/Zubayerskd
// refactoring and Rust conversion by Paradoxdruid

// See https://docs.rs/bio/latest/bio/alignment/pairwise/index.html

lazy_static! {
    static ref ALPHABET: String = "ACGTURYMWSKDHBVNacgturymwskdhbvn".to_string();
}

// struct of quma aligment comparison results
#[pyclass]
struct QumaResult {
    q_ali: String,
    g_ali: String,
    val: String,
    perc: f32,
    pconv: f32,
    gap: i32,
    menum: i32,
    unconv: i32,
    conv: i32,
    quma_match: i32,
    ali_mis: i32,
    ali_len: i32,
}

// struct to to wrap fasta results
#[pyclass]
struct Fasta {
    com: String,
    pos: String,
    seq: String,
}

// struct of quma analysis intermediates.
// includes fasta sequence, quma results, directon of read, genomic direction,
// and whether result meets exclusion criteria.
#[pyclass]
struct Reference {
    fasta: Fasta,
    res: QumaResult,
    dir: i32,
    gdir: i32,
    exc: i32,
}

// Quma methylation analysis parser for bisulfite conversion DNA sequencing.
#[pyclass]
struct Quma {
    gfile_contents: String,
    qfile_contents: String,
    gseq: String,
    qseq: Vec<Fasta>,
    gfilep_f: String,
    data: Vec<Reference>,
    values: String,
}

/// Create new Quma struct
///
/// # Arguments
///
/// * `gfile_contents` - genome fasta file contents
/// * `qfile_contents` - query fasta file contents
///
/// # Returns
///
/// * `Quma` - Quma struct
fn py_new(gfile_contents: String, qfile_contents: String) -> Quma {
    let gseq = parse_genome(gfile_contents);
    let qseq = parse_biseq(qfile_contents);
    let gfilep_f = fasta_make(gseq, String::from("genomeF"));
    let data: Vec<Reference> = process_fasta_output(
        qseq,
        String::from("queryF"),
        String::from("queryR"),
        gfilep_f,
    );
    let values = format_output(gseq, data);
    return Quma {
        gfile_contents: gfile_contents,
        qfile_contents: qfile_contents,
        gseq: gseq,
        qseq: qseq,
        gfilep_f: gfilep_f,
        data: data,
        values: values,
    };
}

/// Parse genome file, removing white spaces and extra returns.
///
/// # Arguments
///
/// * `self` - Quma struct
///
/// # Returns
///
/// * `string` - parsed and curated string of genome sequence
fn parse_genome(gfile_contents: String) -> String {
    let re_one = Regex::new(r"^[\r\s]+").unwrap();
    let out_one = re_one.replace_all(&gfile_contents, "");
    let re_two = Regex::new(r"[\r\s]+$").unwrap();
    let out_two = re_two.replace_all(&out_one, "");
    let re_three = Regex::new(r"(\r|\n|\r\n){2}").unwrap();
    let out_three = re_three.replace_all(&out_two, r"\r|\n|\r\n");

    let out_seq = parse_seq(out_three.to_string());
    return out_seq;
}

/// Remove whitespace and repeated newlines.
///
/// # Arguments
///
/// * `string` - input fasta string
///
/// # Returns
///
/// * `string` - fasta string with whitespace removed
fn scrub_whitespace(string: String) -> String {
    let trim_one = string.trim();
    let re_one = Regex::new(r"\r\n\r\n").unwrap();
    let trim_two = re_one.replace_all(&trim_one, r"\r\n");
    let re_two = Regex::new(r"\n\n").unwrap();
    let trim_three = re_two.replace_all(&trim_two, r"\n");
    let re_three = Regex::new(r"\r\r").unwrap();
    let trim_four = re_three.replace_all(&trim_three, r"\r");
    let re_four = Regex::new(r"\r\n").unwrap();
    let trim_five = re_four.replace_all(&trim_four, r"\n");
    let re_five = Regex::new(r"\r").unwrap();
    let trim_six = re_five.replace_all(&trim_five, r"\n");
    return trim_six.to_string();
}

/// Parse bisulfite sequencing fasta file
///
/// # Returns
///
/// * `vector` - vector of Fasta structs of sequence reads
fn parse_biseq(qfile_contents: String) -> Vec<Fasta> {
    let multi_clean = scrub_whitespace(qfile_contents);
    let mut multi_split = multi_clean.lines();

    let biseq = Vec::<Fasta>::new();
    for line in multi_split {
        let mut fa = Fasta {
            com: String::from(""),
            pos: String::from(""),
            seq: String::from(""),
        };
        if line.contains(">") {
            let re_one = Regex::new(r"^>").unwrap();
            let trim_one = re_one.replace_all(&line, "");
            let re_two = Regex::new(r"\s*$").unwrap();
            let trim_two = re_two.replace_all(&trim_one, "");
            fa.com = trim_two.to_string();
            biseq.push(fa);
        } else {
            let allowed = check_char_in_allowed(line.to_string(), ALPHABET.to_string());
            if allowed == "" {
                continue;
            }
            let seq = allowed.to_uppercase();
            fa.seq = seq;
        }
    }
    return biseq;
}

/// Extract sequence strings from the string of a text file
///
/// # Arguments
///
/// * `seq` - string of text file
///
/// # Returns
///
/// * `string` - sequence string
fn parse_seq(seq: String) -> String {
    let re_one = Regex::new(r"\r\n").unwrap();
    let seq_one = re_one.replace_all(&seq, r"\n").to_string();
    let re_two = Regex::new(r"\r").unwrap();
    let seq_two = re_two.replace_all(&seq_one, r"\n").to_string();
    let seq_three = seq_two.to_uppercase();

    // Different file patterns
    let re_three = Regex::new(r"^\s*>.*?\n").unwrap();
    let seq_four = re_three.replace_all(&seq_three, r"").to_string();

    // TODO: Did not implement unused fasta patterns

    let final_seq = check_char_in_allowed(seq_four, ALPHABET.to_string());

    return final_seq;
}

/// Return only charcters in string present in pattern
///
/// # Arguments
///
/// * `seq` - sequence string
/// * `pattern` - string of allowed characters
///
/// # Returns
///
/// * `string` - sequence string with only allowed characters
fn check_char_in_allowed(seq: String, pattern: String) -> String {
    // let mut seq_out = seq.clone();
    // let pattern_chars = pattern.chars();
    // for each in seq_out.chars() {
    //     if pattern_chars.contains(&each) {
    //         continue;
    //     } else {
    //         seq_out = seq_out.replace(&each, "");
    //     }
    // }
    // return seq_out;
    return ();
}

/// Write a sequence string to a fasta-formatted text file contents
///
/// # Arguments
///
/// * `seq` - sequence string
/// * `seq_name` - name of sequence
///
/// # Returns
///
/// * `string` - fasta-formatted text file contents
fn fasta_make(seq: String, seq_name: String) -> String {
    let re_one = Regex::new(r"[0-9]| |\t|\n|\r|\f").unwrap();
    let seq_one = re_one.replace_all(&seq, "").to_string();
    let re_two = Regex::new(r"(.{1,60})").unwrap();
    // hardcode 60 as nothing else is ever passed
    let seq_two = re_two.replace_all(&seq_one, r"\1\n").to_string();

    let final_str = format!(">{}\n{}", seq_name, seq_two);
    return final_str;
}

/// Process fasta alignment
///
/// # Arguments
///
/// * `qseq` - vector of Fasta structs of query sequence
/// * `qfile_f` - query sequence forward read
/// * `qfile_r` - query sequence reverse complement
/// * `gfilep_f` - genome sequence forward read
///
/// # Returns
///
/// * `vector` - vector of Reference structs
fn process_fasta_output(
    qseq: Vec<Fasta>,
    qfile_f: String,
    qfile_r: String,
    gfilep_f: String,
) -> Vec<Reference> {
    let unconv = 5;
    let pconv = 95.0;
    let mis = 10;
    let perc = 90.0;

    let mut data = Vec::<Reference>::new();
    let mut pos = 0;
    for fa in qseq {
        pos += 1;
        fa.pos = pos.to_string();

        let qfile_f_processed = fasta_make(fa.seq.clone(), qfile_f.clone());
        let qfile_r_processed = fasta_make(rev_comp(fa.seq.clone()), qfile_r.clone());

        let fwd_result = align_seq_and_generate_stats(qfile_f_processed, gfilep_f.clone());
        let rev_result = align_seq_and_generate_stats(qfile_r_processed, gfilep_f.clone());

        let (this_result, final_direction) = find_best_dataset(fwd_result, rev_result);

        let genome_direction = 1;

        let this_ref = Reference {
            fasta: fa,
            res: this_result,
            dir: final_direction,
            gdir: genome_direction,
            exc: 0,
        };

        if this_result.unconv > unconv {
            this_ref.exc = 1;
        } else if this_result.pconv > pconv {
            this_ref.exc = 1;
        } else if this_result.ali_mis > mis {
            this_ref.exc = 1;
        } else if this_result.perc > perc {
            this_ref.exc = 1;
        }

        data.push(this_ref);
    }

    return data;
}

/// Return reverse complement of sequence
///
/// # Arguments
///
/// * `seq` - sequence string
///
/// # Returns
///
/// * `string` - reverse complement of sequence
fn rev_comp(seq: String) -> String {
    let reversed = seq.chars().rev().collect::<String>();

    fn expand_env(value: &str, env_user: &EnvMap) -> String {
        let mut expanded = value.to_string();
        for (env_key, env_value) in env_user {
            expanded = expanded.replace(env_key, env_value);
        }
        expanded.to_string()
    }

    // reverse comp mapping
    let mut mappings = HashMap::from([
        ("A".to_string(), "T".to_string()),
        ("C".to_string(), "G".to_string()),
        ("G".to_string(), "C".to_string()),
        ("T".to_string(), "A".to_string()),
        ("U".to_string(), "A".to_string()),
        ("R".to_string(), "Y".to_string()),
        ("Y".to_string(), "R".to_string()),
        ("M".to_string(), "K".to_string()),
        ("W".to_string(), "W".to_string()),
        ("S".to_string(), "S".to_string()),
        ("K".to_string(), "M".to_string()),
        ("D".to_string(), "H".to_string()),
        ("H".to_string(), "D".to_string()),
        ("B".to_string(), "V".to_string()),
        ("V".to_string(), "B".to_string()),
        ("N".to_string(), "N".to_string()),
        ("a".to_string(), "t".to_string()),
        ("c".to_string(), "g".to_string()),
        ("g".to_string(), "c".to_string()),
        ("t".to_string(), "a".to_string()),
        ("u".to_string(), "a".to_string()),
        ("r".to_string(), "y".to_string()),
        ("y".to_string(), "r".to_string()),
        ("m".to_string(), "k".to_string()),
        ("w".to_string(), "w".to_string()),
        ("s".to_string(), "s".to_string()),
        ("k".to_string(), "m".to_string()),
        ("d".to_string(), "h".to_string()),
        ("h".to_string(), "d".to_string()),
        ("b".to_string(), "v".to_string()),
        ("v".to_string(), "b".to_string()),
        ("n".to_string(), "n".to_string()),
    ]);

    mappings.insert("$A".to_string(), "Alpha".to_string());

    let value = expand_env(&reversed, &mappings);

    return value;
}

/// Find pairwise alignment substrings
///
/// # Arguments
///
/// * `alignment` - alignment object
///
/// # Returns
///
/// * `tuple` - tuple of String, String query and genomic aligned substrings
fn matching_substrings(alignment: Alignment) -> (String, String) {
    return ();
}

/// Run pairwise sequence alignment
///
/// # Arguments
///
/// * `gfile` - genomic sequence file contents
/// * `qfile` - sequencing read(s) file contents
///
/// # Returns
///
/// * `QumaResult` - alignment result struct
fn align_seq_and_generate_stats(gfile: String, qfile: String) -> QumaResult {
    let this_result = QumaResult {
        q_ali: "".to_string(),
        g_ali: "".to_string(),
        val: "".to_string(),
        perc: 0.0,
        pconv: 0.0,
        gap: 0,
        menum: 0,
        unconv: 0,
        conv: 0,
        quma_match: 0,
        ali_mis: 0,
        ali_len: 0,
    };

    // TODO: alignment

    return this_result;
}

fn process_alignment_matches(result: QumaResult) -> QumaResult {
    return ();
}

/// Helper to generate summary statistics in QumaResult struct
///
/// # Arguments
///
/// * `result` - QumaResult struct
///
/// # Returns
///
/// * `QumaResult` - QumaResult struct with summary statistics
fn generate_summary_stats(result: QumaResult) -> QumaResult {
    if result.conv + result.unconv != 0 {
        result.pconv = percentage(result.conv, result.unconv, "sum".to_string())
    } else {
        result.pconv = 0.0;
    }

    result.perc = percentage(result.quma_match, result.ali_len, "total".to_string());
    result.ali_mis = result.ali_len - result.quma_match;
    return result;
}

/// Helper to return percentages
///
/// # Arguments
///
/// * `a` - numerator
/// * `b` - denominator
/// * `calc_type` - type of calculation to perform (`sum` or `total`)
///
/// # Returns
///
/// * `f32` - percentage
fn percentage(a: i32, b: i32, calc_type: String) -> f32 {
    if calc_type == "sum".to_string() {
        return 100.0 * a as f32 / (a as f32 + b as f32);
    } else if calc_type == "total".to_string() {
        return 100.0 * a as f32 / b as f32;
    } else {
        return 0.0;
    }
    // TODO: Implement error behavior
}

/// Helper to find best data returned
///
/// # Arguments
///
/// * `ffres` - quma result from forward alignment
/// * `frres` - quma result from reverse alignment
///
/// # Returns
///
/// * `(QumaResult, i32)` - tuple of best QumaResult and direction
fn find_best_dataset(ffres: QumaResult, frres: QumaResult) -> (QumaResult, i32) {
    // FIXME: Find best dataset better

    if ffres.ali_len > frres.ali_len {
        let fres = ffres;
        let fdir = 1;
    } else {
        let fres = frres;
        let fdir = -1;
    }

    return (fres, fdir);
}

fn format_output(gseq: String, data: Vec<Reference>) -> String {
    return ();
}

/// Run quma and return the quma object
// #[pyfunction]
// fn quma(seq: String) -> PyResult<String> {
//     Ok(Quma::new())
// }

/// A Python module implemented in Rust.
#[pymodule]
fn my_module(_py: Python<'_>, m: &PyModule) -> PyResult<()> {
    m.add_class::<Quma>()?;
    Ok(())
}
// fn rust_quma(_py: Python, m: &PyModule) -> PyResult<()> {
//     m.add_function(wrap_pyfunction!(quma, m)?)?;
//     Ok(())
// }
