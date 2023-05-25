use bio::alignment::pairwise::*;
use bio::alignment::Alignment;
use lazy_static::lazy_static;
use pyo3::prelude::*;
use regex::Regex;

use std::cmp;
use std::collections::HashMap;
extern crate ndarray;

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

lazy_static! {
    // matrix alphabet:  ATGCSWRYKMBVHDNU
    // see https://docs.rs/bio/latest/src/bio/scores/blosum62.rs.html#89-94
    static ref MATRIX: ndarray::Array2<i32> = ndarray::Array::from_shape_vec(
        (16, 16),
        vec![
            5, -4, -4, -4, -4, 1, 1, -4, -4, 1, -4, -1, -1, -1, -2, -4, -4, 5, -4, 5, -4, 1, -4, 1,
            1, -4, -1, -4, -1, -1, -2, 5, -4, -4, 5, -4, 1, -4, 1, -4, 1, -4, -1, -1, -4, -1, -2,
            -4, -4, -4, -4, 5, 1, -4, -4, 1, -4, 1, -1, -1, -1, -4, -2, -4, -4, -4, 1, 1, -1, -4,
            -2, -2, -2, -2, -1, -1, -3, -3, -1, -4, 1, 1, -4, -4, -4, -1, -2, -2, -2, -2, -3, -3,
            -1, -1, -1, 1, 1, -4, 1, -4, -2, -2, -1, -4, -2, -2, -3, -1, -3, -1, -1, -4, -4, 1, -4,
            1, -2, -2, -4, -1, -2, -2, -1, -3, -1, -3, -1, 1, -4, 1, 1, -4, -2, -2, -2, -2, -1, -4,
            -1, -3, -3, -1, -1, 1, 1, -4, -4, 1, -2, -2, -2, -2, -4, -1, -3, -1, -1, -3, -1, -4,
            -4, -1, -1, -1, -1, -3, -3, -1, -1, -3, -1, -2, -2, -2, -1, -1, -1, -4, -1, -1, -1, -3,
            -1, -3, -3, -1, -2, -1, -2, -2, -1, -4, -1, -1, -4, -1, -3, -1, -3, -1, -3, -1, -2, -2,
            -1, -2, -1, -1, -1, -1, -1, -4, -3, -1, -1, -3, -1, -3, -2, -2, -2, -1, -1, -1, -2, -2,
            -2, -2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -2, -4, 5, -4, -4, -4, 1, -4, 1, 1,
            -4, -1, -4, -1, -1, -2, 5
        ]
    )
    .unwrap();
}

#[inline]
fn lookup(a: u8) -> usize {
    // FIXME: These values are wrong, need to account for 16,16 instead of 27,27
    if a == b'Y' {
        23
    } else if a == b'Z' {
        24
    } else if a == b'X' {
        25
    } else if a == b'*' {
        26
    } else {
        (a - 65) as usize
    }
}

// struct of quma aligment comparison results
#[pyclass]
#[derive(Clone)]
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
#[derive(Clone)]
struct Fasta {
    com: String,
    pos: String,
    seq: String,
}

// struct of quma analysis intermediates.
// includes fasta sequence, quma results, directon of read, genomic direction,
// and whether result meets exclusion criteria.
#[pyclass]
#[derive(Clone)]
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
#[pymethods]
impl Quma {
    #[new]
    fn py_new(gfile_contents: String, qfile_contents: String) -> Self {
        let gseq = parse_genome(gfile_contents.clone());
        let qseq = parse_biseq(qfile_contents.clone());
        let gfilep_f = fasta_make(gseq.clone(), String::from("genomeF"));
        let data: Vec<Reference> = process_fasta_output(
            qseq.clone(),
            String::from("queryF"),
            String::from("queryR"),
            gfilep_f.clone(),
        );
        let values = format_output(gseq.clone(), data.clone());
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

    #[getter]
    fn get_values(&self) -> PyResult<String> {
        Ok(self.values.clone())
    }

    #[getter]
    fn get_data(&self) -> PyResult<Vec<Reference>> {
        Ok(self.data.clone())
    }
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
    // println!("Genome sequence parsed: {}", out_seq);
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
    let multi_split = multi_clean.lines();

    let mut biseq = Vec::<Fasta>::new();
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
fn check_char_in_allowed(mut seq: String, pattern: String) -> String {
    let pattern_chars = pattern;
    seq.retain(|p| {
        if !pattern_chars.contains(*&p) {
            false
        } else {
            true
        }
    });
    return seq;
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
    for mut fa in qseq {
        pos += 1;
        fa.pos = pos.to_string();

        let qfile_f_processed = fasta_make(fa.seq.clone(), qfile_f.clone());
        let qfile_r_processed = fasta_make(rev_comp(fa.seq.clone()), qfile_r.clone());

        let fwd_result = align_seq_and_generate_stats(qfile_f_processed, gfilep_f.clone());
        let rev_result = align_seq_and_generate_stats(qfile_r_processed, gfilep_f.clone());

        let (this_result, final_direction) = find_best_dataset(fwd_result, rev_result);

        let genome_direction = 1;

        let mut this_ref = Reference {
            fasta: fa,
            res: this_result.clone(),
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
fn matching_substrings(alignment: Alignment, bio_gseq: &[u8], bio_qseq: &[u8]) -> (String, String) {
    let g_substring =
        String::from_utf8(bio_gseq[alignment.xstart..alignment.xend].to_vec()).unwrap();
    let q_substring =
        String::from_utf8(bio_qseq[alignment.ystart..alignment.yend].to_vec()).unwrap();

    return (g_substring, q_substring);
}

fn quma_score(a: u8, b: u8) -> i32 {
    let a = lookup(a);
    let b = lookup(b);

    MATRIX[(a, b)]
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
    let mut this_result = QumaResult {
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

    let bio_gseq = gfile.lines().next().unwrap().as_bytes();
    let bio_qseq = qfile.lines().next().unwrap().as_bytes();

    // let mut aligner = Aligner::new(-10, -1, &quma_score);
    let score = |a: u8, b: u8| if a == b { 1i32 } else { -1i32 };
    let mut aligner = Aligner::new(-10, -1, &score);
    // TODO: Custom matrix for CpG
    // See https://docs.rs/bio/latest/src/bio/scores/blosum62.rs.html#89-94

    let bio_alignments = aligner.local(bio_gseq, bio_qseq);

    let (query_ali, genome_ali) = matching_substrings(bio_alignments, bio_gseq, bio_qseq);

    let fh_ = format!(">genome\n{}\n>que\n{}\n", genome_ali, query_ali);

    let fh = fh_.lines();

    for (i, el) in fh.clone().enumerate() {
        if el.contains(">que") {
            this_result.q_ali = fh.clone().nth(i + 1).unwrap().to_string();
        } else if el.contains(">genome") {
            this_result.g_ali = fh.clone().nth(i + 1).unwrap().to_string();
        }
    }

    let re_one = Regex::new(r" ").unwrap();
    this_result.q_ali = re_one.replace_all(&this_result.q_ali, "-").to_string();
    this_result.g_ali = re_one.replace_all(&this_result.g_ali, "-").to_string();

    let final_result = process_alignment_matches(this_result);

    return final_result;
}

/// Helper to implement find method for u8 slices
///
/// # Arguments
///
/// * `haystack` - byte slice to search
/// * `needle` - byte slice to search for
///
/// # Returns
///
/// * `Option<usize>` - index of first match
fn find_subsequence(haystack: &[u8], needle: &[u8]) -> Option<usize> {
    haystack
        .windows(needle.len())
        .position(|window| window == needle)
}

/// Process alignment data to populate results dictionary
///
/// # Arguments
///
/// * `result` - QumaResult struct
///
/// # Returns
///
/// * `QumaResult` - QumaResult struct with populated results dictionary
fn process_alignment_matches(mut result: QumaResult) -> QumaResult {
    let g_ali = result.g_ali.as_bytes();
    let q_ali = result.q_ali.as_bytes();

    result.ali_len = q_ali.len() as i32;

    let mut this_sum = 0;
    let it = q_ali.iter().zip(g_ali.iter());
    for (a, b) in it {
        if a == b {
            this_sum += 1;
        } else if *a as char == 'T' && *b as char == 'C' {
            this_sum += 1;
        }
    }

    result.quma_match = this_sum;

    let g_ali_count = g_ali.iter().filter(|&x| x == &b'-').count();
    let q_ali_count = q_ali.iter().filter(|&x| x == &b'-').count();

    result.gap = cmp::max(
        g_ali_count.try_into().unwrap(),
        q_ali_count.try_into().unwrap(),
    );

    let mut exit_cond = 0;
    let mut i = 0;
    while exit_cond == 0 {
        let q_ali_len = q_ali.len();
        let ni = find_subsequence(&q_ali[i..q_ali_len], b"CG");

        if ni != None {
            let ni_value = ni.unwrap();
            if q_ali[ni_value] as char == 'T' {
                result.quma_match += 1;
                result.unconv += 1;
                result.val += "0";
            } else if q_ali[ni_value] as char == 'C' {
                result.conv += 1;
                result.val += "1";
                result.menum += 1;
            } else {
                result.val += &q_ali[ni_value].to_string();
            }

            i = ni_value + 1;
        } else {
            exit_cond = 1;
        }
    }

    if result.val == "" {
        result.val = "-".to_string();
    }

    let results = generate_summary_stats(result);
    return results;
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
fn generate_summary_stats(mut result: QumaResult) -> QumaResult {
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
        return (fres, fdir);
    } else {
        let fres = frres;
        let fdir = -1;
        return (fres, fdir);
    }
}

/// Process program output into quma-formatted string
///
/// # Arguments
///
/// * `gseq` - genomic sequence
/// * `data` - vector of Reference structs
///
/// # Returns
///
/// * `String` - tabular quma-formatted string
fn format_output(gseq: String, data: Vec<Reference>) -> String {
    let header_output = format!("genome\t0\t{}\t1\t0\n", gseq);

    let mut output_holder: Vec<String> = Vec::new();
    for reference in data {
        output_holder.push(format!("{}\t", reference.fasta.pos));
        output_holder.push(format!("{}\t", reference.fasta.com));
        output_holder.push(format!("{}\t", reference.fasta.seq));
        output_holder.push(format!("{}\t", reference.res.q_ali));
        output_holder.push(format!("{}\t", reference.res.g_ali));
        output_holder.push(format!("{}\t", reference.res.ali_len));
        output_holder.push(format!("{}\t", reference.res.ali_mis));
        output_holder.push(format!("{}\t", reference.res.perc));
        output_holder.push(format!("{}\t", reference.res.gap));
        output_holder.push(format!("{}\t", reference.res.menum));
        output_holder.push(format!("{}\t", reference.res.unconv));
        output_holder.push(format!("{}\t", reference.res.conv));
        output_holder.push(format!("{}\t", reference.res.pconv));
        output_holder.push(format!("{}\t", reference.res.val));
        output_holder.push(format!("{}\t", reference.dir));
        output_holder.push(format!("{}\t", reference.gdir));
        output_holder.push(format!("\n"));
    }
    let joined = output_holder.join("");

    return format!("{}{}", header_output, joined);
}

/// Run quma and return the quma object
// #[pyfunction]
// fn quma(seq: String) -> PyResult<String> {
//     Ok(Quma::new())
// }

/// A Python module implemented in Rust.
#[pymodule]
fn rust_quma(_py: Python<'_>, m: &PyModule) -> PyResult<()> {
    m.add_class::<Quma>()?;
    Ok(())
}
// fn rust_quma(_py: Python, m: &PyModule) -> PyResult<()> {
//     m.add_function(wrap_pyfunction!(quma, m)?)?;
//     Ok(())
// }
