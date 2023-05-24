use bio::alignment::pairwise::*;
use bio::alignment::AlignmentOperation::*;
use bio::alphabets;
use bio::scores::blosum62;
use lazy_static::lazy_static;
use pyo3::prelude::*;
use regex::Regex;

// Tools to quantify methylation in reduced representation bisulfite sequencing reads.

// Adapted from QUMA CLI: http://quma.cdb.riken.jp/
// Licensed under GPLv3
// initial Perl conversion by: https://freelancer.com/u/Zubayerskd
// refactoring and Rust conversion by Paradoxdruid

// See https://docs.rs/bio/latest/bio/alignment/pairwise/index.html

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
    qseq: String,
    gfilep_f: String,
    data: Vec<Reference>,
    values: String,
}

#[pymethods]
impl Quma {
    #[new]
    fn py_new(gfile_contents: String, qfile_contents: String) -> PyResult<Self> {
        let gseq = Quma::parse_genome(gfile_contents);
        let qseq = Quma::parse_biseq(qfile_contents);
        let gfilep_f = Quma::fasta_make(gseq, String::from("genomeF"));
        let data: Vec<Reference> = Quma::process_fasta_output(
            qseq,
            String::from("queryF"),
            String::from("queryR"),
            gfilep_f,
        );
        let values = Quma::format_output(gseq, data);
        return Ok(Quma(
            gfile_contents = gfile_contents,
            qfile_contents = qfile_contents,
            gseq = gseq,
            qseq = qseq,
            gfilep_f = gfilep_f,
            data = data,
            values = values,
        ));
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
    fn parse_genome(&self) -> String {
        let seq = self.gfile_contents.clone();
        let re_one = Regex::new(r"^[\r\s]+").unwrap();
        let out_one = re_one.replace_all(seq, "");
        let re_two = Regex::new(r"[\r\s]+$").unwrap();
        let out_two = re_two.replace_all(out_one, "");
        let re_three = Regex::new(r"(\r|\n|\r\n){2}").unwrap();
        let out_three = re_three.replace_all(out_two, r"\r|\n|\r\n");

        let out_seq = self.parse_seq(out_three);
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
        let trim_two = re_one.replace_all(trim_one, r"\r\n");
        let re_two = Regex::new(r"\n\n").unwrap();
        let trim_three = re_two.replace_all(trim_two, r"\n");
        let re_three = Regex::new(r"\r\r").unwrap();
        let trim_four = re_three.replace_all(trim_three, r"\r");
        let re_four = Regex::new(r"\r\n").unwrap();
        let trim_five = re_four.replace_all(trim_four, r"\n");
        let re_five = Regex::new(r"\r").unwrap();
        let trim_six = re_five.replace_all(trim_five, r"\n");
        return trim_six;
    }

    fn parse_biseq(&self) -> Vec<Fasta> {
        let multi = self.qfile_contents.clone();
        let multi_clean = self.scrub_whitespace(multi);
        let mut multi_split = multi_clean.lines();

        for line in multi_split {
            if line.contains(">") {
                let mut fa = Fasta {
                    com: String::from(""),
                    pos: String::from(""),
                    seq: String::from(""),
                };
                let fa: com = line;
            } else {
            }
        }
        return ();
    }

    fn parse_seq(&self, seq: String) -> String {
        return ();
    }

    fn check_char_in_allowed(seq: String, pattern: String) -> String {
        return ();
    }

    fn fasta_make(seq: String, seq_name: String) -> String {
        return ();
    }

    fn process_fasta_output(
        &self,
        qseq: Vec<Fasta>,
        qfile_f: String,
        qfile_r: String,
        gfilep_f: String,
    ) -> Vec<Reference> {
        return ();
    }

    fn rev_comp(seq: String) -> String {
        return ();
    }

    fn matching_substrings(alignment: Bio::Alignment::Pairwise::Alignment) -> (String, String) {
        return ();
    }

    fn align_seq_and_generate_stats(&self, gfile: String, qfile: String) -> QumaResult {
        return ();
    }

    fn process_alignment_matches(&self, result: QumaResult) -> QumaResult {
        return ();
    }

    fn generate_summary_stats(&self, result: QumaResult) -> QumaResult {
        return ();
    }

    fn percentage(a: i32, b: i32, calc_type: String) -> f32 {
        return ();
    }

    fn find_best_dataset(ffres: QumaResult, frres: QumaResult) -> (QumaResult, i32) {
        return ();
    }

    fn format_output(gseq: String, data: Vec<Reference>) -> String {
        return ();
    }
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
