use bio::alignment::pairwise::*;
use bio::alignment::AlignmentOperation::*;
use bio::alphabets;
use bio::scores::blosum62;
use pyo3::prelude::*;

// Tools to quantify methylation in reduced representation bisulfite sequencing reads.

// Adapted from QUMA CLI: http://quma.cdb.riken.jp/
// Licensed under GPLv3
// initial Perl conversion by: https://freelancer.com/u/Zubayerskd
// refactoring and Rust conversion by Paradoxdruid

// See https://docs.rs/bio/latest/bio/alignment/pairwise/index.html

// struct of quma aligment comparison results
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
struct Fasta {
    com: String,
    pos: String,
    seq: String,
}

// struct of quma analysis intermediates.
// includes fasta sequence, quma results, directon of read, genomic direction,
// and whether result meets exclusion criteria.
struct Reference {
    fasta: Fasta,
    res: QumaResult,
    dir: i32,
    gdir: i32,
    exc: i32,
}

// Quma methylation analysis parser for bisulfite conversion DNA sequencing.
struct Quma {
    gfile_contents: String,
    qfile_contents: String,
    gseq: String,
    qseq: String,
    gfilep_f: String,
    data: Vec<Reference>,
    values: String,
}

impl Quma {
    fn new(
        gfile_contents: String,
        qfile_contents: String,
        gseq: String,
        qseq: String,
        gfilep_f: String,
    ) -> Quma {
        gfilep_f = Quma::fasta_make(gseq, String::from("genomeF"));
        data: Vec<Reference> = Quma::process_fasta_output(
            qseq,
            String::from("queryF"),
            String::from("queryR"),
            gfilep_f,
        );
        values = Quma::format_output(gseq, data);
        return Quma(
            gfile_contents = gfile_contents,
            qfile_contents = qfile_contents,
            gseq = gseq,
            qseq = qseq,
            gfilep_f = gfilep_f,
            data = data,
            values = values,
        );
    }

    fn parse_genome(&self) -> String {
        seq = self.gfile_contents.clone();
        return seq;
    }

    fn scrub_whitespace(string: String) -> String {
        return string;
    }

    fn parse_biseq(&self) -> Vec<Fasta> {
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

    fn align_seq_and_generate_stats(&self, gfile: String, qfile: String) -> Result {
        return ();
    }

    fn process_alignment_matches(&self, result: Result) -> Result {
        return ();
    }

    fn generate_summary_stats(&self, result: Result) -> Result {
        return ();
    }

    fn percentage(a: i32, b: i32, calc_type: String) -> f32 {
        return ();
    }

    fn find_best_dataset(ffres: Result, frres: Result) -> (Result, i32) {
        return ();
    }

    fn format_output(gseq: String, data: Vec<Reference>) -> String {
        return ();
    }
}

/// Run quma and return the quma object
#[pyfunction]
fn quma(seq: String) -> PyResult<String> {
    Ok(Quma::new())
}

/// A Python module implemented in Rust.
#[pymodule]
fn rust_quma(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(quma, m)?)?;
    Ok(())
}
