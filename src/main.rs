// operon_finder.rs
use std::{collections::{HashMap, HashSet, BTreeMap}, fs::File, path::PathBuf, io::{BufWriter, Write, BufReader, BufRead}};
use std::fmt::Debug;
use clap::builder::TypedValueParser;
use clap::Parser;
use log::{info, error};
use noodles::gtf;
use noodles::core::Position;
use noodles::gff::record::attributes::field::Value;

#[derive(Parser, Debug)]
#[command(name = "Operon Finder")]
#[command(about = "Detect operons from a GTF file with coverage filtering.", long_about = None)]
struct Args {
    /// Path to the input GTF file
    #[arg(short, long)]
    file: PathBuf,

    /// Coverage threshold multiplier
    #[arg(long, default_value_t = 1.0)]
    threshold: f32,

    /// Output file prefix
    #[arg(short, long)]
    output: Option<String>,

    /// Log file path
    #[arg(long)]
    log: Option<String>,
}

#[derive(Debug, Clone)]
struct Transcript {
    id: String,
    gene_id: String,
    chrom: String,
    start: u64,
    end: u64,
    strand: String,
    coverage: f32,
    exons: Vec<(u64, u64)>,
    raw_lines: Vec<String>,
}

fn transcripts_overlap(t1: &Transcript, t2: &Transcript, tolerance: u64) -> bool {
    t1.start <= t2.start + tolerance && t2.start + tolerance < t1.end + tolerance 
    && t1.end + tolerance >= t2.end && t2.end > t1.start
}

fn operontrans_overlap(t1: &Transcript, t2: &Transcript, tolerance: u64) -> bool {
    t1.start <= t2.end.saturating_sub(tolerance) && t1.end >= t2.start + 250
}

fn main() -> anyhow::Result<()> {
    let args = Args::parse();

    let gtf_path = &args.file;
    let threshold = args.threshold;
    let out_prefix = args.output.clone().unwrap_or_else(|| {
        gtf_path.file_stem().unwrap().to_string_lossy().to_string()
    });
    let log_file = args.log.clone().unwrap_or_else(|| format!("{}_OFv9.log", out_prefix));

    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or("info"))
        .format(move |buf, record| writeln!(buf, "{} - {}", chrono::Local::now().format("%Y-%m-%d %H:%M:%S"), record.args()))
        .target(env_logger::Target::Stdout)
        .init();

    let mut reader = gtf::io::Reader::new(BufReader::new(File::open(gtf_path)?));
    let mut transcripts_by_chrom: BTreeMap<String, Vec<Transcript>> = BTreeMap::new();
    let mut exons_by_transcript: HashMap<String, Vec<(u64, u64)>> = HashMap::new();
    let mut raw_lines_by_id: HashMap<String, Vec<String>> = HashMap::new();

    for result in reader.record_bufs() {
        let record = result?;
        let tid = record.attributes().get("transcript_id".as_ref()).map(|v| v.as_string().unwrap().to_string()).unwrap_or("NA".into()).clone();

        if record.ty() == "transcript" {
            let gid = record.attributes().get("transcript_id".as_ref()).map(|v| v.as_string().unwrap().to_string()).unwrap_or("NA".into()).clone();
            let cov = record.attributes()
                .get("cov".as_ref())
                .and_then(|v| v.as_string())
                .and_then(|s| s.to_string().parse::<f32>().ok())
                .unwrap_or(0.0);
            //let cov = record.attributes().get("cov".as_ref()).and_then(|v| v.parse::<f32>().ok()).unwrap_or(0.0);
            let transcript = Transcript {
                id: tid.clone(),
                gene_id: gid,
                chrom: record.reference_sequence_name().to_string(),
                start: record.start().get() as u64,
                end: record.end().get() as u64,
                strand: format!("{:?}", record.strand()),
                coverage: cov,
                exons: Vec::new(),
                raw_lines: vec![format!("{:?}", record)],
            };
            transcripts_by_chrom.entry(transcript.chrom.clone()).or_default().push(transcript);
            let buf = Vec::new();
            let mut writer = noodles::gff::io::Writer::new(buf);
            writer.write_record(&record).expect("Unable to write GFF record");
            raw_lines_by_id.entry(tid.clone()).or_default().push(String::from_utf8(writer.into_inner())?);
        } else if record.ty() == "exon" {
            let start = record.start().get();
            let end = record.end().get();
            exons_by_transcript.entry(tid.clone()).or_default().push((start as u64, end as u64));
            let buf = Vec::new();
            let mut writer = noodles::gff::io::Writer::new(buf);
            writer.write_record(&record).expect("Unable to write GFF record");
            raw_lines_by_id.entry(tid.clone()).or_default().push(String::from_utf8(writer.into_inner())?);
        }
    }

    for transcripts in transcripts_by_chrom.values_mut() {
        for transcript in transcripts.iter_mut() {
            if let Some(exons) = exons_by_transcript.get(&transcript.id) {
                transcript.exons = exons.clone();
            }
            if let Some(lines) = raw_lines_by_id.get(&transcript.id) {
                transcript.raw_lines = lines.clone();
            }
        }
    }

    let mut operon_to_genes: Vec<(String, String, String)> = Vec::new();
    let mut operon_ids = HashSet::new();
    let mut gene_ids = HashSet::new();
    let mut operon_gene_map: HashMap<String, Vec<String>> = HashMap::new();

    for (chrom, transcripts) in &transcripts_by_chrom {
        info!("Processing chromosome {} with {} transcripts", chrom, transcripts.len());
        for container in transcripts {
            let mut contained = Vec::new();
            for inner in transcripts {
                if container.id == inner.id || container.strand != inner.strand {
                    continue;
                }
                let sufficient_cov = container.coverage * threshold <= inner.coverage;
                let multi_exonic = inner.exons.len() > 1;

                if transcripts_overlap(container,inner, 250)==true && sufficient_cov && multi_exonic {
                    contained.push(inner);
                }
            }

            if contained.len() >= 2 {
                let mut non_overlapping = Vec::new();
                for gene in contained {
                    if non_overlapping.last().map_or(true, |last: &&Transcript| gene.start > last.end.saturating_sub(50)) {
                        non_overlapping.push(gene);
                    } else if gene.coverage > non_overlapping.last().unwrap().coverage {
                        non_overlapping.pop();
                        non_overlapping.push(gene);
                    }
                }

                if non_overlapping.len() >= 2 {
                    for gene in non_overlapping.iter() {
                        operon_to_genes.push((container.gene_id.clone(), container.id.clone(), gene.id.clone()));
                        operon_ids.insert(container.id.clone());
                        gene_ids.insert(gene.id.clone());
                        operon_gene_map.entry(container.gene_id.clone()).or_default().push(gene.id.clone());
                    }
                }
            }
        }
    }

    //I need to add the operon operon overlaping loop to assing same oepron_ID to those oves overlaping in the same strand.

    let mut tsv_path = out_prefix.clone();
    tsv_path.push_str(&format!("_operons_found_v9.t{:.2}.tsv", threshold));
    let mut tsv_file = BufWriter::new(File::create(&tsv_path)?);
    writeln!(tsv_file, "Operon\tOperonTrans\tContained_transcript")?;
    for (operon, op_trans, gene) in &operon_to_genes {
        writeln!(tsv_file, "{}\t{}\t{}", operon, op_trans, gene)?;
    }
    info!("Output written to {}", tsv_path);

    let write_gtf = |filename: &str, ids: &HashSet<String>| -> anyhow::Result<()> {
        let mut file = BufWriter::new(File::create(filename)?);
        for id in ids {
            if let Some(lines) = raw_lines_by_id.get(id) {
                for line in lines {
                    writeln!(file, "{}", line.replace("\"\"", "\"").replace("\";\"", "\";"))?;
                }
            }
        }
        Ok(())
    };

    write_gtf(&format!("{}_Operons_v9.t{:.2}.gtf", out_prefix, threshold), &operon_ids)?;
    write_gtf(&format!("{}_OperonGenes_v9.t{:.2}.gtf", out_prefix, threshold), &gene_ids)?;

    let all_gene_ids: HashSet<String> = raw_lines_by_id
        .keys()
        .filter(|id| !operon_ids.contains(*id) && gene_ids.contains(*id))
        .cloned()
        .collect();
    write_gtf(&format!("{}_OperonGenesALL_v9.t{:.2}.gtf", out_prefix, threshold), &all_gene_ids)?;

    let clean_ids: HashSet<String> = raw_lines_by_id
        .keys()
        .filter(|id| !operon_ids.contains(*id) && !gene_ids.contains(*id))
        .cloned()
        .collect();
    write_gtf(&format!("{}_opCLEAN_v9.t{:.2}.gtf", out_prefix, threshold), &clean_ids)?;

    info!("GTF files written successfully.");

    // Summary
    let mut summary = HashMap::from([
        ("2 genes", 0),
        ("3 genes", 0),
        ("4 genes", 0),
        ("5 genes", 0),
        (">5 genes", 0),
    ]);

    for (_operon, genes) in operon_gene_map {
        match genes.len() {
            2 => *summary.get_mut("2 genes").unwrap() += 1,
            3 => *summary.get_mut("3 genes").unwrap() += 1,
            4 => *summary.get_mut("4 genes").unwrap() += 1,
            5 => *summary.get_mut("5 genes").unwrap() += 1,
            n if n > 5 => *summary.get_mut(">5 genes").unwrap() += 1,
            _ => {},
        }
    }

    println!("Summary of operons by gene number:");
    info!("Summary of operons by gene number:");
    for (category, count) in &summary {
        println!("{}: {}", category, count);
        info!("{}: {}", category, count);
    }

    Ok(())
}
