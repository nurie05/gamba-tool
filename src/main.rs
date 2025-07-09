// operon_finder.rs
use std::{cmp::Ordering, collections::{BTreeMap, HashMap, HashSet}, fs::File, io::{BufRead, BufReader, BufWriter, Write}, path::PathBuf};
use std::fmt::Debug;
//use clap::builder::TypedValueParser;
use clap::Parser;
use itertools::Itertools;
use log::{info, error};
use noodles::gtf;
//use noodles::core::Position;
//use noodles::gff::record::attributes::field::Value;

type GeneId = String;

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
    fpkm_val: f32,
    exons: Vec<(u64, u64)>,
    raw_lines: Vec<String>,
}

fn transcripts_inside_op(t1: &Transcript, t2: &Transcript, tolerance: u64, threshold: f32) -> bool {
    t1.start <= t2.start + tolerance && t2.start + tolerance < t1.end + tolerance 
    && t1.end + tolerance >= t2.end && t2.end > t1.start
    && t1.coverage * threshold < t2.coverage
    && (t2.exons.len() > 1 || (t1.coverage * threshold * 10.0 < t2.coverage))
}
fn transcripts_inside(t1: &Transcript, t2: &Transcript, tolerance: u64, threshold: f32) -> bool {
    t1.start <= t2.start + tolerance && t2.start + tolerance < t1.end + tolerance 
    && t1.end + tolerance >= t2.end && t2.end > t1.start
    && t1.coverage > t2.coverage * threshold
    && (t2.exons.len() > 1 || (t1.coverage > t2.coverage * threshold * 10.0 ))
}

fn transcripts_no_overlap(t1: &Transcript, t2: &Transcript, tolerance: u64) -> bool {
    t1.start > t2.end.saturating_sub(tolerance)
}

fn operontrans_overlap(t1: &Transcript, t2: &Transcript, tolerance: u64) -> bool {
    t1.start <= t2.end.saturating_sub(tolerance) && t1.end >= t2.start + tolerance
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
            let gid = record.attributes().get("gene_id".as_ref()).map(|v| v.as_string().unwrap().to_string()).unwrap_or("NA".into()).clone();
            let cov = record.attributes()
                .get("cov".as_ref())
                .and_then(|v| v.as_string())
                .and_then(|s| s.to_string().parse::<f32>().ok())
                .unwrap_or(0.0);
            let fpkm = record.attributes()
                .get("FPKM".as_ref())
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
                fpkm_val: fpkm,
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

    let mut operon_to_genes: Vec<(GeneId, Transcript, &Transcript)> = Vec::new();

    for (chrom, transcripts) in &transcripts_by_chrom {
        info!("Processing chromosome {} ({} transcripts)...", chrom, transcripts.len());
        for container in transcripts {
            let mut contained = Vec::new();
            let mut counter=0;
            for inner in transcripts {
                if container.id == inner.id || container.strand != inner.strand {
                    continue;
                }
                if transcripts_inside_op(container,inner, 250, threshold) {
                    contained.push(inner);
                }
                if transcripts_inside(inner,container, 250, threshold) {
                    counter += 1;
                }
            }
            if contained.len() >= 2 && counter == 0 {
                let mut non_overlapping = Vec::new();
                contained.sort_by(|&e1 , &e2| {
                    let id1 = e1.start;
                    let id2 = e2.start;
                    if id1 > id2 { Ordering::Greater } else if id1 < id2 { Ordering::Less } else { Ordering::Equal }
                });
                for gene in contained {
                    if non_overlapping.last().map_or(true, |last: &&Transcript| transcripts_no_overlap(gene,last,50) ) {
                        non_overlapping.push(gene);
                    } else {
                        let last = non_overlapping.last().unwrap();
                        if gene.fpkm_val > last.fpkm_val || (gene.fpkm_val == last.fpkm_val && gene.exons.len() > non_overlapping.last().unwrap().exons.len()) {
                            non_overlapping.pop();
                            non_overlapping.push(gene);
                        }
                    }
                }

                if non_overlapping.len() >= 2 {
                    for gene in non_overlapping {
                        operon_to_genes.push((container.gene_id.clone(), container.clone(), gene));
                    }
                }
            }
        }
    }

    let mut chr_to_operons: HashMap<(String, String), Vec<(String, &Transcript, &&Transcript)>> = HashMap::new();
    for (op_gene_id, op_id, trans_id) in operon_to_genes.iter() {
        chr_to_operons.entry((op_id.chrom.clone(),op_id.strand.clone())).or_default().push((op_gene_id.clone(), op_id, trans_id));
    }

    let mut overlapping: Vec<(String, Transcript, &Transcript)> = Vec::new();
    let mut seen_transcripts = HashSet::new();
    let mut counter = 1;
    for ((_chrom, _strand), mut op_list) in chr_to_operons {
        println!("Chromosome {} strand {}: {} putative operons", &_chrom, &_strand, &op_list.len());
        op_list.sort_by_key(|(_, p, _)| p.start);
        for (_ , current_op, inner_trans) in op_list {
            if seen_transcripts.contains(&inner_trans.id) {
                continue;
            }
            if let Some((_, last, _ ) ) = overlapping.last() {
                if operontrans_overlap(current_op, last, 250) {
                    overlapping.push((format!("OPRN.{}", counter), current_op.clone(), inner_trans.clone()));
                } else {
                    counter += 1;
                    overlapping.push((format!("OPRN.{}", counter), current_op.clone(), inner_trans.clone()));
                }
            } else {
                overlapping.push((format!("OPRN.{}", counter), current_op.clone(), inner_trans.clone()));
            }
            seen_transcripts.insert(inner_trans.id.clone());
        }
    }

    let mut operon_to_trans: HashMap<String, Vec<(Transcript,&Transcript)>> = HashMap::new();
    for (op_id, operon, inner_trans) in overlapping.iter() {
        operon_to_trans.entry(op_id.clone()).or_default().push((operon.clone(), inner_trans.clone()));
    }
    let mut operon_to_trans_def: Vec<(String, String, String)> = Vec::new();
    let mut operon_ids = HashSet::new();
    let mut gene_ids = HashSet::new();
    let mut operon_gene_map: HashMap<String, Vec<String>> = HashMap::new();

    for (operon_id, transcripts_list) in operon_to_trans {
        let mut non_overlapping_def: Vec<&Transcript> = Vec::new();
        let transcripts_list_ordered = transcripts_list.iter().sorted_by(|(_, e1 ), (_, e2)| {
            let id1 = e1.start;
            let id2 = e2.start;
            let id1_strand = &e1.strand;
            let id2_strand = &e2.strand;
            if id1_strand == id2_strand && id1 > id2 { Ordering::Greater } else if id1_strand == id2_strand && id1 < id2 { Ordering::Less } else { Ordering::Equal }
        });
        for (_, gene) in transcripts_list_ordered {
            if non_overlapping_def.last().map_or(true, |last: &&Transcript| transcripts_no_overlap(gene,last,50) ) {
                non_overlapping_def.push(gene);
            } else if gene.fpkm_val > non_overlapping_def.last().unwrap().fpkm_val {
                non_overlapping_def.pop();
                non_overlapping_def.push(gene);
            }
        } 
        
        if non_overlapping_def.len() >= 2 {
            for (operon, gene) in transcripts_list {
                if non_overlapping_def.iter().any(|&i| i.id == gene.id ) {
                    operon_to_trans_def.push((operon_id.clone(), operon.id.clone(), gene.id.clone()));
                    operon_ids.insert(operon.id.clone());
                    gene_ids.insert(gene.id.clone());
                    operon_gene_map.entry(operon_id.clone()).or_default().push(gene.id.clone());
                }
            }
        }
    }

    let mut tsv_path = out_prefix.clone();
    tsv_path.push_str(&format!("_operons_found_v9.t{:.2}.tsv", threshold));
    let mut tsv_file = BufWriter::new(File::create(&tsv_path)?);
    writeln!(tsv_file, "Operon\tOperonTrans\tContained_transcript")?;
    operon_to_trans_def.sort_by(|(_,e1,_) , (_, e2, _)| {
            let e1_els = e1.split(".").collect::<Vec<_>>();
            let e2_els = e2.split(".").collect::<Vec<_>>();
            let id1 = format!("{}.{}", e1_els[1], e1_els[2]).parse::<f32>().unwrap();
            let id2 = format!("{}.{}", e2_els[1], e2_els[2]).parse::<f32>().unwrap();
            if id1 > id2 { Ordering::Greater } else if id1 < id2 { Ordering::Less } else { Ordering::Equal }
        });
    for (operon_id, operon, inner_trans) in &operon_to_trans_def {
        writeln!(tsv_file, "{}\t{}\t{}", operon_id, operon, inner_trans)?;
    }
    info!("Output written to {}", tsv_path);

    let write_gtf = |filename: &str, ids: &HashSet<String>| -> anyhow::Result<()> {
        let mut file = BufWriter::new(File::create(filename)?);
        let mut ids_ordered = Vec::from_iter(ids);
        ids_ordered.sort_by(|e1, e2| {
            let e1_els = e1.split(".").collect::<Vec<_>>();
            let e2_els = e2.split(".").collect::<Vec<_>>();
            let id1 = e1_els[1].parse::<u32>().unwrap();
            let id2 = e2_els[1].parse::<u32>().unwrap();
            if id1 > id2 { Ordering::Greater } else if id1 < id2 { Ordering::Less } else { Ordering::Equal }
        });

        for id in ids_ordered {
            if let Some(lines) = raw_lines_by_id.get(id) {
                for line in lines {
                    //writeln!(file, "{}", line.replace("\"\"", "\"").replace("\";\"", "\";"))?;
                    writeln!(file, "{}", line.replace("\n", ";"))?;
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
    
    println!("Total number of OPRNs found: {}", &operon_gene_map.keys().len());
    info!("Total number of OPRNs found: {}", &operon_gene_map.keys().len());
    println!("Total number of OpGs found: {}", &operon_to_trans_def.len());
    info!("Total number of OpGs found: {}", operon_to_trans_def.len());

    // Summary
    let mut summary = HashMap::from([
        //("1 gene", 0),
        ("2 genes", 0),
        ("3 genes", 0),
        ("4 genes", 0),
        ("5 genes", 0),
        (">5 genes", 0),
    ]);

    for (_operon, genes) in operon_gene_map {
        match genes.len() {
            //1 => *summary.get_mut("1 genes").unwrap() += 1,
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
