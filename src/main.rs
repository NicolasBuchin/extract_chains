use clap::{Parser, ValueHint};
use plotters::{
    chart::{ChartBuilder, SeriesLabelPosition},
    prelude::{BitMapBackend, Cross, IntoDrawingArea, PathElement},
    series::{LineSeries, PointSeries},
    style::{
        full_palette::{ORANGE, PURPLE},
        Color, BLACK, BLUE, GREEN, RED, WHITE,
    },
};
use rayon::prelude::*;
use std::{
    fs::{create_dir_all, read_to_string},
    io::Result,
    path::Path,
};

#[derive(Debug, Clone)]
struct Anchor {
    ref_start: u32,
    query_start: u32,
}

#[derive(Debug, Clone)]
struct Chain {
    ref_id: u32,
    score: f64,
    qspan: [u32; 2],
    rspan: [u32; 2],
    is_revcomp: bool,
    anchors: Vec<Anchor>,
    cigar: String,
    ref_start: u32,
    considered: bool,
    ssw_cigar: String,
    ssw_ref_start: u32,
}

#[derive(Debug)]
struct Read {
    name: String,
    read_len: u32,
    k: u32,
    fwd_anchors: Vec<Anchor>,
    rev_anchors: Vec<Anchor>,
    chains: Vec<Chain>,
}

#[derive(Parser, Debug)]
struct Args {
    #[arg(value_hint = ValueHint::FilePath)]
    file: String,

    #[arg(short = 'n')]
    n: Option<usize>,

    #[arg(short = 'o', default_value = "plots")]
    output: String,

    #[arg(short = 'x')]
    mapping_only: bool,
}

fn parse_anchors(bytes: &[u8], i: &mut usize) -> Vec<Anchor> {
    let mut anchors = Vec::new();
    while bytes[*i] != b']' {
        *i += 1;
        let start = *i;
        while bytes[*i] != b',' {
            *i += 1;
        }
        let ref_start = std::str::from_utf8(&bytes[start..*i])
            .unwrap()
            .parse()
            .unwrap();
        *i += 1;
        let start = *i;
        while bytes[*i] != b'}' {
            *i += 1;
        }
        let query_start = std::str::from_utf8(&bytes[start..*i])
            .unwrap()
            .parse()
            .unwrap();
        *i += 1;

        anchors.push(Anchor {
            ref_start,
            query_start,
        });
    }
    anchors
}

fn parse_chains(bytes: &[u8], i: &mut usize) -> Vec<Chain> {
    let mut chains = Vec::new();
    while bytes[*i] != b']' {
        *i += 8;
        let start = *i;
        while bytes[*i] != b',' {
            *i += 1;
        }
        let ref_id = std::str::from_utf8(&bytes[start..*i])
            .unwrap()
            .parse()
            .unwrap();
        *i += 7;
        let start = *i;
        while bytes[*i] != b',' {
            *i += 1;
        }
        let score = std::str::from_utf8(&bytes[start..*i])
            .unwrap()
            .parse()
            .unwrap();
        *i += 13;
        let start = *i;
        while bytes[*i] != b',' {
            *i += 1;
        }
        let query_start = std::str::from_utf8(&bytes[start..*i])
            .unwrap()
            .parse()
            .unwrap();
        *i += 11;
        let start = *i;
        while bytes[*i] != b',' {
            *i += 1;
        }
        let query_end = std::str::from_utf8(&bytes[start..*i])
            .unwrap()
            .parse()
            .unwrap();
        *i += 11;
        let start = *i;
        while bytes[*i] != b',' {
            *i += 1;
        }
        let ref_start = std::str::from_utf8(&bytes[start..*i])
            .unwrap()
            .parse()
            .unwrap();
        *i += 9;
        let start = *i;
        while bytes[*i] != b',' {
            *i += 1;
        }
        let ref_end = std::str::from_utf8(&bytes[start..*i])
            .unwrap()
            .parse()
            .unwrap();
        *i += 12;
        let start = *i;
        while bytes[*i] != b',' {
            *i += 1;
        }
        let is_revcomp = std::str::from_utf8(&bytes[start..*i])
            .unwrap()
            .parse()
            .unwrap();
        *i += 10;
        let anchors = parse_anchors(bytes, i);
        *i += 2;

        chains.push(Chain {
            ref_id,
            score,
            qspan: [query_start, query_end],
            rspan: [ref_start, ref_end],
            is_revcomp,
            anchors,
            cigar: "".to_owned(),
            ref_start: 0,
            considered: false,
            ssw_cigar: "".to_owned(),
            ssw_ref_start: 0,
        });
    }
    chains
}

fn parse_cigars(bytes: &[u8], i: &mut usize, chains: &mut [Chain]) {
    let mut n = 0;
    while bytes[*i] != b']' {
        *i += 1;
        let start = *i;
        while bytes[*i] != b',' {
            *i += 1;
        }
        let cigar = std::str::from_utf8(&bytes[start..*i])
            .unwrap()
            .parse()
            .unwrap();
        *i += 16;
        let considered = bytes[*i] == b'1';
        *i += 9;
        let start = *i;
        while bytes[*i] != b',' {
            *i += 1;
        }
        let ref_start = std::str::from_utf8(&bytes[start..*i])
            .unwrap()
            .parse()
            .unwrap();
        *i += 5;
        let start = *i;
        while bytes[*i] != b',' {
            *i += 1;
        }
        let ssw_cigar = std::str::from_utf8(&bytes[start..*i])
            .unwrap()
            .parse()
            .unwrap();
        *i += 12;
        let start = *i;
        while bytes[*i] != b')' {
            *i += 1;
        }
        let ssw_ref_start = std::str::from_utf8(&bytes[start..*i])
            .unwrap()
            .parse()
            .unwrap();
        *i += 1;
        chains[n].cigar = cigar;
        chains[n].ref_start = ref_start;
        chains[n].ssw_cigar = ssw_cigar;
        chains[n].ssw_ref_start = ssw_ref_start;
        chains[n].considered = considered;
        n += 1;
    }
}

fn parse_reads(bytes: &[u8], i: &mut usize, mapping_only: bool) -> Option<Read> {
    *i += 7;
    let start = *i;
    while bytes[*i] != b'\n' {
        *i += 1;
    }
    let name = std::str::from_utf8(&bytes[start..*i])
        .unwrap()
        .parse()
        .unwrap();
    *i += 3;
    let start = *i;
    while bytes[*i] != b',' {
        *i += 1;
    }
    let read_len = std::str::from_utf8(&bytes[start..*i])
        .unwrap()
        .parse()
        .unwrap();
    *i += 3;
    let start = *i;
    while bytes[*i] != b'\n' {
        *i += 1;
    }
    let k = std::str::from_utf8(&bytes[start..*i])
        .unwrap()
        .parse()
        .unwrap();
    *i += 29;
    let fwd_anchors = parse_anchors(bytes, i);
    *i += 30;
    let rev_anchors = parse_anchors(bytes, i);
    *i += 9;
    let mut chains = parse_chains(bytes, i);
    *i += 2;
    if chains.is_empty() {
        return None;
    }

    if mapping_only {
        for (idx, chain) in chains.iter_mut().enumerate() {
            chain.considered = idx == 0;
        }
    } else {
        *i += 8;
        parse_cigars(bytes, i, &mut chains);
    }

    Some(Read {
        name,
        read_len,
        k,
        fwd_anchors,
        rev_anchors,
        chains,
    })
}

fn parse_file(f: &str, n: Option<usize>, mapping_only: bool) -> Vec<Read> {
    let bytes = f.as_bytes();
    let mut reads = Vec::new();

    let mut i = 0;

    while i + 7 < bytes.len() {
        if &bytes[i..i + 7] == b"Query: " {
            if let Some(read) = parse_reads(bytes, &mut i, mapping_only) {
                reads.push(read);
                if let Some(max) = n {
                    if reads.len() >= max {
                        break;
                    }
                }
            }
        } else {
            i += 1;
        }
    }

    println!("parsed {} reads", reads.len());
    reads
}

fn sanitize_filename<S: AsRef<str>>(name: S) -> String {
    name.as_ref()
        .chars()
        .map(|c| {
            if c.is_ascii_alphanumeric() || c == '-' || c == '_' {
                c
            } else {
                '_'
            }
        })
        .collect()
}

fn plot_reads(reads: Vec<Read>, output: &str, mapping_only: bool) {
    create_dir_all(output).unwrap();

    reads.par_iter().for_each(|read| {
        let safe_name = sanitize_filename(&read.name);
        let read_dir = Path::new(output).join(safe_name);
        create_dir_all(&read_dir).unwrap();
        read.chains
            .par_iter()
            .enumerate()
            .for_each(|(chain_idx, chain)| {
                plot_chain(read, chain, chain_idx, &read_dir, mapping_only)
            });
    });
}

fn parse_cigar_to_path(cigar: &str, ref_start: u32) -> Vec<(u32, u32)> {
    let mut path = Vec::new();
    let mut ref_pos = ref_start;
    let mut query_pos = 0u32;

    path.push((ref_pos, query_pos));

    let mut i = 0;
    let chars: Vec<char> = cigar.chars().collect();

    while i < chars.len() {
        let mut num_str = String::new();
        while i < chars.len() && chars[i].is_ascii_digit() {
            num_str.push(chars[i]);
            i += 1;
        }

        if i >= chars.len() {
            break;
        }

        let count: u32 = num_str.parse().unwrap_or(0);
        let operation = chars[i];
        i += 1;

        match operation {
            'M' | '=' | 'X' => {
                ref_pos += count;
                query_pos += count;
                path.push((ref_pos, query_pos));
            }
            'I' => {
                query_pos += count;
                path.push((ref_pos, query_pos));
            }
            'D' => {
                ref_pos += count;
                path.push((ref_pos, query_pos));
            }
            'S' => {
                query_pos += count;
                path.push((ref_pos, query_pos));
            }
            'N' => {
                ref_pos += count;
                path.push((ref_pos, query_pos));
            }
            _ => {}
        }
    }

    path
}

fn plot_chain(read: &Read, chain: &Chain, chain_idx: usize, read_dir: &Path, mapping_only: bool) {
    let ref_start = chain.rspan[0];
    let ref_end = chain.rspan[1];
    let padding = read.read_len / 10;

    let (ref_plot_start, ref_plot_end) = if mapping_only {
        let ref_plot_start = ref_start.saturating_sub(padding);
        let ref_plot_end = ref_end + padding;
        (ref_plot_start, ref_plot_end)
    } else {
        let min_ref_start = ref_start.min(chain.ssw_ref_start);
        let max_ref_start = ref_start.max(chain.ssw_ref_start);
        let ref_plot_start = min_ref_start.saturating_sub(padding);
        let ref_plot_end = max_ref_start + padding + read.read_len;
        (ref_plot_start, ref_plot_end)
    };

    let filename = format!("chain_score={}_cigar={:.2}.png", chain_idx, chain.score);
    let filepath = read_dir.join(filename.clone());

    let root = BitMapBackend::new(&filepath, (1600, 1600)).into_drawing_area();
    root.fill(&WHITE).unwrap();

    let title = format!(
        "Score: {:.2}, Ref ID: {}, Ref Span: {}-{}",
        chain.score, chain.ref_id, ref_start, ref_end
    );

    let mut chart = ChartBuilder::on(&root)
        .caption(&title, ("Arial", 20))
        .margin(50)
        .x_label_area_size(60)
        .y_label_area_size(40)
        .build_cartesian_2d(ref_plot_start..ref_plot_end, 0u32..read.read_len)
        .unwrap();

    chart
        .configure_mesh()
        .x_desc("Reference")
        .y_desc("Query")
        .draw()
        .unwrap();

    let anchors_to_plot = if chain.is_revcomp {
        &read.rev_anchors
    } else {
        &read.fwd_anchors
    };

    let filtered_anchors: Vec<&Anchor> = anchors_to_plot
        .iter()
        .filter(|anchor| {
            anchor.ref_start >= ref_plot_start && anchor.ref_start + read.k <= ref_plot_end
        })
        .collect();

    for anchor in &filtered_anchors {
        let query_end = anchor.query_start + read.k;
        let ref_end = anchor.ref_start + read.k;

        chart
            .draw_series(LineSeries::new(
                vec![(anchor.ref_start, anchor.query_start), (ref_end, query_end)],
                BLUE.stroke_width(2),
            ))
            .unwrap();

        chart
            .draw_series(PointSeries::of_element(
                vec![(anchor.ref_start, anchor.query_start)],
                10,
                &BLUE,
                &|c, s, st| Cross::new(c, s, st.filled()),
            ))
            .unwrap();

        chart
            .draw_series(PointSeries::of_element(
                vec![(ref_end, query_end)],
                10,
                &BLUE,
                &|c, s, st| Cross::new(c, s, st.filled()),
            ))
            .unwrap();
    }

    let chain_color = if chain.considered {
        GREEN.mix(0.5)
    } else {
        RED.mix(0.5)
    };

    for anchor in &chain.anchors {
        let query_end = anchor.query_start + read.k;
        let ref_end = anchor.ref_start + read.k;

        chart
            .draw_series(LineSeries::new(
                vec![(anchor.ref_start, anchor.query_start), (ref_end, query_end)],
                chain_color.stroke_width(4),
            ))
            .unwrap();
    }

    for i in 0..chain.anchors.len().saturating_sub(1) {
        let current_anchor = &chain.anchors[i];
        let next_anchor = &chain.anchors[i + 1];

        let current_end_query = current_anchor.query_start + read.k;
        let current_end_ref = current_anchor.ref_start + read.k;

        chart
            .draw_series(LineSeries::new(
                vec![
                    (current_end_ref, current_end_query),
                    (next_anchor.ref_start, next_anchor.query_start),
                ],
                chain_color.stroke_width(4),
            ))
            .unwrap();
    }

    if !mapping_only {
        let piecewise_path = parse_cigar_to_path(&chain.cigar, chain.ref_start);
        if piecewise_path.len() > 1 {
            chart
                .draw_series(LineSeries::new(
                    piecewise_path.clone(),
                    PURPLE.mix(0.5).stroke_width(4),
                ))
                .unwrap();
        }

        let ssw_path = parse_cigar_to_path(&chain.ssw_cigar, chain.ssw_ref_start);
        if ssw_path.len() > 1 {
            chart
                .draw_series(LineSeries::new(
                    ssw_path.clone(),
                    ORANGE.mix(0.5).stroke_width(4),
                ))
                .unwrap();
        }
    }

    chart
        .draw_series(std::iter::once(PathElement::new(
            [(ref_plot_start, 0), (ref_plot_start + 1, 0)],
            BLUE,
        )))
        .unwrap()
        .label("Blue: Background anchors")
        .legend(|(x, y)| PathElement::new([(x, y), (x + 20, y)], BLUE));

    let chain_label = format!(
        "{}: Chain (considered: {})",
        if chain.considered { "Green" } else { "Red" },
        chain.considered
    );

    chart
        .draw_series(std::iter::once(PathElement::new(
            [(ref_plot_start, 0), (ref_plot_start + 1, 0)],
            chain_color,
        )))
        .unwrap()
        .label(&chain_label)
        .legend(move |(x, y)| PathElement::new([(x, y), (x + 30, y)], chain_color.stroke_width(4)));

    if !mapping_only {
        let ssw_label = format!("Orange: SSW path: {}", chain.ssw_cigar);
        chart
            .draw_series(std::iter::once(PathElement::new(
                [(ref_plot_start, 0), (ref_plot_start + 1, 0)],
                ORANGE.mix(0.5),
            )))
            .unwrap()
            .label(&ssw_label)
            .legend(|(x, y)| {
                PathElement::new([(x, y), (x + 30, y)], ORANGE.mix(0.5).stroke_width(4))
            });

        let piecewise_label = format!("Purple: Piecewise path: {}", chain.cigar);
        chart
            .draw_series(std::iter::once(PathElement::new(
                [(ref_plot_start, 0), (ref_plot_start + 1, 0)],
                PURPLE.mix(0.5),
            )))
            .unwrap()
            .label(&piecewise_label)
            .legend(|(x, y)| {
                PathElement::new([(x, y), (x + 30, y)], PURPLE.mix(0.5).stroke_width(4))
            });
    }

    chart
        .configure_series_labels()
        .background_style(WHITE.mix(0.9))
        .border_style(BLACK)
        .label_font(("Arial", 22))
        .position(SeriesLabelPosition::UpperLeft)
        .draw()
        .unwrap();

    root.present().unwrap();

    println!("{}", filename);
}

fn main() -> Result<()> {
    let args = Args::parse();
    let file = read_to_string(args.file)?;
    let reads = parse_file(&file, args.n, args.mapping_only);
    plot_reads(reads, &args.output, args.mapping_only);
    Ok(())
}
