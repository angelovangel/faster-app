
use std::sync::Arc;
use std::io;
use std::io::prelude::*;

use bio::io::fastq;
use flate2::read::MultiGzDecoder;

use dioxus::prelude::*;
use dioxus::prelude::dioxus_elements::FileEngine;
use dioxus::desktop::{Config, LogicalSize, WindowBuilder};

use indicatif::HumanCount;
use human_repr::HumanCount as HC;

mod modules;

fn main() {
    LaunchBuilder::desktop()
        .with_cfg(
            Config::new()
                .with_window(WindowBuilder::new()
                    .with_focused(true)
                    .with_title("fasterx")
                    .with_inner_size(LogicalSize::new(1200, 700))
                )
        )
        .launch(app)
}
struct UploadedFile {
    name: String,
    reads: u64,
    bases: u64,
    nx: u64,
    gc: String,
    q20: String,
    q30: String
}

async fn decode_reader(bytes: Vec<u8>, filename: &String) -> io::Result<String> {
   if filename.ends_with(".gz") {
       let mut gz = MultiGzDecoder::new(&bytes[..]);
       let mut s = String::new();
       gz.read_to_string(&mut s)?;
       Ok(s)   
   } else {
       let s = String::from_utf8(bytes).unwrap();
       Ok(s)
   }
}

fn maketable(
    entries: Signal<Vec<UploadedFile>>, 
    numbers_type: String,
    treads: Signal<u64>,
    tbases: Signal<u64>
) -> Element {
    rsx! {
        for f in entries.read().iter() {
            tr {
            td {"{f.name}"}
            if numbers_type =="comma" {
                td {"{HumanCount(f.reads)}"}
                td {"{HumanCount(f.bases)}"}
                td {"{HumanCount(f.nx)}"}
            } else if numbers_type == "human" {
                td {"{f.reads.human_count_bare()}"}
                td {"{f.bases.human_count_bare()}"}
                td {"{f.nx.human_count_bare()}"}
            } else {
                td {"{f.reads}"}
                td {"{f.bases}"}
                td {"{f.nx}"}
            }
            td {"{f.gc}"}
            td {"{f.q20}"}
            td {"{f.q30}"}
            }
        }
        tr {
            td {"Total"}
            if numbers_type == "comma" {
                td {"{HumanCount(treads())}"}
                td {"{HumanCount(tbases())}"}
            } else if numbers_type == "human" {
                td {"{treads().human_count_bare()}"}
                td {"{tbases().human_count_bare()}"}
            } else {
                td {"{treads()}"}
                td {"{tbases()}"}
            }
        }
    }
}

fn app() -> Element {
    let mut reads_processed = use_signal(|| 0);
    //let mut enable_directory_upload = use_signal(|| false);
    let mut numbers = use_signal(|| String::new());
    //let mut human_numbers = use_signal(|| false);
    let mut files_uploaded = use_signal(|| Vec::new() as Vec<UploadedFile>);
    let mut total_reads = use_signal(|| 0);
    let mut total_bases = use_signal(|| 0);
    //let mut hovered = use_signal(|| false);
    let mut busy = use_signal(|| false);

    let read_files = move |file_engine: Arc<dyn FileEngine>| async move {
        let files = file_engine.files();
        //to_owned![files];
        for file in &files {
            let mut nreads: u64 = 0;
            let mut nbases: u64 = 0;
            let mut gcbases: u64 = 0;
            let mut qual20: i64 = 0;
            let mut qual30: i64 = 0;
            let mut len_vector: Vec<i64> = Vec::new();

            if let Some(bytes) = file_engine.read_file(&file).await {
                let strings = decode_reader(bytes, file).await.unwrap();
                let mut recs = fastq::Reader::new(strings.as_bytes()).records();
                    
                while let Some(Ok(rec)) = recs.next() {
                    nreads += 1;
                    gcbases += modules::get_gc_bases(rec.seq());
                    reads_processed += 1;
                    nbases += rec.seq().len() as u64;
                    qual20 += modules::get_qual_bases(rec.qual(), 53); // 33 offset
                    qual30 += modules::get_qual_bases(rec.qual(), 63); // 33 offset
                    len_vector.push(rec.seq().len() as i64);   
                }
                let n50 = modules::get_nx(&mut len_vector, 0.5);

                files_uploaded.write().push(UploadedFile {
                    name: file.clone(),
                    reads: nreads,
                    bases: nbases,
                    q20: format!("{:.2}", qual20 as f64 / nbases as f64 * 100.0),
                    q30: format!("{:.2}", qual30 as f64 / nbases as f64 * 100.0),
                    nx: n50 as u64,
                    gc: format!("{:.2}", gcbases as f64 / nbases as f64 * 100.0)
                });
                total_bases += nbases;
                total_reads += nreads;
            } 
        }
    };

    let upload_files = move |evt: FormEvent| async move {
        if let Some(file_engine) = evt.files() {
            busy.set(true);
            read_files(file_engine).await;
            busy.set(false);
        }
    };
    
//    let downloadtable = "<a> Download table </a>";

    rsx! {
        style { 
            {include_str!("../assets/custom.css")} 
        }
        div {
            id: "title",
            h2 { "Simple fastx file analysis" }
        }
        div {
            p{
            "This application runs basic analysis on sequencing files in fastq format. 
            Select or drop .fastq or .fastq.gz files to analyse" 
            }
        }
        div {
            label { r#for: "textreader", "" }
            input {
                class: "usercontrols",
                r#type: "file",
                accept: ".fastq,.gz",
                multiple: true,
                name: "textreader",
                //directory: enable_directory_upload,
                directory: false,
                onchange: upload_files
            }

            label {""}
            button {
                class: "usercontrols",
                onclick: move |_| {
                    files_uploaded.write().clear(); 
                    total_bases.set(0); 
                    total_reads.set(0)
                },
                "Clear" 
            }
            
            if files_uploaded.len() > 0{
                label {r#for: "numbers", "Number format: "}
                select {
                    r#name: "numbers", id: "numbers",
                    class: "usercontrols",
                    multiple: false,
                    oninput: move |ev| {
                        numbers.set(ev.value())
                    },
                    option {value: "none", "Plain numeric"},
                    option {value: "comma", "Thousands sep"},
                    option {value: "human", "SI suffix"}
                }
            }
        }
        // div {
        //     id: "drop-zone",
        //     prevent_default: "ondragover ondrop",
        //     background_color: if hovered() { "#3498DB" } else { "#D6EAF8" },
        //     ondragover: move |_| hovered.set(true),
        //     ondragleave: move |_| hovered.set(false),
        //     ondrop: move |evt| async move {
        //         hovered.set(false);
        //         if let Some(file_engine) = evt.files() {
        //             busy.set(true);
        //             read_files(file_engine).await;
        //             busy.set(false);
        //         }
        //     },
        //     "Drop files here"
        // }
        if busy() {
            div {
                // this loader should run on separate thread
                div {
                    class: "loader",
                }
                p {
                    "Please wait... {files_uploaded.len()} files processed"
                }
            }
        }
        if files_uploaded.len() > 0 {
            table {
                id: "resultstable",
                thead {
                tr {
                    th {"File"}
                    th {"Reads"}
                    th {"Bases"}
                    th {"N50"}
                    th {"GC%"}
                    th {"Q20%"}
                    th {"Q30%"}
                }
                }
                tbody {
                    {maketable(files_uploaded, numbers(), total_reads, total_bases)}
                }
            }
        }
    }
}