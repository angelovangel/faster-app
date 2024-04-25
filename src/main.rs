
use std::sync::Arc;
use std::io;
use std::io::prelude::*;

use bio::io::fastq;
use flate2::read::MultiGzDecoder;

use dioxus::prelude::*;
use dioxus::prelude::dioxus_elements::FileEngine;
use dioxus::desktop::{Config, WindowBuilder};

use indicatif::HumanCount;

mod modules;


fn main() {
    LaunchBuilder::desktop()
        .with_cfg(Config::new().with_window(WindowBuilder::new().with_focused(true).with_title("faster-app"))
        )
        .launch(app)
}
struct UploadedFile {
    name: String,
    reads: u64,
    bases: u64,
    nx: u64,
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

fn app() -> Element {
    let mut reads_processed = use_signal(|| 0);
    //let mut enable_directory_upload = use_signal(|| false);
    let mut pretty_numbers = use_signal(|| true);
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
            let mut qual20: i64 = 0;
            let mut qual30: i64 = 0;
            let mut len_vector: Vec<i64> = Vec::new();

            if let Some(bytes) = file_engine.read_file(&file).await {
                let recs2 = decode_reader(bytes, file).await.unwrap();
                let mut recs = fastq::Reader::new(recs2.as_bytes()).records();
                    
                while let Some(Ok(rec)) = recs.next() {
                    nreads += 1;
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
                    nx: n50 as u64
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
    
    let downloadtable = "<a> Download table </a>";

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
            
        //}
        //div {
            if files_uploaded.len() > 0{
            label { r#for: "pretty-numbers", "Use thousands separator" }
            input {
                class: "usercontrols",
                r#type: "checkbox",
                checked: pretty_numbers,
                oninput: move |evt| pretty_numbers.set(evt.checked()),
            }
            label {""}
            button {
                class: "usercontrols",
                dangerous_inner_html: "{downloadtable}",
                //onclick: open_compose_window
            
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
                    th {"Q20%"}
                    th {"Q30%"}
                }
                }
                tbody {
                    if pretty_numbers() {
                        for file in files_uploaded.read().iter() {
                            tr {
                            td {"{file.name}"}
                            td {"{HumanCount(file.reads)}"}
                            td {"{HumanCount(file.bases)}"}
                            td {"{HumanCount(file.nx)}"}
                            td {"{file.q20}"}
                            td {"{file.q30}"}
                            }
                        }
                        tr {
                        td {"Total"}
                        td {"{HumanCount(total_reads())}"}
                        td {"{HumanCount(total_bases())}"}
                        }
                    } else {
                        for file in files_uploaded.read().iter() {
                            tr {
                            td {"{file.name}"}
                            td {"{file.reads}"}
                            td {"{file.bases}"}
                            td {"{file.nx}"}
                            td {"{file.q20}"}
                            td {"{file.q30}"}
                            }
                        }
                        tr {
                        td {"Total"}
                        td {"{total_reads}"}
                        td {"{total_bases}"}
                        }
                    }
                }
        }
        }
    }
}