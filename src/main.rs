
use std::sync::Arc;
use std::io;
use std::io::prelude::*;

use bio::io::fastq;
use flate2::read::MultiGzDecoder;

use dioxus::prelude::*;
use dioxus::{html::HasFileData, prelude::dioxus_elements::FileEngine};

use indicatif::HumanCount;



fn main() {
    launch(app);
}

struct UploadedFile {
    name: String,
    reads: u64,
    bases: u64,
}

 fn decode_reader(bytes: Vec<u8>, filename: &String) -> io::Result<String> {
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
    let mut enable_directory_upload = use_signal(|| false);
    //let mut pretty_numbers = use_signal(|| true);
    let mut files_uploaded = use_signal(|| Vec::new() as Vec<UploadedFile>);
    let mut hovered = use_signal(|| false);


    let read_files = move |file_engine: Arc<dyn FileEngine>| async move {
        let files = file_engine.files();
        for file_name in &files {
            let mut nreads: u64 = 0;
            let mut nbases: u64 = 0;
            // decide which reader to use based on extension

            if let Some(bytes) = file_engine.read_file(&file_name).await {
                let recs2 = decode_reader(bytes, file_name).unwrap();
                let mut recs = fastq::Reader::new(recs2.as_bytes()).records();
                    
                while let Some(Ok(rec)) = recs.next() {
                    nreads += 1;
                    nbases += rec.seq().len() as u64;
                }

                files_uploaded.write().push(UploadedFile {
                    name: file_name.clone(),
                    //contents,
                    reads: nreads,
                    bases: nbases
                });
            } 
        }
    };

    let upload_files = move |evt: FormEvent| async move {
        if let Some(file_engine) = evt.files() {
            read_files(file_engine).await;
        }
    };

    rsx! {
        style { {include_str!("../assets/file_upload.css")} }
        div {
            id: "title",
            h2 { "Fastq analysis app in web assembly" }
        }
        p { 
        "This is a Wasm application that runs basic analysis on sequencing files in fastq format. 
        Drop fastq files to analyse. Analysis is done in the browser, no data is sent out." 
        }

        div {
            label { r#for: "directory-upload", "Enable directory upload" }
            input {
                r#type: "checkbox",
                id: "directory-upload",
                checked: enable_directory_upload,
                oninput: move |evt| enable_directory_upload.set(evt.checked()),
            }
        }

        div {
            label { r#for: "textreader", "Upload fastx files" }
            input {
                r#type: "file",
                accept: ".txt,.fastq,.gz",
                multiple: true,
                name: "textreader",
                directory: enable_directory_upload,
                onchange: upload_files,
            }
        }
        div {
            label {"Clear files"}
            button { onclick: move |_| files_uploaded.write().clear(), "Clear" }
        }
        
        div {
            id: "drop-zone",
            prevent_default: "ondragover ondrop",
            background_color: if hovered() { "lightblue" } else { "lightgray" },
            ondragover: move |_| hovered.set(true),
            ondragleave: move |_| hovered.set(false),
            ondrop: move |evt| async move {
                hovered.set(false);
                if let Some(file_engine) = evt.files() {
                    read_files(file_engine).await;
                }
            },
            "Drop files here"
        }
        

        table {
            id: "resultstable",
            thead {
                tr {
                    th {"File"}
                    th {"Reads"}
                    th {"Bases"}
                }
            }
            tbody {
                for file in files_uploaded.read().iter().rev() {
                    tr {
                        td {"{file.name}"}
                        td {"{HumanCount(file.reads)}"}
                        td {"{HumanCount(file.bases)}"}
                    }
                }
            }
        }
    }
}