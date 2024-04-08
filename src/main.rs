//! This example shows how to use the `file` methods on FormEvent and DragEvent to handle file uploads and drops.
//!
//! Dioxus intercepts these events and provides a Rusty interface to the file data. Since we want this interface to
//! be crossplatform,
//! 
use std::sync::Arc;

use bio::io::fastq;

use dioxus::prelude::*;
use dioxus::{html::HasFileData, prelude::dioxus_elements::FileEngine};



fn main() {
    launch(app);
}

struct UploadedFile {
    name: String,
    reads: u64,
    bases: u64,
}

fn app() -> Element {
    let mut enable_directory_upload = use_signal(|| false);
    //let mut files_uploaded = use_signal(|| Vec::new() as Vec<UploadedFile>);
    let mut files_uploaded = use_signal(|| Vec::new() as Vec<UploadedFile>);
    let mut hovered = use_signal(|| false);


    let read_files = move |file_engine: Arc<dyn FileEngine>| async move {
        let files = file_engine.files();
        for file_name in &files {
            let mut nreads: u64 = 0;
            let mut nbases: u64 = 0;

            if let Some(contents) = file_engine.read_file_to_string(file_name).await {
                let mut recs = fastq::Reader::new(contents.as_bytes()).records();
                //let mut recs = fastq::Reader::new(decoder).records();
                //lines += contents.lines().count() as u64;
                while let Some(Ok(rec)) = recs.next() {
                    nreads += 1;
                    nbases += rec.seq().len() as u64;
                }
                //chrs += contents.chars().count() as u64;
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

        h1 { "File Upload Example" }
        //p { "Drop a .txt, .rs, or .js file here to read it" }
        button { onclick: move |_| files_uploaded.write().clear(), "Clear files" }

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
            label { r#for: "textreader", "Upload text/rust files and read them" }
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
            thead {
                tr {
                    th {"File"}
                    th {"Reads"}
                    th {"Bases"}
                }
            }
            tbody {
                for file in files_uploaded.read().iter() {
                    tr {
                        td {"{file.name}"}
                        td {"{file.reads}"}
                        td {"{file.bases}"}
                    }
                }
            }
        }
    }
}