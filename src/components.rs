use dioxus::prelude::*;

pub fn app_title(filescount: Signal<usize>) -> Element {
    rsx! {
        div {
            id: "title",
            if *filescount.read() == 0 {
                h3 { "Simple fastx file analysis" }
            } else {
                h5 { "Simple fastx file analysis" }
            }
            
        }
        div {
            if *filescount.read() == 0 {
                p{
                "This application runs basic analysis on sequencing files in fastq format.
                Select fastq or fastq.gz files to analyse (no data leaves the browser)."
                }
            }
        }
    }
}