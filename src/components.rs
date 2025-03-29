use dioxus::prelude::*;

pub fn app_title() -> Element {
    rsx! {
        div {
            id: "title",
            h2 { "Simple fastx file analysis" }
        }
        div {
            p{
            "This application runs basic analysis on sequencing files in fastq format.
            Select fastq or fastq.gz files to analyse."
            }
        }
    }
}