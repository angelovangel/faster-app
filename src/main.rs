use std::io::{Write, Cursor};
use std::fs::File;
use std::path::Path;
use std::sync::Arc;
//use std::io::Cursor;
use std::time::{Instant, Duration};
use arboard::Clipboard;
use rfd::FileDialog;

use bio::io::fastq;
use flate2::read::MultiGzDecoder;

use dioxus::desktop::{Config, LogicalSize, WindowBuilder};
//use dioxus::prelude::dioxus_elements::FileEngine;
use dioxus::prelude::*;

use human_repr::HumanCount as HC;
use indicatif::{HumanCount, HumanDuration};
use tokio::task;

mod modules;
mod components;

fn main() {
    let config = Config::new()
    .with_window(
        WindowBuilder::new()
            .with_focused(true)
            .with_title("fasterx")
            .with_inner_size(LogicalSize::new(1300, 700)),
    );
    //.with_custom_head(custom_head.to_string());

    LaunchBuilder::desktop()
        .with_cfg(config)
        .launch(app)
}
#[derive(Clone)]
struct UploadedFile {
    name: String,
    basename: String,
    reads: u64,
    bases: u64,
    nx: u64,
    l_vector: Vec<i64>,
    gc: String,
    q20: String,
    q30: String,
    m_qscore: u8, // median q score
    q_vector: Vec<u8>, // Add this field to store the quality scores
}

async fn decode_reader(bytes: Vec<u8>, filename: &String) -> std::io::Result<Box<dyn std::io::Read + Send>> {
    if filename.ends_with(".gz") {
        // Use a Cursor to wrap the bytes and pass it to MultiGzDecoder
        let gz = MultiGzDecoder::new(Cursor::new(bytes));
        Ok(Box::new(gz))
    } else {
        // Wrap the bytes in a Cursor for non-gzipped files
        Ok(Box::new(Cursor::new(bytes)))
    }
}

fn maketable(
    entries: Signal<Vec<UploadedFile>>,
    name_type: String,
    numbers_type: String,
    binsize: Signal<usize>,
    treads: Signal<u64>,
    tbases: Signal<u64>,
    sort_by: Signal<(String, bool)>, // Track column and sort direction
) -> Element {
    let mut sorted_entries = entries.read().clone();

    // Sort entries based on the current column and direction
    let (column, ascending) = sort_by.read().clone();
    sorted_entries.sort_by(|a, b| {
        let order = match column.as_str() {
            "name" => a.name.cmp(&b.name),
            "reads" => a.reads.cmp(&b.reads),
            "bases" => a.bases.cmp(&b.bases),
            "nx" => a.nx.cmp(&b.nx),
            "gc" => a.gc.cmp(&b.gc),
            "q20" => a.q20.cmp(&b.q20),
            "q30" => a.q30.cmp(&b.q30),
            "m_qscore" => a.m_qscore.cmp(&b.m_qscore),
            _ => std::cmp::Ordering::Equal,
        };
        if ascending {
            order
        } else {
            order.reverse()
        }
    });

    rsx! {
        for f in sorted_entries.iter() {
            tr {
                if name_type == "fullpath" {
                    td {"{f.name}"}
                } else {
                    td {"{f.basename}"}
                }
                if numbers_type == "comma" {
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
                td {
                    class: "histogram-cell",
                    dangerous_inner_html: "{generate_l_histogram(&f.l_vector, binsize())}"
                }
                td {"{f.gc}"}
                //td {"{f.q20}"}
                td {"{f.q30}"}
                td {"{f.m_qscore}"}
                td {
                    class: "histogram-cell",
                    dangerous_inner_html: "{generate_q_histogram(&f.q_vector)}" // Render the histogram as HTML
                }
            }
        }
        tr {
            class: "total-row",
            td {"Total ({entries.len()} files)"}
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
            td {  }
            td {  }
            td {  }
            td {  }
            td {  }
        }
    }
}

fn copy_to_clipboard(f_uploaded: Signal<Vec<UploadedFile>>) {
    let mut csv_data = String::new();
    csv_data.push_str("File,Reads,Bases,N50,GC%,Q20%,Q30%,Median_Qscore\n");

    for file in f_uploaded.read().iter() {
        csv_data.push_str(&format!(
            "{},{},{},{},{},{},{},{}\n",
            //file.name,
            file.basename,
            file.reads,
            file.bases,
            file.nx,
            file.gc,
            file.q20,
            file.q30,
            file.m_qscore
        ));
    }

    // Use the `arboard` crate to copy the data to the clipboard
    let mut clipboard = Clipboard::new().unwrap();
    clipboard.set_text(csv_data).unwrap();
}

fn save_html(f_uploaded: Signal<Vec<UploadedFile>>, numbers_type: String, name_type: String, binsize: usize) {
    let html_data = {
        let mut html_data = String::new();

        // Add basic HTML structure
        html_data.push_str("<!DOCTYPE html>\n<html>\n<head>\n");
        html_data.push_str("<title>FasterX Results</title>\n");
        html_data.push_str("<style>\n");
        html_data.push_str(include_str!("../assets/custom.css")); // Include your CSS styles
        html_data.push_str("</style>\n");
        //html_data.push_str("</head>\n<body>\n");
        //html_data.push_str("<h1>FasterX Results</h1>\n");
        html_data.push_str("<h2 id='title'>fasterX app results</h2>\n");
        html_data.push_str("<p>Generated on: ");
        html_data.push_str(&format!("{}", chrono::Local::now().format("%Y-%m-%d %H:%M:%S")));
        html_data.push_str("</p>\n");
        
        // Add table structure
        html_data.push_str("<table id='resultstable'>\n<thead>\n<tr>\n");
        html_data.push_str("<th>File</th><th>Reads</th><th>Bases</th><th>N50</th><th class='histogram-header'>Length Histogram</th><th>GC%</th><th>Q30%</th><th>Median Qscore</th><th class='histogram-header'>Qscore Histogram</th>\n");
        html_data.push_str("</tr>\n</thead>\n<tbody>\n");

        // Add table rows
        for file in f_uploaded.read().iter() {
            html_data.push_str("<tr>\n");
            if name_type == "fullpath" {
                html_data.push_str(&format!("<td>{}</td>\n", file.name));
            } else {
                html_data.push_str(&format!("<td>{}</td>\n", file.basename));
            }

            if numbers_type == "comma" {
                html_data.push_str(&format!("<td>{}</td>\n", HumanCount(file.reads)));
                html_data.push_str(&format!("<td>{}</td>\n", HumanCount(file.bases)));
                html_data.push_str(&format!("<td>{}</td>\n", HumanCount(file.nx)));
            } else if numbers_type == "human" {
                html_data.push_str(&format!("<td>{}</td>\n", file.reads.human_count_bare()));
                html_data.push_str(&format!("<td>{}</td>\n", file.bases.human_count_bare()));
                html_data.push_str(&format!("<td>{}</td>\n", file.nx.human_count_bare()));
            } else {
                html_data.push_str(&format!("<td>{}</td>\n", file.reads));
                html_data.push_str(&format!("<td>{}</td>\n", file.bases));
                html_data.push_str(&format!("<td>{}</td>\n", file.nx));
            }

            // Embed the length histogram as raw HTML
            html_data.push_str(&format!(
                "<td class='histogram-cell'>{}</td>\n",
                generate_l_histogram(&file.l_vector, binsize)
            ));

            html_data.push_str(&format!("<td>{}</td>\n", file.gc));
            html_data.push_str(&format!("<td>{}</td>\n", file.q30));
            html_data.push_str(&format!("<td>{}</td>\n", file.m_qscore));

            // Embed the Qscore histogram as raw HTML
            html_data.push_str(&format!(
                "<td class='histogram-cell'>{}</td>\n",
                generate_q_histogram(&file.q_vector)
            ));

            html_data.push_str("</tr>\n");
        }

        html_data.push_str("</tbody>\n</table>\n");
        html_data.push_str("</body>\n</html>");

        html_data
    };

    // Use `spawn_blocking` to avoid blocking the main thread
    task::spawn_blocking(move || {
        let path = std::env::current_dir().unwrap();

        if let Some(mypath) = FileDialog::new()
            .set_title("Save HTML File")
            .set_directory(&path)
            .set_file_name("fasterx_results.html")
            .save_file()
        {
            // Save the HTML to the selected file
            if let Ok(mut file) = File::create(mypath) {
                if let Err(e) = file.write_all(html_data.as_bytes()) {
                    eprintln!("Failed to write HTML file: {}", e);
                }
            } else {
                eprintln!("Failed to create HTML file");
            }
        } else {
            eprintln!("Save operation canceled");
        }
    });
}

fn format_thead(sortby: Signal<(String, bool)>, sortcol: &str) -> String {
    if sortby.read().0 == sortcol {
        if sortby.read().1 {
            "↑".to_string()
        } else {
            "↓".to_string()
        }
    } else {
        "↕".to_string()
    }
}

fn generate_q_histogram(q_vector: &[u8]) -> String {
    let mut bins = [0; 30]; // Create 30 bins for the histogram
    let max_bin_index = bins.len() - 1; // Index of the last bin

    for &q in q_vector {
        let bin = (q / 2) as usize; // bin spans 2 qvalues e.g. 8-10
        if bin < bins.len() {
            bins[bin] += 1;
        } else {
            bins[max_bin_index] += 1; // Increment the last bin for out-of-range values
        }
    }

    // Find the maximum count to normalize the bar heights
    let max_count = *bins.iter().max().unwrap_or(&1);
    let sum_count = bins.iter().sum::<u64>();

    // Generate HTML for the bar chart
    bins.iter()
        .enumerate()
        .map(|(i, &count)| {
            let height = (count as f64 / max_count as f64) * 100.0; // Normalize height to a percentage
            let percent = format!("{:.0}", count as f64 / sum_count as f64 * 100.0);
            format!(
                r#"<div class="bar" style="height: {height}%;">
                    <div class="tooltip">Q {range_start}-{range_end}: {count} reads ({percent} %)</div>
                </div>"#,
                height = height,
                range_start = i * 2,
                range_end = i * 2 + 2,
                count = count.human_count_bare()
            )
        })
        .collect::<Vec<_>>()
        .join("") // Combine all bars into a single string
}

fn generate_l_histogram(l_vector: &[i64], binsize: usize) -> String {
    let mut bins = [0; 30]; // Create 30 bins for the histogram
    let max_bin_index = bins.len() - 1; // Index of the last bin
    
    for &l in l_vector {
        let bin = l as usize/ binsize; // 1 bin spans binsize bp
        if bin < bins.len() {
            bins[bin] += 1;
        } else {
            bins[max_bin_index] += 1; // Increment the last bin for out-of-range values
        }
    }

    // Find the maximum count to normalize the bar heights
    let max_count = *bins.iter().max().unwrap_or(&1);
    let sum_count = bins.iter().sum::<u64>();

    // Generate HTML for the bar chart
    bins.iter()
        .enumerate()
        .map(|(i, &count)| {
            let height = (count as f64 / max_count as f64) * 100.0; // Normalize height to a percentage
            let percent = format!("{:.0}", count as f64 / sum_count as f64 * 100.0);
            let range_start = i * binsize;
            let range_end = if i == max_bin_index {
                format!(">{}", range_start.human_count_bare()) // Use ">range_start" for the last bin
            } else {
                format!("{}-{}", range_start.human_count_bare(), (range_start + binsize).human_count_bare())
            };

            format!(
                r#"<div class="bar" style="height: {height}%;">
                    <div class="tooltip">Length {range_end}: {count} reads ({percent} %) </div>
                </div>"#,
                height = height,
                range_end = range_end,
                count = count.human_count_bare()
            )
        })
        .collect::<Vec<_>>()
        .join("") // Combine all bars into a single string
}

fn app() -> Element {
    //let mut enable_directory_upload = use_signal(|| false);
    let mut numbers = use_signal(|| String::new());
    let mut basesperbin = use_signal(|| 1000);
    let mut name_type_sig = use_signal(|| String::new());
    let mut files_uploaded = use_signal(|| Vec::new() as Vec<UploadedFile>);
    let mut files_count = use_signal(|| 0);
    
    let mut total_reads = use_signal(|| 0);
    let mut total_bases = use_signal(|| 0);
    
    let mut busy = use_signal(|| false);
    let mut ready = use_signal(|| false);
    let mut cancel_processing = use_signal(|| false);
    let mut start_time = use_signal(|| Instant::now());
    let mut myduration = use_signal(|| String::new());
    let mut progress_percentage = use_signal(|| 0.0);
    let mut show_popup = use_signal(|| false);
    let mut sort_by = use_signal(|| ("name".to_string(), true)); // Default sort by name ascending

    let read_files = move |file_engine: Arc<dyn dioxus_elements::FileEngine>| async move {
        let files = file_engine.files();
        for file in &files {
            if *cancel_processing.read() {
                break; // Exit the loop if processing is canceled
            }

            let mut nreads: u64 = 0;
            let mut nbases: u64 = 0;
            let mut gcbases: u64 = 0;
            let mut qual20: i64 = 0;
            let mut qual30: i64 = 0;
            let mut len_vector: Vec<i64> = Vec::new();
            let mut qual_vector: Vec<u8> = Vec::new();

            if let Some(bytes) = file_engine.read_file(&file).await {
                let filepath = Path::new(&file);
                let basename = filepath.file_name().unwrap().to_str().unwrap();
                let reader = decode_reader(bytes, file).await.unwrap();
                let mut recs = fastq::Reader::new(reader).records();

                while let Some(Ok(rec)) = recs.next() {
                    if *cancel_processing.read() {
                        break; // Exit the loop if processing is canceled
                    }
                    nreads += 1;
                    gcbases += modules::get_gc_bases(rec.seq());
                    nbases += rec.seq().len() as u64;
                    qual20 += modules::get_qual_bases(rec.qual(), 53); // 33 offset
                    qual30 += modules::get_qual_bases(rec.qual(), 63); // 33 offset
                    len_vector.push(rec.seq().len() as i64);
                    qual_vector.push(modules::qscore_mean(rec.qual()));

                    tokio::task::yield_now().await; // Yield to allow the UI to update
                }
                let n50 = modules::get_nx(&mut len_vector, 0.5);
                let median_qscore = modules::median(&mut qual_vector);

                files_uploaded.write().push(UploadedFile {
                    name: file.clone(),
                    basename: basename.to_string(),
                    reads: nreads,
                    bases: nbases,
                    q20: format!("{:.2}", qual20 as f64 / nbases as f64 * 100.0),
                    q30: format!("{:.2}", qual30 as f64 / nbases as f64 * 100.0),
                    nx: n50 as u64,
                    gc: format!("{:.2}", gcbases as f64 / nbases as f64 * 100.0),
                    m_qscore: median_qscore,
                    q_vector: qual_vector.clone(),
                    l_vector: len_vector.clone(),
                });

                progress_percentage.set((files_uploaded.len() as f64 / *files_count.read() as f64) * 100.0);
                tokio::task::yield_now().await; // Yield to allow the UI to update
                tokio::time::sleep(std::time::Duration::from_millis(500)).await; // this makes sure that the progress is drawn
                total_bases += nbases;
                total_reads += nreads;
            }
        }

        if *cancel_processing.read() {
            busy.set(false); // Reset the busy state
            progress_percentage.set(0.0); // Reset progress
        }
    };


    let upload_files = move |evt: FormEvent| {
        if let Some(file_engine) = evt.files() {
            busy.set(true);
            files_count.set(file_engine.files().len());
            start_time.set(Instant::now()); // Record the start time
            
            spawn(async move { 
                read_files(file_engine).await; // Process all files
                myduration.set(HumanDuration(Instant::now() - start_time()).to_string()); 

                // Show the "ready" popup once
                ready.set(true);
                busy.set(false); 
                progress_percentage.set(0.0);
                tokio::time::sleep(std::time::Duration::from_secs(3)).await; // Wait for 3 seconds
                ready.set(false); // Hide the popup
            });
        } 
    };
    
    rsx! {
        style { 
            {include_str!("../assets/custom.css")} 
        }
        
        components::app_title {}
        
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

            if files_uploaded.len() > 0{
                button {
                    class: "usercontrols",
                    onclick: move |_| {
                        files_uploaded.write().clear();
                        total_bases.set(0);
                        total_reads.set(0);
                        files_count.set(0);
                        name_type_sig.set("basename".to_string());
                        progress_percentage.set(0.0);
                    },
                    "Clear"
                }
                button {
                    class: "usercontrols",
                    onclick: move |_| {
                        copy_to_clipboard(files_uploaded.clone());
                        show_popup.set(true);

                        // Use an async task to handle the delay, dioxus::prelude::spawn()
                        spawn(async move {
                            tokio::time::sleep(Duration::from_secs(3)).await;
                            show_popup.set(false);
                        });
                    },
                    "Copy to clipboard"
                }

                button {
                    class: "usercontrols",
                    onclick: move |_| {
                        save_html(files_uploaded.clone(), numbers(), name_type_sig(), basesperbin());
                    },
                    "Save as HTML"
                }
                
                select {
                    r#name: "name_type_sig",
                    class: "usercontrols",
                    multiple: "false",
                    oninput: move |ev| {
                        name_type_sig.set(ev.value())
                    },
                    option {value: "basename", "Base filename"}
                    option {value: "fullpath", "Full path"},
                }

                //label {r#for: "numbers", ""}
                div {
                    class: "tooltip-container",
                    div {
                        class: "tooltip",
                        "Choose number format", 
                    }
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
                
                div {
                    class: "tooltip-container",
                    div {
                        class: "tooltip",
                        "Length histogram: bases per bin", 
                        br {}, 
                        "Showing {basesperbin().human_count_bare()}-{(basesperbin * 30).human_count_bare()}", " bp",
                    }
                    input {
                        r#type: "number",
                        id: "basesperbin",
                        //title: "Bases per bin for the length histogram",
                        class: "usercontrols",
                        value: "{basesperbin}", // Bind the current value of basesperbin
                        min: "50",
                        max: "5000",
                        step: "50", // Set the step size for increment/decrement
                        oninput: move |ev| {
                            if let Ok(value) = ev.value().parse::<usize>() {
                                if value > 0 {
                                    basesperbin.set(value); // Update the basesperbin signal
                                }
                            }
                        }
                    }
                }
            }
        }

        if files_uploaded.len() > 0 {
            table {
                id: "resultstable",
                thead {
                    tr {
                        th {
                            class: "sortable-header",
                            onclick: {
                                let current_sort = sort_by.read().1;
                                move |_| sort_by.set(("name".to_string(), !current_sort))
                            },
                            "File ",
                            {format_thead(sort_by, "name")}
                            // if sort_by.read().0 == "name" {
                            //     if sort_by.read().1 {
                            //         "↑" // Ascending
                            //     } else {
                            //         "↓" // Descending
                            //     }
                            // } else {
                            //     "↕" // Default indicator for unsorted columns
                            // }
                        }
                        th {
                            class: "sortable-header",
                            onclick: {
                                let current_sort = sort_by.read().1;
                                move |_| sort_by.set(("reads".to_string(), !current_sort))
                            },
                            "Reads ",
                            {format_thead(sort_by, "reads")}
                        }
                        th {
                            class: "sortable-header",
                            onclick: {
                                let current_sort = sort_by.read().1;
                                move |_| sort_by.set(("bases".to_string(), !current_sort))
                            },
                            "Bases ",
                            {format_thead(sort_by, "bases")}
                        }
                        th {
                            class: "sortable-header",
                            onclick: {
                                let current_sort = sort_by.read().1;
                                move |_| sort_by.set(("nx".to_string(), !current_sort))
                            },
                            "N50 ",
                            {format_thead(sort_by, "nx")}
                        }
                        th {
                            class: "histogram-header",
                            "Length histogram"
                        }
                        th {
                            class: "sortable-header",
                            onclick: {
                                let current_sort = sort_by.read().1;
                                move |_| sort_by.set(("gc".to_string(), !current_sort))
                            },
                            "GC% ",
                            {format_thead(sort_by, "gc")}
                        }
                        // th {
                        //     class: "sortable-header",
                        //     onclick: {
                        //         let current_sort = sort_by.read().1;
                        //         move |_| sort_by.set(("q20".to_string(), !current_sort))
                        //     },
                        //     "Q20% ",
                        //     {format_thead(sort_by, "q20")}
                        // }
                        th {
                            class: "sortable-header",
                            onclick: {
                                let current_sort = sort_by.read().1;
                                move |_| sort_by.set(("q30".to_string(), !current_sort))
                            },
                            "Q30% ",
                            {format_thead(sort_by, "q30")}
                        }
                        th {
                            class: "sortable-header",
                            onclick: {
                                let current_sort = sort_by.read().1;
                                move |_| sort_by.set(("m_qscore".to_string(), !current_sort))
                            },
                            "Median Qscore ",
                            {format_thead(sort_by, "m_qscore")}
                        }
                        th {
                            class: "histogram-header",
                            "Qscore histogram"
                        }
                    }
                }
                tbody {
                    {maketable(files_uploaded, name_type_sig(), numbers(), basesperbin, total_reads, total_bases, sort_by)}
                }
            }
        }

        // Footer with app version info
        footer {
            class: "app-footer",
            "fasterX app - v0.2.4 © 2025"
        }

        if *busy.read() {
            button {
                class: "usercontrols usercontrols-cancel",
                onclick: move |_| {
                    cancel_processing.set(true);

                },
                "Cancel processing ✕"
            }
            div {
                class: "popup",
                "Please wait... {files_uploaded.len()} of {files_count} files processed",
                div {
                    class: "progress-bar-container",
                    div {
                        class: "progress-bar",
                        style: "width: {progress_percentage()}%;",
                    }
                }
            }
        }

        if *ready.read() {
            div {
                class: "popup",
                "Finished processing {files_count} files in {myduration}!"
            }
        }

        if *show_popup.read() {
            div {
                class: "popup",
                "{files_uploaded.len()} entries copied to clipboard!"
            }
        }
    }
}
