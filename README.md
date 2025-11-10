## GAMBA: Gene Aggregation tool for Multicistronic Block Annotation

GAMBA is a tool coded in **_Rust_** that identifies polycistronic transcriptional units (operons) among the annotated transcripts of a GTF file.
The detection is based on transcript position, structure and expression level.
It requires a GTF files with the values of 'cov' and 'fpkm', used during operon selection.

### Installation

Source and binary files can be directly downloaded from the [**Releases**](https://github.com/nurie05/gamba-tool/releases) page on this repository.
Other option is to clone the repository and build the binary with **cargo**:

    git clone https://github.com/nurie05/gamba-tool.git
    cd gamba-tool
    cargo build --release
    target/release/gamba --help

Additionally it is possible to install it by using _Homebrew_:

    brew tap nurie05/homebrew-tools
    brew install nurie05/homebrew-tools/gamba

### Usage

This tool has few parameters that can be specified (use --help to get the usage message):

    -f, --file <FILE>                Path to the input GTF file
    -t, --threshold <THRESHOLD>      Coverage threshold multiplier [default: 1]
    -m, --min-overlap <MIN_OVERLAP>  Minimum percentage of exonic overlap to be considered 'contained transcript' [default: 0.5]
    -b, --bp-overlap <BP_OVERLAP>    Minimum bp overlap to consider exonic overlap [default: 50]
    -p, --prefix <PREFIX>            Output file prefix
    -o, --outdir <OUTDIR>            Output directory
        --log <LOG>                  Log file path
    -h, --help                       Print help
    -V, --version                    Print version

The generic command line for the default usage is:

    gamba --file FILE.gtf --threshold 1 --output OUTDIR -p FILE-PREFIX --log FILE-PREFIX_gamba.log

Notice that only the '--file' parameter is mandatory.  In case '--output' is not specified, it will output all the files in the same directory where the GTF file is.
In case '--prefix' is not specified, it will use the file name as output prefix for all files, including the log file.
