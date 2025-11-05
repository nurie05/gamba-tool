## GAMBA: Gene Aggregation tool for Multicistronic Block Annotation

GAMBA is a tool coded in **_Rust_** that identifies polycistronic transcriptional units (operons) among the annotated transcripts of a GTF file.
The detection is based on transcript position, structure and expression level.
It requires a GTF files with the values of 'cov' and 'fpkm', used during operon selection.

### Usage

This tool has few parameters that can be specified (use --help to get the usage message):

    -f, --file <FILE>            Path to the input GTF file
    -t, --threshold <THRESHOLD>  Coverage threshold multiplier [default: 1]
    -o, --output <OUTPUT>        Output file prefix
        --log <LOG>              Log file path
    -h, --help                   Print help

The generic command line for the default usage is:

    gamba-tool --file FILE.gtf --threshold 1 --output FILE-PREFIX --log FILE-PREFIX_gamba.log

Notice that only the '--file' parameter is mandatory. In case the '--output' is not specified it will use the file name as output prefix for all files, including the log file.
