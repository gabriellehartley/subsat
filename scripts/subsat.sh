#!/bin/bash

echo "launching subsat. Last Updated JAN 21 2026. Requires bedtools, TRF, trftobed.py, and difflib.py in path. contact: gabrielle.hartley@uconn.edu"

module load bedtools
module load TRF

export TMPDIR=./

# Exit on error
set -e

# Function to check if a command exists
require_cmd() {
    if ! command -v "$1" &> /dev/null; then
        echo "Error: Required command '$1' not found in PATH." >&2
        exit 1
    fi
}

# Check required commands
require_cmd complementBed
require_cmd bedtools
require_cmd trf

# Input validation
if [ $# -ne 4 ]; then
    echo "Usage: $0 <rmout_file> <genesgtf_file> <genomesizes> <genomefasta>" >&2
    exit 1
fi

rmout_file="$1"
genesgtf_file="$2"
genomesizes="$3"
genomefasta="$4"

# Check input files exist
for file in "$rmout_file" "$genomesizes" "$genomefasta"; do
    if [ ! -f "$file" ]; then
        echo "Error: File '$file' not found." >&2
        exit 1
    fi
done

# Output filenames
filename=$(basename "$genomefasta")
basename="${filename%%.*}"
rmbed_file="${basename}.rmbed"
rmbedsort_file="${basename}.rmsort.bed"
rmgapbed_file="${basename}.rmgaps.bed"
exon_file="${basename}.exons.bed"
exonsort_file="${basename}.exons.sort.bed"
gapbed_file="${basename}.gaps.bed"
gapbedfiltered_file="${basename}.gap.filtered.bed"
gapfasta_file="${basename}.gap.fasta"
trf_results="${basename}.trf.results"
trf_bed_temp="${basename}.trf.temp.bed"
trf_bed_split="${basename}.trf.split.bed"
trf_bed_final="${basename}.trf.final.bed"

echo "Generating $rmbed_file..."
tail -n +4 "$rmout_file" | awk '{
    print $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $15, $1, $2, $3
}' OFS='\t' > "$rmbed_file"

echo "sorting $rmbed_file..."
bedtools sort -i "$rmbed_file" -g "$genomesizes" > "$rmbedsort_file"

echo "Running complementBed -> $rmgapbed_file..."
complementBed -i "$rmbedsort_file" -g "$genomesizes" > "$rmgapbed_file"

echo "Generating $exon_file..."
awk '$3=="exon" {print $1, $4, $5, $3}' OFS='\t' "$genesgtf_file" > "$exon_file"

echo "Sorting $exon_file... "
bedtools sort -i "$exon_file" -g "$genomesizes" > "$exonsort_file"

echo "Generating gap file based on $rmgapbed_file and $exonsort_file..."
bedtools subtract -a "$rmgapbed_file" -b "$exonsort_file" > "$gapbed_file"

echo "Filtering regions >20000 bp into $gapbedfiltered_file..."
awk '{
    len = $3 - $2;
    if (len > 20000)
        print $1, $2, $3, len
}' OFS='\t' "$gapbed_file" > "$gapbedfiltered_file"

echo "Running bedtools getfasta on $gapbedfiltered_file -> $gapfasta_file..."
bedtools getfasta -fi "$genomefasta" -bed "$gapbedfiltered_file" > "$gapfasta_file"

echo "Running trf on $gapfasta_file -> $trf_results..."
trf "$gapfasta_file" 2 7 7 80 10 50 2000 -f -d -m > "$trf_results" &

# Wait for TRF to finish
wait

echo "Running trf to bed ..."

dat_file="${gapfasta_file}.2.7.7.80.10.50.2000.dat"

if [[ ! -f "$dat_file" ]]; then
    echo "ERROR: TRF .dat file not found: $dat_file" >&2
    exit 1
fi

module load python

trftobed.py --dat "$dat_file" --bed "$trf_bed_temp"

echo "Reformating bed file intervals..."

awk 'BEGIN {OFS="\t"} {sub(/:/, "\t", $1); print}' "$trf_bed_temp" | awk 'BEGIN {OFS="\t"} {sub(/-/, "\t", $2); print}' > "$trf_bed_split"

awk '{
    # fix intervals
    print $1, $2 + $4, $3 + $4, $6
}'  "$trf_bed_split" > "$trf_bed_final"

mkdir  trf_results
mv *html trf_results
mv *dat trf_results
mv *mask trf_results

echo "Completed."



