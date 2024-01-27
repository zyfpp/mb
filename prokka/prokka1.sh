#!/bin/bash

genome_file=$(basename "$1")

input_dir="/home/zczl463/Scratch/mb/prok"
full_path_to_genome="${input_dir}/${genome_file}"

output_dir="/home/zczl463/Scratch/mb/presults/${genome_file}_prokka"

prokka --outdir "$output_dir" --prefix "${genome_file}" "$full_path_to_genome"
