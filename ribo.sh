#!/usr/bin/env bash
search_dir=data/Berlin_AKI_tubuloid_bulk_RNA_seq
arr=()
arr_out=()
for entry in "$search_dir"/*.fq.gz
do
    arr+=($entry)
    arr_out+=($entry'output.fq')
done
echo "$arr_out"
