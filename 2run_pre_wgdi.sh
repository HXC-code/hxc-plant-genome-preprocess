#!/usr/bin/env bash

set -euo pipefail
cd "$(dirname "$0")"

[[ -s folder_list.txt ]] || { echo "错误：folder_list.txt 不存在或为空。"; exit 1; }
[[ -s species_map.tsv ]] || { echo "错误：species_map.tsv 不存在或为空。"; exit 1; }

declare -A source_dir=()
declare -A pep_name=()
declare -A cds_name=()
declare -A gff_name=()
declare -A lens_name=()

while IFS=$'\t' read -r species_id directory pep cds gff lens; do
    source_dir[$species_id]="$directory"
    pep_name[$species_id]="$pep"
    cds_name[$species_id]="$cds"
    gff_name[$species_id]="$gff"
    lens_name[$species_id]="$lens"
done < species_map.tsv

mapfile -t species_ids < folder_list.txt
mkdir -p Auto_wgdi

copy_species_files() {
    local species_id="$1"
    local destination="$2"

    cp -f "${source_dir[$species_id]}/${pep_name[$species_id]}" "$destination/"
    cp -f "${source_dir[$species_id]}/${cds_name[$species_id]}" "$destination/"
    cp -f "${source_dir[$species_id]}/${gff_name[$species_id]}" "$destination/"
    cp -f "${source_dir[$species_id]}/${lens_name[$species_id]}" "$destination/"
}

comparison_count=0
for base_id in "${species_ids[@]}"; do
    base_folder="Auto_wgdi/${base_id}_base"
    mkdir -p "$base_folder"

    for other_id in "${species_ids[@]}"; do
        sub_folder="$base_folder/${base_id}_${other_id}"
        mkdir -p "$sub_folder"

        copy_species_files "$base_id" "$sub_folder"
        if [[ "$base_id" != "$other_id" ]]; then
            copy_species_files "$other_id" "$sub_folder"
        fi

        echo "已准备比较目录：$sub_folder"
        comparison_count=$((comparison_count + 1))
    done
done

echo "WGDI 输入准备完成：$comparison_count 个比较目录。"
