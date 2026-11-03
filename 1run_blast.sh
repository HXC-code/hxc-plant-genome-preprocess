#!/usr/bin/env bash

set -euo pipefail
cd "$(dirname "$0")"
shopt -s nullglob
project_root="$(pwd -P)"

blastp_path="D:/ProgramFiles/bioinfo/bin/blastp.exe"
makeblastdb_path="D:/ProgramFiles/bioinfo/bin/makeblastdb.exe"
blast_threads=16

[[ -x "$blastp_path" ]] || { echo "错误：找不到 blastp：$blastp_path"; exit 1; }
[[ -x "$makeblastdb_path" ]] || { echo "错误：找不到 makeblastdb：$makeblastdb_path"; exit 1; }

mkdir -p "Auto_pep/db"

map_tmp="$(mktemp)"
skipped_file="skipped_species.txt"
: > "$skipped_file"

declare -A used_ids=()
complete_count=0

for species_dir in */; do
    species_dir="${species_dir%/}"

    case "$species_dir" in
        Auto_*|_pipeline_*|.*) continue ;;
    esac

    pep_files=("$species_dir"/*.pep)
    if [[ ${#pep_files[@]} -ne 1 ]]; then
        echo "$species_dir：需要且只能有一个 .pep 文件，实际 ${#pep_files[@]} 个" >> "$skipped_file"
        continue
    fi

    pep_name="$(basename "${pep_files[0]}")"
    species_id="${pep_name%.pep}"

    if [[ ! "$species_id" =~ ^[A-Za-z0-9._-]+$ ]]; then
        echo "$species_dir：物种代码“$species_id”含有空格或不安全字符" >> "$skipped_file"
        continue
    fi
    if [[ -n "${used_ids[$species_id]:-}" ]]; then
        echo "$species_dir：物种代码 $species_id 与 ${used_ids[$species_id]} 重复" >> "$skipped_file"
        continue
    fi

    cds_name="${species_id}.cds"
    gff_name="${species_id}.new.gff"
    lens_name="${species_id}.lens"
    missing=()
    for required in "$pep_name" "$cds_name" "$gff_name" "$lens_name"; do
        [[ -s "$species_dir/$required" ]] || missing+=("$required")
    done
    if [[ ${#missing[@]} -gt 0 ]]; then
        echo "$species_dir：缺少 ${missing[*]}" >> "$skipped_file"
        continue
    fi

    used_ids[$species_id]="$species_dir"
    printf '%s\t%s\t%s\t%s\t%s\t%s\n' \
        "$species_id" "$species_dir" "$pep_name" "$cds_name" "$gff_name" "$lens_name" >> "$map_tmp"
    complete_count=$((complete_count + 1))
done

if [[ $complete_count -eq 0 ]]; then
    echo "错误：没有任何物种同时具备 .pep、.cds、.new.gff 和 .lens。"
    echo "详见：$skipped_file"
    rm -f "$map_tmp"
    exit 1
fi

LC_ALL=C sort -t $'\t' -k1,1 "$map_tmp" > species_map.tsv
cut -f1 species_map.tsv > folder_list.txt
rm -f "$map_tmp"

declare -A source_dir_by_id=()
declare -A pep_name_by_id=()

while IFS=$'\t' read -r species_id species_dir pep_name cds_name gff_name lens_name; do
    source_dir_by_id[$species_id]="$species_dir"
    pep_name_by_id[$species_id]="$pep_name"
done < species_map.tsv

mapfile -t species_ids < folder_list.txt

echo "准备执行 ${#species_ids[@]} 个完整物种的 BLAST。"
echo "不完整物种已写入：$skipped_file"

for species_id in "${species_ids[@]}"; do
    species_dir="${source_dir_by_id[$species_id]}"
    pep_name="${pep_name_by_id[$species_id]}"
    pep_path="$species_dir/$pep_name"
    db_prefix="Auto_pep/db/$species_id"
    db_prefix_abs="$project_root/Auto_pep/db/$species_id"

    echo "建立数据库：$species_id"
    if ! (
        cd "$species_dir"
        "$makeblastdb_path" \
            -in "$pep_name" \
            -dbtype prot \
            -out "$db_prefix_abs" \
            -logfile "$project_root/Auto_pep/db/${species_id}.makeblastdb.log"
    ); then
        echo "错误：$species_id 建库失败。"
        exit 1
    fi
done

pair_count=0

for base_id in "${species_ids[@]}"; do
    base_dir="${source_dir_by_id[$base_id]}"
    base_pep="$base_dir/${pep_name_by_id[$base_id]}"
    mkdir -p "Auto_pep/${base_id}_pep"

    for target_id in "${species_ids[@]}"; do
        target_dir="${source_dir_by_id[$target_id]}"
        target_pep="$target_dir/${pep_name_by_id[$target_id]}"
        db_prefix="Auto_pep/db/$target_id"
        output="Auto_pep/${base_id}_pep/${base_id}_${target_id}.blast"
        partial="${output}.part"
        output_abs="$project_root/$output"
        partial_abs="$project_root/$partial"
        db_prefix_abs="$project_root/Auto_pep/db/$target_id"

        echo "BLAST：$base_id -> $target_id"
        rm -f "$partial"
        if ! (
            cd "$base_dir"
            "$blastp_path" \
                -query "${pep_name_by_id[$base_id]}" \
                -db "$db_prefix_abs" \
                -out "$partial_abs" \
                -outfmt 6 \
                -evalue 1e-5 \
                -max_target_seqs 20 \
                -num_threads "$blast_threads"
        ); then
            echo "错误：BLAST 失败：$base_id -> $target_id"
            rm -f "$partial"
            exit 1
        fi
        mv -f "$partial" "$output"
        pair_count=$((pair_count + 1))
    done
done

echo "BLAST 完成：本轮重新生成 $pair_count 个。"
echo "完整物种数：${#species_ids[@]}；物种映射：species_map.tsv"
