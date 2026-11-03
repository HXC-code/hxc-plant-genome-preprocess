#!/usr/bin/env bash

set -euo pipefail
cd "$(dirname "$0")"

[[ -s folder_list.txt ]] || { echo "错误：folder_list.txt 不存在或为空。"; exit 1; }
mapfile -t species_ids < folder_list.txt

moved=0
already_there=0
missing=0

for base_id in "${species_ids[@]}"; do
    for target_id in "${species_ids[@]}"; do
        source_blast="Auto_pep/${base_id}_pep/${base_id}_${target_id}.blast"
        target_folder="Auto_wgdi/${base_id}_base/${base_id}_${target_id}"
        target_blast="$target_folder/${base_id}_${target_id}.blast"

        if [[ ! -d "$target_folder" ]]; then
            echo "错误：缺少目标文件夹：$target_folder"
            missing=$((missing + 1))
            continue
        fi
        if [[ -f "$source_blast" ]]; then
            mv -f "$source_blast" "$target_blast"
            echo "已写入最新 BLAST：$base_id -> $target_id"
            moved=$((moved + 1))
            continue
        fi
        if [[ -f "$target_blast" ]]; then
            echo "源文件已归位，跳过：$base_id -> $target_id"
            already_there=$((already_there + 1))
            continue
        else
            echo "错误：缺少 BLAST 文件：$source_blast"
            missing=$((missing + 1))
            continue
        fi
    done
done

echo "BLAST 归位完成：移动 $moved 个，已存在 $already_there 个，缺失 $missing 个。"
[[ $missing -eq 0 ]]
