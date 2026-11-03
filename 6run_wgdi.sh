#!/usr/bin/env bash

set -euo pipefail
cd "$(dirname "$0")"
export MPLBACKEND=Agg

[[ -s folder_list.txt ]] || { echo "错误：folder_list.txt 不存在或为空。"; exit 1; }
mapfile -t species_ids < folder_list.txt

check_outputs() {
    local folder="$1"
    local base_id="$2"
    local target_id="$3"
    local prefix
    local required=()

    if [[ "$base_id" == "$target_id" ]]; then
        prefix="$base_id"
        required=(
            "$prefix.png"
            "$prefix.icl"
            "$prefix.ks"
            "$prefix.bi.csv"
            "$prefix.bi.c.csv"
            "$prefix.bk.png"
            "$prefix.kp.png"
            "$prefix.pf.png"
        )
    else
        prefix="${base_id}_${target_id}"
        required=(
            "$prefix.png"
            "$prefix.icl"
            "$prefix.ks"
            "$prefix.bi.csv"
            "$prefix.bk.png"
            "$prefix.kp.png"
            "$prefix.pf.png"
        )
    fi

    local missing=()
    local name
    for name in "${required[@]}"; do
        [[ -f "$folder/$name" ]] || missing+=("$name")
    done
    if [[ ${#missing[@]} -gt 0 ]]; then
        printf '%s\n' "${missing[@]}" > "$folder/missing_outputs.txt"
        return 1
    fi
    rm -f "$folder/missing_outputs.txt"
    return 0
}

ran=0
skipped=0
failed=0

for base_id in "${species_ids[@]}"; do
    for target_id in "${species_ids[@]}"; do
        folder="Auto_wgdi/${base_id}_base/${base_id}_${target_id}"

        if [[ ! -f "$folder/run.sh" || ! -f "$folder/total.conf" ]]; then
            echo "错误：配置缺失：$folder"
            failed=$((failed + 1))
            continue
        fi

        if [[ -f "$folder/.wgdi.done" ]] && check_outputs "$folder" "$base_id" "$target_id"; then
            echo "WGDI 已完成，跳过：$base_id / $target_id"
            skipped=$((skipped + 1))
            continue
        fi

        rm -f "$folder/.wgdi.done"
        echo "运行 WGDI：$base_id / $target_id"

        if ! (
            cd "$folder"
            bash run.sh
        ) > "$folder/run.log" 2>&1; then
            echo "错误：WGDI 运行失败：$base_id / $target_id"
            tail -n 10 "$folder/run.log" || true
            failed=$((failed + 1))
            continue
        fi

        if ! check_outputs "$folder" "$base_id" "$target_id"; then
            echo "错误：WGDI 返回成功，但结果文件不完整：$base_id / $target_id"
            cat "$folder/missing_outputs.txt"
            failed=$((failed + 1))
            continue
        fi

        touch "$folder/.wgdi.done"
        ran=$((ran + 1))
        echo "WGDI 完成：$base_id / $target_id"
    done
done

echo "WGDI 汇总：运行 $ran 个，已跳过 $skipped 个，失败 $failed 个。"
[[ $failed -eq 0 ]]
