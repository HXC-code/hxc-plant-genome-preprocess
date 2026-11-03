#!/usr/bin/env bash

set -euo pipefail
cd "$(dirname "$0")"

[[ -s folder_list.txt ]] || { echo "错误：folder_list.txt 不存在或为空。"; exit 1; }
mapfile -t species_ids < folder_list.txt

configured=0
for species1 in "${species_ids[@]}"; do
    base_folder="Auto_wgdi/${species1}_base"

    for species2 in "${species_ids[@]}"; do
        [[ "$species1" == "$species2" ]] && continue

        sub_folder="$base_folder/${species1}_${species2}"
        [[ -d "$sub_folder" ]] || { echo "错误：缺少目录 $sub_folder"; exit 1; }

        cat > "$sub_folder/run.sh" <<EOF
#!/usr/bin/env bash
set -euo pipefail
export MPLBACKEND=Agg

mkdir -p .wgdi_steps
cat "${species1}.cds" "${species2}.cds" > all.cds
cat "${species1}.pep" "${species2}.pep" > all.pep

run_stage() {
    local stage_name="\$1"
    local output_file="\$2"
    shift 2
    local marker=".wgdi_steps/\${stage_name}.done"
    if [[ -f "\$marker" && -f "\$output_file" ]]; then
        echo "跳过已完成阶段：\$stage_name"
        return 0
    fi
    "\$@"
    [[ -f "\$output_file" ]] || {
        echo "错误：阶段 \$stage_name 未生成 \$output_file"
        return 1
    }
    touch "\$marker"
}

run_stage dotplot "${species1}_${species2}.png" wgdi -d total.conf
run_stage collinearity "${species1}_${species2}.icl" wgdi -icl total.conf
run_stage ks "${species1}_${species2}.ks" wgdi -ks total.conf
run_stage blockinfo "${species1}_${species2}.bi.csv" wgdi -bi total.conf
run_stage blockks "${species1}_${species2}.bk.png" wgdi -bk total.conf
run_stage kspeaks "${species1}_${species2}.kp.png" wgdi -kp total.conf
run_stage peaksfit "${species1}_${species2}.pf.png" wgdi -pf total.conf
EOF

        cat > "$sub_folder/total.conf" <<EOF
[dotplot]
blast = ${species1}_${species2}.blast
gff1 = ${species1}.new.gff
gff2 = ${species2}.new.gff
lens1 = ${species1}.lens
lens2 = ${species2}.lens
genome1_name = ${species1}
genome2_name = ${species2}
multiple = 1
score = 100
evalue = 1e-5
repeat_number = 10
position = order
blast_reverse = false
ancestor_left = none
ancestor_top = none
markersize = 0.5
figsize = 10,10
savefig = ${species1}_${species2}.png

[collinearity]
blast = ${species1}_${species2}.blast
gff1 = ${species1}.new.gff
gff2 = ${species2}.new.gff
lens1 = ${species1}.lens
lens2 = ${species2}.lens
blast_reverse = false
multiple = 1
process = 8
evalue = 1e-5
score = 100
grading = 50,25,10
mg = 25,25
pvalue = 1
repeat_number = 10
position = order
savefile = ${species1}_${species2}.icl

[ks]
cds_file = all.cds
pep_file = all.pep
align_software = muscle
pairs_file = ${species1}_${species2}.icl
ks_file = ${species1}_${species2}.ks

[blockinfo]
blast = ${species1}_${species2}.blast
gff1 = ${species1}.new.gff
gff2 = ${species2}.new.gff
lens1 = ${species1}.lens
lens2 = ${species2}.lens
collinearity = ${species1}_${species2}.icl
score = 100
evalue = 1e-5
repeat_number = 20
position = order
ks = ${species1}_${species2}.ks
ks_col = ks_NG86
savefile = ${species1}_${species2}.bi.csv

[blockks]
lens1 = ${species1}.lens
lens2 = ${species2}.lens
genome1_name = ${species1}
genome2_name = ${species2}
blockinfo = ${species1}_${species2}.bi.csv
pvalue = 0.2
tandem = true
tandem_length = 200
markersize = 1
area = 0,2
block_length = 5
figsize = 8,8
savefig = ${species1}_${species2}.bk.png

[kspeaks]
blockinfo = ${species1}_${species2}.bi.csv
pvalue = 0.2
tandem = true
block_length = 5
ks_area = 0,10
multiple = 1
homo = 0,1
fontsize = 9
area = 0,3
figsize = 10,6.18
savefig = ${species1}_${species2}.kp.png
savefile = ${species1}_${species2}.kp

[peaksfit]
blockinfo = ${species1}_${species2}.bi.csv
mode = median
bins_number = 200
ks_area = 0,10
fontsize = 9
area = 0,3
figsize = 10,6.18
shadow = true
savefig = ${species1}_${species2}.pf.png

[ksfigure]
ksfit = ks_fit_result.csv
labelfontsize = 15
legendfontsize = 15
xlabel = none
ylabel = none
title = none
area = 0,2
figsize = 10,6.18
shadow = false
savefig = ${species1}_${species2}.kf.png
EOF

        find "$sub_folder/.wgdi_steps" -maxdepth 1 -type f -name '*.done' -delete 2>/dev/null || true
        rm -f "$sub_folder/.wgdi.done"
        configured=$((configured + 1))
        echo "已生成不同物种配置：$species1 / $species2"
    done
done

echo "不同物种 WGDI 配置完成：$configured 个。"
