#!/usr/bin/env bash

set -euo pipefail
cd "$(dirname "$0")"

[[ -s folder_list.txt ]] || { echo "错误：folder_list.txt 不存在或为空。"; exit 1; }
mapfile -t species_ids < folder_list.txt

configured=0
for species in "${species_ids[@]}"; do
    sub_folder="Auto_wgdi/${species}_base/${species}_${species}"
    [[ -d "$sub_folder" ]] || { echo "错误：缺少目录 $sub_folder"; exit 1; }

    cat > "$sub_folder/run.sh" <<EOF
#!/usr/bin/env bash
set -euo pipefail
export MPLBACKEND=Agg

mkdir -p .wgdi_steps

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

run_stage dotplot "${species}.png" wgdi -d total.conf
run_stage collinearity "${species}.icl" wgdi -icl total.conf
run_stage ks "${species}.ks" wgdi -ks total.conf
run_stage blockinfo "${species}.bi.csv" wgdi -bi total.conf
run_stage correspondence "${species}.bi.c.csv" wgdi -c total.conf
run_stage blockks "${species}.bk.png" wgdi -bk total.conf
run_stage kspeaks "${species}.kp.png" wgdi -kp total.conf
run_stage peaksfit "${species}.pf.png" wgdi -pf total.conf
EOF

    cat > "$sub_folder/total.conf" <<EOF
[dotplot]
blast = ${species}_${species}.blast
gff1 = ${species}.new.gff
gff2 = ${species}.new.gff
lens1 = ${species}.lens
lens2 = ${species}.lens
genome1_name = ${species}
genome2_name = ${species}
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
savefig = ${species}.png

[collinearity]
blast = ${species}_${species}.blast
gff1 = ${species}.new.gff
gff2 = ${species}.new.gff
lens1 = ${species}.lens
lens2 = ${species}.lens
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
savefile = ${species}.icl

[ks]
cds_file = ${species}.cds
pep_file = ${species}.pep
align_software = muscle
pairs_file = ${species}.icl
ks_file = ${species}.ks

[blockinfo]
blast = ${species}_${species}.blast
gff1 = ${species}.new.gff
gff2 = ${species}.new.gff
lens1 = ${species}.lens
lens2 = ${species}.lens
collinearity = ${species}.icl
score = 100
evalue = 1e-5
repeat_number = 20
position = order
ks = ${species}.ks
ks_col = ks_NG86
savefile = ${species}.bi.csv

[correspondence]
blockinfo = ${species}.bi.csv
lens1 = ${species}.lens
lens2 = ${species}.lens
tandem = true
tandem_length = 200
pvalue = 0.2
block_length = 5
tandem_ratio = 0.5
multiple = 1
homo = -1,1
savefile = ${species}.bi.c.csv

[blockks]
lens1 = ${species}.lens
lens2 = ${species}.lens
genome1_name = ${species}
genome2_name = ${species}
blockinfo = ${species}.bi.csv
pvalue = 0.2
tandem = true
tandem_length = 200
markersize = 1
area = 0,2
block_length = 5
figsize = 8,8
savefig = ${species}.bk.png

[kspeaks]
blockinfo = ${species}.bi.csv
pvalue = 0.2
tandem = true
block_length = 5
ks_area = 0,10
multiple = 1
homo = 0,1
fontsize = 9
area = 0,3
figsize = 10,6.18
savefig = ${species}.kp.png
savefile = ${species}.kp

[peaksfit]
blockinfo = ${species}.bi.csv
mode = median
bins_number = 200
ks_area = 0,10
fontsize = 9
area = 0,3
figsize = 10,6.18
shadow = true
savefig = ${species}.pf.png

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
savefig = ${species}.kf.png
EOF

    find "$sub_folder/.wgdi_steps" -maxdepth 1 -type f -name '*.done' -delete 2>/dev/null || true
    rm -f "$sub_folder/.wgdi.done"
    configured=$((configured + 1))
    echo "已生成同物种配置：$species"
done

echo "同物种 WGDI 配置完成：$configured 个。"
