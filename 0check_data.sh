#!/usr/bin/env bash
# Thin orchestration wrapper for the Python preprocessor and validator.
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"
PROJECT_ROOT="$SCRIPT_DIR"
PROJECT_ROOT_EXPLICIT=false
EXPECTED_COUNT=""
PREPROCESS=false
PREPROCESS_SPECIES=""
REFRESH_LINKS=false
REQUIRE_COMPAT=true
VERIFY_HASHES=true
POSITIONAL_ROOT=""

usage() {
    cat <<'EOF'
Usage:
  bash 0check_data.sh [PROJECT_ROOT] [options]

Options:
  --project-root PATH    Project root containing plant-genome-resources/
  --expected-count N     Fail unless the manifest contains exactly N species
  --preprocess           Run the Python standardization pipeline first
  --species LIST         With --preprocess, process comma-separated show_ids
  --refresh-links        Refresh legacy root s<show_id> compatibility links
  --no-compat-links      Validate canonical analysis_v2 data only
  --skip-hashes          Skip output SHA-256 verification
  -h, --help             Show this help

Exit codes:
  0 all requested steps and validation passed
  1 preprocessing or validation failed
  2 project layout or arguments are invalid
  3 Python is unavailable
EOF
}

die_usage() {
    printf '错误: %s\n' "$1" >&2
    usage >&2
    exit 2
}

while (($# > 0)); do
    case "$1" in
        --project-root)
            (($# >= 2)) || die_usage "--project-root 缺少参数"
            PROJECT_ROOT="$2"
            PROJECT_ROOT_EXPLICIT=true
            shift 2
            ;;
        --expected-count)
            (($# >= 2)) || die_usage "--expected-count 缺少参数"
            EXPECTED_COUNT="$2"
            shift 2
            ;;
        --preprocess)
            PREPROCESS=true
            shift
            ;;
        --species)
            (($# >= 2)) || die_usage "--species 缺少参数"
            PREPROCESS_SPECIES="$2"
            shift 2
            ;;
        --refresh-links)
            REFRESH_LINKS=true
            shift
            ;;
        --no-compat-links)
            REQUIRE_COMPAT=false
            shift
            ;;
        --skip-hashes)
            VERIFY_HASHES=false
            shift
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        --*)
            die_usage "未知参数: $1"
            ;;
        *)
            [[ -z "$POSITIONAL_ROOT" ]] || die_usage "只能指定一个位置参数"
            POSITIONAL_ROOT="$1"
            shift
            ;;
    esac
done

if [[ -n "$POSITIONAL_ROOT" ]]; then
    if [[ "$PROJECT_ROOT_EXPLICIT" == true ]]; then
        die_usage "--project-root 与位置参数不能同时使用"
    fi
    PROJECT_ROOT="$POSITIONAL_ROOT"
fi

[[ -d "$PROJECT_ROOT" ]] || {
    printf '错误: 项目目录不存在: %s\n' "$PROJECT_ROOT" >&2
    exit 2
}
PROJECT_ROOT="$(cd -- "$PROJECT_ROOT" && pwd -P)"

if [[ -f "$SCRIPT_DIR/plant-genome-resources/scripts/validate_wgdi_dataset.py" ]]; then
    TOOL_ROOT="$SCRIPT_DIR/plant-genome-resources"
elif [[ -f "$SCRIPT_DIR/scripts/validate_wgdi_dataset.py" ]]; then
    TOOL_ROOT="$SCRIPT_DIR"
else
    printf '错误: 未找到工具目录 scripts/validate_wgdi_dataset.py\n' >&2
    exit 2
fi

if [[ -f "$PROJECT_ROOT/plant-genome-resources/resources/wgdi_v2_index.tsv" ]]; then
    TARGET_CANONICAL_ROOT="$PROJECT_ROOT/plant-genome-resources"
elif [[ -f "$PROJECT_ROOT/resources/wgdi_v2_index.tsv" ]]; then
    TARGET_CANONICAL_ROOT="$PROJECT_ROOT"
else
    printf '错误: 目标目录缺少 resources/wgdi_v2_index.tsv\n' >&2
    exit 2
fi

PYTHON_BIN="${PYTHON:-python}"
if ! command -v "$PYTHON_BIN" >/dev/null 2>&1; then
    printf '错误: 找不到 Python: %s\n' "$PYTHON_BIN" >&2
    exit 3
fi

if [[ -n "$EXPECTED_COUNT" && ! "$EXPECTED_COUNT" =~ ^[1-9][0-9]*$ ]]; then
    die_usage "--expected-count 必须是正整数"
fi
if [[ "$REFRESH_LINKS" == true && "$REQUIRE_COMPAT" == false ]]; then
    die_usage "--refresh-links 与 --no-compat-links 不能同时使用"
fi
if [[ -n "$PREPROCESS_SPECIES" && "$PREPROCESS" == false ]]; then
    die_usage "--species 只能与 --preprocess 一起使用"
fi

if [[ "$PREPROCESS" == true ]]; then
    [[ -f "$TARGET_CANONICAL_ROOT/scripts/preprocess_wgdi_v2.py" ]] || {
        printf '错误: 目标目录缺少 scripts/preprocess_wgdi_v2.py\n' >&2
        exit 2
    }
    echo "运行数据预处理: $TARGET_CANONICAL_ROOT"
    (
        cd -- "$TARGET_CANONICAL_ROOT"
        preprocess_args=()
        if [[ -n "$PREPROCESS_SPECIES" ]]; then
            preprocess_args+=(--species "$PREPROCESS_SPECIES")
        fi
        "$PYTHON_BIN" scripts/preprocess_wgdi_v2.py "${preprocess_args[@]}"
    )
fi

validator_args=(
    --project-root "$PROJECT_ROOT"
)
if [[ -n "$EXPECTED_COUNT" ]]; then
    validator_args+=(--expected-count "$EXPECTED_COUNT")
fi
if [[ "$REFRESH_LINKS" == true ]]; then
    validator_args+=(--refresh-links)
fi
if [[ "$REQUIRE_COMPAT" == false ]]; then
    validator_args+=(--no-compat-links)
fi
if [[ "$VERIFY_HASHES" == false ]]; then
    validator_args+=(--skip-hashes)
fi

echo "验收项目: $PROJECT_ROOT"
"$PYTHON_BIN" "$TOOL_ROOT/scripts/validate_wgdi_dataset.py" "${validator_args[@]}"
