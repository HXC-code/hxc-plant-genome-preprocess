#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"
if (($# < 1)); then
    echo "Usage: run_pipeline.sh DATA_ROOT [options]" >&2
    echo "Options: --species LIST --expected-count N --refresh-links --skip-hashes" >&2
    exit 2
fi
if [[ "$1" == "-h" || "$1" == "--help" ]]; then
    echo "Usage: run_pipeline.sh DATA_ROOT [options]"
    echo "Options: --species LIST --expected-count N --refresh-links --skip-hashes"
    exit 0
fi

DATA_ROOT="$1"
shift
SPECIES=""
EXPECTED_COUNT=""
REFRESH_LINKS=false
SKIP_HASHES=false

while (($# > 0)); do
    case "$1" in
        --species)
            (($# >= 2)) || { echo "--species requires a value" >&2; exit 2; }
            SPECIES="$2"
            shift 2
            ;;
        --expected-count)
            (($# >= 2)) || { echo "--expected-count requires a value" >&2; exit 2; }
            EXPECTED_COUNT="$2"
            shift 2
            ;;
        --refresh-links)
            REFRESH_LINKS=true
            shift
            ;;
        --skip-hashes)
            SKIP_HASHES=true
            shift
            ;;
        -h|--help)
            echo "Usage: run_pipeline.sh DATA_ROOT [options]"
            echo "Options: --species LIST --expected-count N --refresh-links --skip-hashes"
            exit 0
            ;;
        *)
            echo "Unknown option: $1" >&2
            exit 2
            ;;
    esac
done

[[ -d "$DATA_ROOT/data/downloads" ]] || {
    echo "Missing data/downloads under: $DATA_ROOT" >&2
    exit 2
}
[[ -f "$DATA_ROOT/resources/download_inventory.tsv" ]] || {
    echo "Missing resources/download_inventory.tsv under: $DATA_ROOT" >&2
    exit 2
}
DATA_ROOT="$(cd -- "$DATA_ROOT" && pwd -P)"
PYTHON_BIN="${PYTHON:-python}"

preprocess_args=(--project-root "$DATA_ROOT")
if [[ -n "$SPECIES" ]]; then
    preprocess_args+=(--species "$SPECIES")
fi
"$PYTHON_BIN" "$SCRIPT_DIR/preprocess.py" "${preprocess_args[@]}"

validate_args=(--project-root "$DATA_ROOT")
if [[ -n "$EXPECTED_COUNT" ]]; then
    validate_args+=(--expected-count "$EXPECTED_COUNT")
fi
if [[ "$REFRESH_LINKS" == true ]]; then
    validate_args+=(--refresh-links)
else
    validate_args+=(--no-compat-links)
fi
if [[ "$SKIP_HASHES" == true ]]; then
    validate_args+=(--skip-hashes)
fi
"$PYTHON_BIN" "$SCRIPT_DIR/validate_dataset.py" "${validate_args[@]}"

