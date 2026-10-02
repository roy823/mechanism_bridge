#!/usr/bin/env bash
# Run inside an allocated compute job with the project environment activated.
# Usage: bash run_pytest_comparison.sh NEW_OUTPUT_DIR TREE [TREE ...]
# Run every tree, retain its exit code, and fail if any suite fails.
set -eo pipefail

if [ "$#" -lt 2 ]; then
    echo "Usage: $0 NEW_OUTPUT_DIR TREE [TREE ...]" >&2
    exit 2
fi
output=$1
shift
for tree in "$@"; do
    test -d "$tree/tests" || { echo "Missing tests directory: $tree/tests" >&2; exit 2; }
    git -C "$tree" rev-parse --verify HEAD >/dev/null
done
# Refuse to overwrite evidence from an earlier job.
mkdir "$output"
output=$(cd "$output" && pwd -P)
printf 'index\tcommit\tpytest_exit\ttree\n' > "$output/results.tsv"
overall=0
index=0
for tree in "$@"; do
    index=$((index + 1))
    tree=$(cd "$tree" && pwd -P)
    commit=$(git -C "$tree" rev-parse HEAD)
    git -C "$tree" status --porcelain --untracked-files=normal -- \
        src scripts configs tests pyproject.toml > "$output/tree_${index}.status"
    git -C "$tree" diff HEAD --binary -- \
        src scripts configs tests pyproject.toml > "$output/tree_${index}.patch"
    printf 'Testing %s at %s\n' "$tree" "$commit"
    if (cd "$tree" && python -m pytest -q -rfE -p no:cacheprovider \
        --junitxml="$output/tree_${index}.xml" tests) > "$output/tree_${index}.log" 2>&1; then
        code=0
    else
        code=$?
        overall=1
    fi
    cat "$output/tree_${index}.log"
    printf '%s\t%s\t%s\t%s\n' "$index" "$commit" "$code" "$tree" >> "$output/results.tsv"
    printf 'pytest exit: %s\n' "$code"
done
cat "$output/results.tsv"
exit "$overall"
