#!/usr/bin/env bash
#
# Build one Docker image per commit in the perf lineage, run the test suites,
# golden scenarios and benchmarks inside each, and write a report.
#
#   docker/lineage.sh            # build + run + report
#   docker/lineage.sh build      # only build images
#   docker/lineage.sh run        # only run (images must exist)
#   docker/lineage.sh report     # only rebuild the report from existing output
#
# Environment:
#   BASE=main                    first commit of the lineage (the baseline)
#   TIP=perf/baseline-tests      last commit of the lineage
#   OUT=lineage-out              where results go (relative to the repo root)
#   IMAGE=hypervolume-tools      image name; tags are :deps and :<short sha>
#   BENCH_TIMES=20               iterations per bench scenario
#   BENCH_ALL=0                  1 = also time commits that don't change package code
#   HYPERVOLUME_TEST_SLOW / HYPERVOLUME_TEST_OMP   passed through to the tests
#   DOCKER=docker                container CLI (podman works)
#
# The harness -- tests/, tools/bench.R and docker/R/ -- is taken from the
# working tree this script runs in, and mounted read-only into every image.

set -euo pipefail

BASE=${BASE:-main}
TIP=${TIP:-perf/baseline-tests}
OUT=${OUT:-lineage-out}
IMAGE=${IMAGE:-hypervolume-tools}
BENCH_TIMES=${BENCH_TIMES:-20}
BENCH_ALL=${BENCH_ALL:-0}
DOCKER=${DOCKER:-docker}
STEP=${1:-all}

repo=$(git rev-parse --show-toplevel)
cd "$repo"
mkdir -p "$OUT"/{tests,golden,bench,logs}
out_abs=$(cd "$OUT" && pwd)

# --- the lineage -------------------------------------------------------------

# commits.tsv: idx, short sha, whether the installed package changed relative
# to the previous commit, subject. Package code = what R CMD INSTALL reads.
# Written by build/run; `report` reuses the one in $OUT so it always matches
# the results it is reporting on.
pkg_paths=(DESCRIPTION NAMESPACE R src man data demo)
write_manifest() {
  printf 'idx\tsha\tcode_changed\tsubject\n'
  prev=""
  i=0
  for sha in $(git rev-parse --short=7 "$BASE") \
             $(git rev-list --reverse --abbrev=7 --abbrev-commit "$BASE..$TIP"); do
    if [[ -z $prev ]]; then
      changed=baseline
    elif git diff --quiet "$prev" "$sha" -- "${pkg_paths[@]}"; then
      changed=no
    else
      changed=yes
    fi
    printf '%02d\t%s\t%s\t%s\n' "$i" "$sha" "$changed" "$(git log -1 --format=%s "$sha")"
    prev=$sha
    i=$((i + 1))
  done
}

# rows of commits.tsv minus the header, as "idx sha code_changed"
lineage() { tail -n +2 "$OUT/commits.tsv" | cut -f1-3; }

# --- build -------------------------------------------------------------------

build() {
  echo "==> deps image $IMAGE:deps"
  "$DOCKER" build -f docker/Dockerfile.deps -t "$IMAGE:deps" .

  local idx sha changed src
  while read -r idx sha changed; do
    echo "==> [$idx] $IMAGE:$sha  $(git log -1 --format=%s "$sha")"
    src=$(mktemp -d)
    git archive "$sha" | tar -x -C "$src"
    # The repo-root .dockerignore is for the deps build (DESCRIPTION only);
    # commits that contain it would otherwise strip their own build context.
    rm -f "$src/.dockerignore"
    "$DOCKER" build -f docker/Dockerfile.commit \
      --build-arg BASE_IMAGE="$IMAGE:deps" \
      --build-arg COMMIT="$(git rev-parse "$sha")" \
      --build-arg SUBJECT="$(git log -1 --format=%s "$sha")" \
      -t "$IMAGE:$sha" "$src"
    rm -rf "$src"
  done < <(lineage)
}

# --- run ---------------------------------------------------------------------

in_image() {
  local tag=$1; shift
  # Run as the calling user so results aren't root-owned on Linux hosts.
  "$DOCKER" run --rm \
    --user "$(id -u):$(id -g)" -e HOME=/tmp \
    -e HV_COMMIT="$tag" -e BENCH_TIMES="$BENCH_TIMES" \
    -e HYPERVOLUME_TEST_SLOW="${HYPERVOLUME_TEST_SLOW:-}" \
    -e HYPERVOLUME_TEST_OMP="${HYPERVOLUME_TEST_OMP:-}" \
    -v "$repo":/harness:ro -v "$out_abs":/out \
    "$IMAGE:$tag" "$@" < /dev/null
}

run() {
  local harness_rev
  harness_rev=$(git rev-parse --short=7 HEAD)
  git diff --quiet HEAD -- tests tools docker || harness_rev="$harness_rev + uncommitted changes"
  {
    echo "harness (tests, bench, scripts): \`$harness_rev\`"
    echo "lineage: \`$BASE\` .. \`$TIP\`"
    echo "run on: $(date '+%Y-%m-%d %H:%M %Z'), $(uname -sm)"
  } > "$OUT/run-info.txt"

  local idx sha changed log
  while read -r idx sha changed; do
    echo "==> [$idx] $sha"
    log="$OUT/logs/$idx-$sha"

    # A crash at one commit is recorded in its log and the lineage carries on.
    if in_image "$sha" test -d /pkg/tests/testthat; then
      in_image "$sha" Rscript /harness/docker/R/lineage_tests.R \
        /pkg/tests/testthat "/out/tests/$idx-$sha-own.csv" own > "$log-tests-own.log" 2>&1 \
        || echo "    own tests crashed, see $log-tests-own.log"
      echo "    $(tail -n 1 "$log-tests-own.log")"
    fi
    in_image "$sha" Rscript /harness/docker/R/lineage_tests.R \
      /harness/tests/testthat "/out/tests/$idx-$sha-tip.csv" tip > "$log-tests-tip.log" 2>&1 \
      || echo "    tip tests crashed, see $log-tests-tip.log"
    echo "    $(tail -n 1 "$log-tests-tip.log")"

    if [[ $changed != no || $BENCH_ALL == 1 ]]; then
      in_image "$sha" Rscript /harness/docker/R/lineage_golden.R \
        /harness "/out/golden/$idx-$sha.rds" > "$log-golden.log" 2>&1 \
        || echo "    golden crashed, see $log-golden.log"
      in_image "$sha" Rscript /harness/docker/R/lineage_bench.R \
        /harness "/out/bench/$idx-$sha.rds" > "$log-bench.log" 2>&1 \
        || echo "    bench crashed, see $log-bench.log"
      echo "    golden + bench done"
    else
      echo "    no package code change; golden + bench skipped"
    fi
  done < <(lineage)
}

# --- report ------------------------------------------------------------------

report() {
  in_image deps Rscript /harness/docker/R/lineage_report.R /out
  echo "==> $OUT/lineage-report.md"
}

case $STEP in
  build|run|all) write_manifest > "$OUT/commits.tsv" ;;
  report) [[ -f $OUT/commits.tsv ]] || { echo "no $OUT/commits.tsv; run first" >&2; exit 1; } ;;
  *) echo "unknown step '$STEP' (build | run | report | all)" >&2; exit 2 ;;
esac

case $STEP in
  build)  build ;;
  run)    run ;;
  report) report ;;
  all)    build; run; report ;;
esac
