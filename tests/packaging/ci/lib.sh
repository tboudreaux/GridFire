# Shared helpers for the packaging CI scripts (source this file). Every check is run, whatever happened to the ones
# before it, and its outcome is recorded so that the workflow can report which checks pass and which fail:
#
#   ci_init <platform> <output dir> [job] results go to <output dir>/results-<platform>[-<job>].json (jobs on the same
#                                         platform are merged by report.py), logs to <output dir>/logs
#   ci_check <name> <command...>          run a check; record pass or fail
#   ci_skip <name> <reason>               record a check that could not run (e.g. the wheel it needs failed to build)
#   ci_finish                             print the table, write the JSON, exit non-zero if anything failed
#
# On GitHub Actions each check's output is a collapsible group, and the table is added to the job summary.

CI_PLATFORM=""
CI_OUT=""
CI_RESULTS=""
CI_NAMES=()
CI_STATUSES=()
CI_SECONDS=()

ci_init() {
    CI_PLATFORM="$1"
    CI_OUT="$2"
    CI_RESULTS="${CI_OUT}/results-${CI_PLATFORM}${3:+-$3}.json"
    mkdir -p "${CI_OUT}/logs"
}

_ci_record() {
    CI_NAMES+=("$1")
    CI_STATUSES+=("$2")
    CI_SECONDS+=("$3")
}

ci_check() {
    local name="$1"
    shift
    local log="${CI_OUT}/logs/${CI_PLATFORM}-${name}.log"
    local start=${SECONDS}
    [[ -n "${GITHUB_ACTIONS:-}" ]] && echo "::group::${name}"
    echo "=== ${name}: $*"
    # In a subshell, so a failing `set -e` or `exit` inside the check stops only that check
    if ( set -o pipefail; "$@" 2>&1 | tee "${log}" ); then
        _ci_record "${name}" pass $((SECONDS - start))
        [[ -n "${GITHUB_ACTIONS:-}" ]] && echo "::endgroup::"
        echo "=== ${name}: PASS ($((SECONDS - start)) s)"
    else
        _ci_record "${name}" fail $((SECONDS - start))
        [[ -n "${GITHUB_ACTIONS:-}" ]] && echo "::endgroup::"
        echo "=== ${name}: FAIL ($((SECONDS - start)) s); last lines of ${log}:"
        tail -40 "${log}"
        [[ -n "${GITHUB_ACTIONS:-}" ]] && echo "::error title=${CI_PLATFORM} ${name}::${name} failed on ${CI_PLATFORM}"
    fi
    return 0
}

ci_skip() {
    _ci_record "$1" skipped 0
    echo "=== $1: SKIPPED ($2)"
}

# True if the named check has been recorded as passing
ci_passed() {
    local i
    for i in "${!CI_NAMES[@]}"; do
        [[ "${CI_NAMES[$i]}" == "$1" && "${CI_STATUSES[$i]}" == pass ]] && return 0
    done
    return 1
}

ci_finish() {
    local i failed=0 json table
    table="| Check | Result | Time |"$'\n'"|---|---|---|"
    json="{\"platform\": \"${CI_PLATFORM}\", \"checks\": {"
    for i in "${!CI_NAMES[@]}"; do
        [[ "${CI_STATUSES[$i]}" == pass ]] || failed=1
        table+=$'\n'"| ${CI_NAMES[$i]} | ${CI_STATUSES[$i]} | ${CI_SECONDS[$i]} s |"
        [[ $i -gt 0 ]] && json+=", "
        json+="\"${CI_NAMES[$i]}\": \"${CI_STATUSES[$i]}\""
    done
    json+="}}"
    echo "${json}" > "${CI_RESULTS}"
    echo
    echo "Packaging checks on ${CI_PLATFORM}:"
    echo "${table}"
    if [[ -n "${GITHUB_STEP_SUMMARY:-}" && -w "$(dirname "${GITHUB_STEP_SUMMARY}")" ]]; then
        printf '### Packaging checks: %s\n\n%s\n\n' "${CI_PLATFORM}" "${table}" >> "${GITHUB_STEP_SUMMARY}"
    fi
    return ${failed}
}
