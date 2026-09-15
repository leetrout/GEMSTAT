#!/usr/bin/env python3
"""Golden-output regression runner for GEMSTAT.

A *case* is a directory containing ``case.json``::

    {
      "args": ["-s", "@DATA@/seqs.fa", "-fo", "out.txt", "-po", "par.txt", "..."],
      "exit_code": 0,                      # optional, default 0
      "timeout": 300,                      # optional, seconds
      "compare": {                         # files written by the run
        "out.txt": {"rtol": 1e-6, "atol": 1e-9},
        "par.txt": {}
      },
      "stdout_match": ["^Performance"],    # optional: stdout lines to keep and compare
      "stderr_match": [],                  # optional: same for stderr
      "labels": ["fast"]                   # optional: ctest labels
    }

Placeholders in ``args``: ``@CASE@`` (the case directory), ``@ROOT@`` (the
repository root), ``@DATA@`` (the repository ``data/`` directory) and
``@TESTDATA@`` (``tests/data``).

The program runs in a scratch working directory. Every file named in
``compare`` is checked against ``<case>/expected/<name>``. Numeric tokens are
compared with a tolerance, everything else must match exactly. Run with
``--update`` to (re)generate the expected files from the current binary.
"""
import argparse
import json
import math
import os
import re
import shutil
import subprocess
import sys
import tempfile

TOKEN_SPLIT = re.compile(r'[\s,:\[\]{}"]+')
DEFAULT_RTOL = 1e-6
DEFAULT_ATOL = 1e-9


def tokens(text):
    return [t for t in TOKEN_SPLIT.split(text) if t]


def as_float(tok):
    try:
        return float(tok)
    except ValueError:
        return None


def compare_tokens(expected, actual, rtol, atol, label):
    """Return a list of human readable differences (empty when equal)."""
    diffs = []
    e_toks, a_toks = tokens(expected), tokens(actual)
    if len(e_toks) != len(a_toks):
        diffs.append("%s: token count differs: expected %d, got %d"
                     % (label, len(e_toks), len(a_toks)))
    for i, (e, a) in enumerate(zip(e_toks, a_toks)):
        ef, af = as_float(e), as_float(a)
        if ef is not None and af is not None:
            if math.isnan(ef) and math.isnan(af):
                continue
            if math.isinf(ef) or math.isinf(af):
                if ef != af:
                    diffs.append("%s: token %d: expected %s, got %s" % (label, i, e, a))
                continue
            if abs(ef - af) > atol + rtol * abs(ef):
                diffs.append("%s: token %d: expected %s, got %s (rtol=%g atol=%g)"
                             % (label, i, e, a, rtol, atol))
        elif e != a:
            diffs.append("%s: token %d: expected %r, got %r" % (label, i, e, a))
        if len(diffs) >= 10:
            diffs.append("%s: (further differences suppressed)" % label)
            break
    return diffs


def filter_lines(text, patterns):
    regs = [re.compile(p) for p in patterns]
    return "".join(line + "\n" for line in text.splitlines()
                   if any(r.search(line) for r in regs))


def expand(arg, mapping):
    for key, val in mapping.items():
        arg = arg.replace(key, val)
    return arg


def run_case(case_dir, binary, root_dir, data_dir, testdata_dir, workdir, update, keep):
    case_dir = os.path.abspath(case_dir)
    name = os.path.basename(case_dir)
    with open(os.path.join(case_dir, "case.json")) as fh:
        spec = json.load(fh)

    mapping = {"@CASE@": case_dir, "@ROOT@": os.path.abspath(root_dir),
               "@DATA@": os.path.abspath(data_dir),
               "@TESTDATA@": os.path.abspath(testdata_dir)}
    argv = [os.path.abspath(binary)] + [expand(a, mapping) for a in spec["args"]]

    own_tmp = workdir is None
    if own_tmp:
        workdir = tempfile.mkdtemp(prefix="gemstat_case_")
    else:
        shutil.rmtree(workdir, ignore_errors=True)
        os.makedirs(workdir)

    failures = []
    try:
        try:
            proc = subprocess.run(argv, cwd=workdir, capture_output=True, text=True,
                                  timeout=spec.get("timeout", 600))
        except subprocess.TimeoutExpired:
            return ["%s: timed out after %ss" % (name, spec.get("timeout", 600))]

        # The usage message echoes argv[0]; keep the goldens independent of
        # where the binary was built.
        stdout = proc.stdout.replace(argv[0], "seq2expr")
        stderr = proc.stderr.replace(argv[0], "seq2expr")
        with open(os.path.join(workdir, "stdout.txt"), "w") as fh:
            fh.write(stdout)
        with open(os.path.join(workdir, "stderr.txt"), "w") as fh:
            fh.write(stderr)

        want_rc = spec.get("exit_code", 0)
        if proc.returncode != want_rc:
            failures.append("%s: exit code %d, expected %d\n--- stderr tail ---\n%s"
                            % (name, proc.returncode, want_rc, proc.stderr[-2000:]))

        expected_dir = os.path.join(case_dir, "expected")
        checks = []  # (expected file name, actual text, rtol, atol)
        for fname, opts in spec.get("compare", {}).items():
            path = os.path.join(workdir, fname)
            if not os.path.exists(path):
                failures.append("%s: output file %s was not written" % (name, fname))
                continue
            with open(path) as fh:
                checks.append((fname, fh.read(), opts.get("rtol", DEFAULT_RTOL),
                               opts.get("atol", DEFAULT_ATOL)))
        for stream, key in (("stdout", "stdout_match"), ("stderr", "stderr_match")):
            pats = spec.get(key)
            if pats:
                text = stdout if stream == "stdout" else stderr
                checks.append((stream + ".txt", filter_lines(text, pats),
                               DEFAULT_RTOL, DEFAULT_ATOL))

        if update:
            if failures:
                return failures
            os.makedirs(expected_dir, exist_ok=True)
            for fname, text, _, _ in checks:
                with open(os.path.join(expected_dir, fname), "w") as fh:
                    fh.write(text)
            print("%s: updated %d expected file(s)" % (name, len(checks)))
            return []

        for fname, text, rtol, atol in checks:
            epath = os.path.join(expected_dir, fname)
            if not os.path.exists(epath):
                failures.append("%s: no expected file %s (run with --update)" % (name, fname))
                continue
            with open(epath) as fh:
                failures.extend(compare_tokens(fh.read(), text, rtol, atol, name + ":" + fname))
        return failures
    finally:
        if own_tmp and not keep:
            shutil.rmtree(workdir, ignore_errors=True)
        elif failures or keep:
            print("%s: working directory kept at %s" % (name, workdir))


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--binary", required=True, help="path to seq2expr")
    ap.add_argument("--case", action="append", default=[], help="case directory (repeatable)")
    ap.add_argument("--all-cases", help="run every case directory found under this directory")
    ap.add_argument("--root", required=True, help="repository root directory")
    ap.add_argument("--data", required=True, help="repository data/ directory")
    ap.add_argument("--testdata", required=True, help="tests/data directory")
    ap.add_argument("--workdir", help="working directory for a single case (default: temp)")
    ap.add_argument("--update", action="store_true", help="regenerate expected files")
    ap.add_argument("--keep", action="store_true", help="keep the working directory")
    args = ap.parse_args()

    cases = list(args.case)
    if args.all_cases:
        for entry in sorted(os.listdir(args.all_cases)):
            d = os.path.join(args.all_cases, entry)
            if os.path.isfile(os.path.join(d, "case.json")):
                cases.append(d)
    if not cases:
        ap.error("no cases given")
    if args.workdir and len(cases) > 1:
        ap.error("--workdir only makes sense with a single case")

    all_failures = []
    for case in cases:
        fails = run_case(case, args.binary, args.root, args.data, args.testdata, args.workdir,
                         args.update, args.keep)
        if fails:
            all_failures.extend(fails)
            print("FAIL %s" % os.path.basename(case))
        elif not args.update:
            print("ok   %s" % os.path.basename(case))
    if all_failures:
        print("\n".join(all_failures))
        sys.exit(1)


if __name__ == "__main__":
    main()
