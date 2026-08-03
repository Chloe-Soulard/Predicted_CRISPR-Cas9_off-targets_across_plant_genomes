"""
Step 5 driver — run `step5_crispor.py --all --resume` until nothing is pending.

CRISPOR queues jobs, so a single pass leaves some windows "pending". This loop
re-runs the step until every window is resolved (or a safety cap of passes is
reached), waiting out any pass already in progress so it never clashes with a
manual run.

Run detached and check the log:
    python step5_loop.py > step5_loop.log 2>&1
"""

from __future__ import annotations

import argparse
import subprocess
import sys
import time

import crispor
from paths import DEFAULT_PAM, ROOT, load_json, step5_dir

MAX_PASSES      = 40
SLEEP_BETWEEN_S = 900      # 15 min between passes: CRISPOR jobs take minutes
LOCK_POLL_S     = 60


def status_counts(directory) -> dict[str, int]:
    """Tally window statuses across every gene's step 5 results."""
    counts: dict[str, int] = {}
    for path in sorted(directory.glob("*.json")):
        for record in load_json(path).values():
            status = record.get("status", "unknown")
            counts[status] = counts.get(status, 0) + 1
    return counts


def describe(counts: dict[str, int]) -> str:
    return " ".join(f"{k}:{v}" for k, v in sorted(counts.items())) or "no results yet"


def log(message: str) -> None:
    print(time.strftime("%H:%M:%S"), message, flush=True)


def wait_for_lock(directory) -> None:
    """Block while another CRISPOR pass holds the write lock."""
    while crispor.lock_holder(directory) is not None:
        log("another CRISPOR pass is running — waiting ...")
        time.sleep(LOCK_POLL_S)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.strip().splitlines()[0])
    parser.add_argument("--max-passes", type=int, default=MAX_PASSES,
                        help=f"safety cap on the number of passes (default {MAX_PASSES})")
    parser.add_argument("--sleep", type=int, default=SLEEP_BETWEEN_S,
                        help=f"seconds between passes (default {SLEEP_BETWEEN_S})")
    parser.add_argument("--pam", default=DEFAULT_PAM,
                        help=f"CRISPOR PAM code to drive (default {DEFAULT_PAM})")
    args = parser.parse_args()

    pam = args.pam.strip().upper()
    directory = step5_dir(pam)
    directory.mkdir(parents=True, exist_ok=True)
    log(f"driving PAM {pam} -> {directory.name}")

    for pass_number in range(1, args.max_passes + 1):
        wait_for_lock(directory)

        before = status_counts(directory)
        log(f"--- pass {pass_number}: before = {describe(before)}")
        if pass_number > 1 and not before.get("pending"):
            log("nothing pending — done")
            return

        subprocess.run([sys.executable, str(ROOT / "step5_crispor.py"),
                        "--all", "--resume", "--pam", pam],
                       cwd=str(ROOT), check=False)

        after = status_counts(directory)
        log(f"--- pass {pass_number}: after  = {describe(after)}")
        if not after.get("pending"):
            log("all windows resolved — done")
            return

        if pass_number < args.max_passes:
            log(f"sleeping {args.sleep // 60} min before the next pass ...")
            time.sleep(args.sleep)

    log(f"reached the {args.max_passes}-pass cap with windows still pending")


if __name__ == "__main__":
    main()
