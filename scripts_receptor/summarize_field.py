#!/usr/bin/env python3
"""Summarize field_diff.csv by the difference in total charge between the two tools.

  python3 summarize_field.py <field_diff.csv>

measure_field.py prints the median, quartiles and range over the whole table. This script adds what
the manuscript reports beyond that:

  - the number of structures in which the two totals differ, the largest difference, and its sign
  - the median RMS of the potential difference in classes of the absolute difference in total charge
    (0, 1, 2 to 4, 5 or more)
  - the Spearman rank correlation between the absolute difference in total charge and the RMS
  - for the structures whose totals agree, the quartiles of the RMS and the median and range of the
    largest difference at any surface point

d_sum_q is rounded to the nearest integer before it is classed. Both tools return an integral total,
and the stored values carry floating-point noise.

The Spearman coefficient is computed here, with average ranks for ties, so that no package beyond the
standard library is needed.
"""

import csv
import statistics
import sys


def ranks(v):
    """Average ranks, 1-based, ties sharing the mean of their positions."""
    order = sorted(range(len(v)), key=lambda i: v[i])
    r = [0.0] * len(v)
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and v[order[j + 1]] == v[order[i]]:
            j += 1
        for k in range(i, j + 1):
            r[order[k]] = (i + j) / 2 + 1
        i = j + 1
    return r


def pearson(a, b):
    ma, mb = statistics.fmean(a), statistics.fmean(b)
    num = sum((x - ma) * (y - mb) for x, y in zip(a, b, strict=True))
    den = (sum((x - ma) ** 2 for x in a) * sum((y - mb) ** 2 for y in b)) ** 0.5
    return num / den


def quart(v):
    q = statistics.quantiles(v, n=4, method="inclusive")
    return q[0], q[2]


def main():
    if len(sys.argv) != 2 or sys.argv[1] in ("-h", "--help"):
        sys.exit(__doc__)
    rows = list(csv.DictReader(open(sys.argv[1])))
    dq = [round(float(r["d_sum_q"])) for r in rows]
    rms = [float(r["rms_dv"]) for r in rows]
    mx = [float(r["max_abs_dv"]) for r in rows]
    n = len(rows)

    differ = [d for d in dq if d]
    print(f"{n} structures")
    print(f"total charge differs in {len(differ)}/{n}, by up to {max(abs(d) for d in differ)}")
    print(f"  the first tool's total is the more negative in {sum(1 for d in differ if d < 0)},"
          f" the more positive in {sum(1 for d in differ if d > 0)}")

    classes = (("0", lambda d: d == 0), ("1", lambda d: d == 1), ("2 to 4", lambda d: 2 <= d <= 4),
               ("5 or more", lambda d: d >= 5))
    print("\n|difference in total charge|   structures   median RMS of dV")
    for label, test in classes:
        v = [rms[i] for i in range(n) if test(abs(dq[i]))]
        print(f"  {label:<26}{len(v):>10}{statistics.median(v):>15.2f}")

    rho = pearson(ranks([abs(d) for d in dq]), ranks(rms))
    print(f"\nSpearman correlation between |difference in total charge| and RMS of dV: {rho:.3f}")

    same = [i for i in range(n) if dq[i] == 0]
    r0 = [rms[i] for i in same]
    m0 = [mx[i] for i in same]
    q1, q3 = quart(r0)
    print(f"\n{len(same)} structures whose totals agree")
    print(f"  RMS of dV: median {statistics.median(r0):.2f}, quartiles {q1:.2f} to {q3:.2f}")
    print(f"  maximum of dV: median {statistics.median(m0):.2f}, range {min(m0):.2f} to {max(m0):.2f}")


if __name__ == "__main__":
    main()
