#!/usr/bin/env python3
"""Exact finite check for teaching; this is not an unbounded proof."""


def check_sum_of_cubes(limit: int = 10_000) -> None:
    if not isinstance(limit, int) or isinstance(limit, bool) or limit < 0:
        raise ValueError("limit must be a nonnegative integer")
    running_sum = 0
    for n in range(limit + 1):
        if n:
            running_sum += n**3
        triangular = n * (n + 1) // 2
        expected = triangular**2
        if running_sum != expected:
            raise AssertionError(
                f"Counterexample at n={n}: sum={running_sum}, formula={expected}"
            )
    print(f"PASS: exact integer equality for n = 0 through {limit:,} ({limit + 1:,} cases).")
    print("Finite computational evidence only; see proof.md for the general argument.")


if __name__ == "__main__":
    check_sum_of_cubes()
