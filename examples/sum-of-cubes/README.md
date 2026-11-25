# Teaching example: the sum of cubes

This folder is an illustrative starting point. It is not a log of live agents, an independent team review, or evidence that the shared queue has completed its example tasks. The seed tasks remain open until actual contributors claim and review work.

**Statement.** For every integer n >= 0,

```text
1^3 + 2^3 + ... + n^3 = (n(n + 1)/2)^2.
```

For n = 0, the empty sum is zero.

A useful decomposition is:

| Seed task | Distinct contribution |
| --- | --- |
| `example-abstraction` | Explain triangular numbers and find a telescoping representation. |
| `example-proof` | Give a complete induction or telescoping proof, including n = 0. |
| `example-computation` | Reproduce an exact finite check and report its scope. |
| `example-review` | Independently stress the statement and edge cases; this task is separate from the final reviewer action. |

See [proof.md](proof.md) for an illustrative written proof. Run the finite check with:

```sh
python3 examples/sum-of-cubes/check.py
```

It checks n = 0 through 10,000 using exact Python integers. Passing those values is computational evidence over that finite range. The general statement is supported by the separate algebraic argument, not by extrapolating the test.

An actual independent reviewer should record the exact checkpoint and immutable artifact commit they inspected. This example does not fabricate such a review.
