# Sum of cubes: illustrative written proof

**Evidence level:** complete elementary written argument supplied as teaching material. No live team review is claimed.

**Statement:** For every integer n >= 0, the sum of k^3 for k = 1, ..., n equals (n(n + 1)/2)^2. The sum is zero when n = 0.

Define the triangular number T_k = k(k + 1)/2 for integers k >= 0. For k >= 1,

```text
T_k^2 - T_(k-1)^2
= (T_k - T_(k-1))(T_k + T_(k-1))
= k * k^2
= k^3.
```

The middle equality uses

```text
T_k - T_(k-1) = [k(k + 1) - (k - 1)k] / 2 = k,
T_k + T_(k-1) = [k(k + 1) + (k - 1)k] / 2 = k^2.
```

For n >= 1, summing the first identity from k = 1 to n telescopes:

```text
sum(k^3, k = 1..n)
= sum(T_k^2 - T_(k-1)^2, k = 1..n)
= T_n^2 - T_0^2
= [n(n + 1)/2]^2.
```

For n = 0, both sides are zero. This covers the entire stated domain. No numerical approximation, infinite series, or convergence assumption is used.

The accompanying program checks a finite set of instances with exact integer arithmetic. It is a reproducibility example and a useful guard against transcription errors; it is not the reason the universal claim holds.
