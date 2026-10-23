# random: Random Numbers

*Status: draft (2026-10-23).*

## 1. Purpose

Simulations draw millions of random numbers in parallel and must give the same result for the same seed, whatever the number of threads and the order of work. A counter-based generator makes every draw a pure function of (key, counter), so a draw is identified by what it is for (a purpose, an object, an entity, a step) rather than by who makes it.

## 2. Content

| Header | Content |
|---|---|
| `philox.hpp` | Philox4x32-10 (Salmon et al. 2011), evaluable at compile time |
| `stream.hpp` | `Stream(seed, a, b, c)`: three identity words of the counter and a walking block; uniform on (0, 1) with 53 bits, bounded integers (Lemire 2019), normal (Box–Muller), lognormal, exponential, geometric, Bernoulli, Poisson (multiplication below a mean of 10, PTRS above) |

## 3. Checks

The known-answer vectors of Random123; the moments of each distribution over a million draws, as z-scores; a chi-square of bounded integers; the independence of identities and seeds; `seek`.

## 4. Dependencies

None.

## 5. References

- Salmon, J. K., Moraes, M. A., Dror, R. O. and Shaw, D. E. (2011), Parallel random numbers: as easy as 1, 2, 3, *Proc. SC11*, 16:1–16:12.
- Lemire, D. (2019), Fast random integer generation in an interval, *ACM TOMACS* 29, 3:1–3:12.
- Box, G. E. P. and Muller, M. E. (1958), A note on the generation of random normal deviates, *Ann. Math. Stat.* 29, 610–611.
- Hörmann, W. (1993), The transformed rejection method for generating Poisson random variables, *Insurance: Mathematics and Economics* 12, 39–45.

## 6. Open Questions

1. Further distributions (binomial, gamma, beta, multivariate normal) as callers need them.
2. Quasi-random sequences (Sobol, Halton) for quadrature in many dimensions.
