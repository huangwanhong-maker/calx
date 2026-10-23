# verification: Checks of Numerical Results

*Status: draft (2026-10-23).*

## 1. Purpose

A numerical result is believed when it agrees with a reference within a stated limit. This module states that agreement and reports it in the form every check of The Science Infrastructures prints, whatever its language, so that one runner adds up the checks of all the libraries:

```
ok    name: 1.030e-2 unit (limit 0.02 unit)
FAIL  name: detail
skip  what

41 of 41 checks passed.
```

## 2. Content

- `Suite`: `check(name, value, limit)` passes when the value is finite and its magnitude within the limit; `require(name, condition)`; `skip(what)`; `done()` prints the count and returns the exit status.
- Measures of error: `relative`, `max_difference`, and `order`, the observed order of convergence from errors at two steps.
- CMake: `calx_add_check(<name> SOURCES ... LIBRARIES ...)` makes a check program and registers it with CTest; physis and the domains use it too.

## 3. Dependencies

None.

## 4. Open Questions

1. Checks of floating-point results across backends, with the tolerance each computation declares (calx invariant 3).
