---
name: run-tests
description: Run or add analytical-solution tests for ROMS structure drag/production and body force. Use when validating model or template changes.
---

# Tests

```bash
python tests/run_tests.py                          # use existing output (runs model if missing)
python tests/run_tests.py --run-model              # fresh model run first
python tests/run_tests.py --tests test_UV_BODYFORCE
```

| Test | Check |
|---|---|
| `test_UV_BODYFORCE` | u = F·t, rtol 1e-5 |
| `test_STRUCTURE_DRAG` | u = √(F/α)·tanh(t√(Fα)), rtol 1e-4 |
| `test_STRUCTURE_PRODUCTION` | integrated ε matches Pd within 2% |

Adding a test: create `tests/test_<NAME>.py` (auto-discovered) and a matching `configs/variants/test_<NAME>.yaml`; compare against an analytical solution.
