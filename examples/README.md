# examples/

Worked examples for the Repo Rules in `../CLAUDE.md`.

| File | Rule(s) | Shows |
| --- | --- | --- |
| `type_hints_demo.py` | 2 | Annotating scalars, containers, `Optional`/union, `Callable`, `*args`, keyword-only args, and dataclass fields. |
| `running_stats.py` | 1, 2, 3 | The unit under test: terse academic docstrings, full type hints, composition via `merge`/`from_iterable`. |
| `test_running_stats.py` | 4, 5 | pytest structure — test classes by scope, terse docstrings, fixtures, `parametrize`, `approx`. |

## Running the tests

```bash
uv run pytest examples/
```

## pytest in 60 seconds

- **Discovery.** pytest collects files named `test_*.py`, then classes
  named `Test*` (no `__init__`), then functions/methods named `test_*`.
- **Assertions.** Use a bare `assert`; pytest rewrites it to show both
  operands on failure. No `self.assertEqual`.
- **Fixtures.** A function decorated with `@pytest.fixture` is requested
  by naming it as a test argument. Its return value is injected; code
  after a `yield` runs as teardown. Put shared fixtures in `conftest.py`
  to share them across files without importing.
- **Parametrize.** `@pytest.mark.parametrize("name", [v1, v2, ...])` runs
  the test once per value, each as its own reported node.
- **Floats.** Compare with `pytest.approx(expected)`; never `==`.
- **Expected exceptions.** `with pytest.raises(ValueError): ...`.
- **Node ids.** A failure prints as
  `path::TestClass::test_method[param]` — which is why classes and test
  names should read as a sentence about the behaviour under test.
