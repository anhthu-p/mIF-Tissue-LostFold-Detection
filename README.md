# HTAN_Breast_Cancer_cyCIF_2026

Spatial proteomics (cyCIF) analysis of the tumor microenvironment in Invasive
Lobular Carcinoma (ILC) vs. Invasive Ductal Carcinoma/No-Special Type (NST)
breast cancer, correlated against tumor proliferation.

STEP 1: AUTOMATED 



## uv cheat sheet

[uv](https://docs.astral.sh/uv/) manages the Python version, virtual
environment, and dependencies together.

```bash
# Create/refresh the .venv from pyproject.toml + uv.lock
uv sync

# Add a runtime dependency (writes to pyproject.toml and uv.lock)
uv add <package>

# Add a dev-only dependency
uv add --dev <package>

# Remove a dependency
uv remove <package>

# Run a command inside the project's venv without activating it
uv run <command>
uv run python some_script.py
uv run pytest examples/

# Activate the venv directly (optional; uv run usually makes this unnecessary)
source .venv/bin/activate

# Upgrade all locked dependencies
uv lock --upgrade

# Pin/change the Python version this project uses
uv python pin 3.12
```

`.python-version` pins the interpreter version; `uv.lock` pins exact
dependency versions and should be committed. `.venv/` is git-ignored and
disposable — delete it and re-run `uv sync` any time.
