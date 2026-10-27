# DL_Summer_2026

## GOAL
I am a Ph.D student in Cancer Biology. My background is more lab work and I have been coding for a year, mainly for exploratory for bioinformatics tasks so I haven't had a full project at individual industry level before. I want to learn how to do that with the high priority goal to publish a reusable project on github and use it for my internship application this upcoming fall 2026/spring 2027. I want you to teach me new things so that I can be proficient and competitive for the internship application, I don't want you to code for me unless I ask you to help me directly.

For this project, I want to use analyze the tumor microenvironment using spatial proteomics (method: cyCIF) in Invasive Lobular Carcinoma (ILC) and Invasive Ductal Carcinoma/ No-Special Type (NST) to see how it correlates with the tumor proliferation. 

- Trying to learn PyTorch


## Data
- Location for raw images:
* ILC: "D:\thu\HTAN\images\ILC\level_2"
* NST: "D:\thu\HTAN\images\NST"
- Format: ome.tif
- cyCIF

## Stack
- Python 
- Trying to learn PyTorch
- Dependency/venv management via `uv` (see `pyproject.toml`, `uv.lock`, `.python-version`)

## Repo Rules

1. Always adhere to best code practices for Python. Write terse docstrings for methods and properties in an academic, informative tone, using the same format shown in `examples/running_stats.py`.

2. Always use type hints for function arguments wherever possible.

3. Prefer composition to inheritance wherever possible.

4. Write tests that capture fundamental behavior using pytest, so regressions get caught before they become a problem.

5. Tests should always contain terse docstrings, and should be grouped into test classes by code proximity and test scope wherever possible.

6. Never attempt to prove assertions with random tests or experiments without validating them with me first. No examples.

7. Code should always remain tightly scoped, with single responsibility as much as possible. If scope begins to creep, suggest pausing to isolate scope into small, testable components before generating massive scripts. Keep scripts under 700-800 lines of code in general.

See `examples/` for worked examples of rules 1, 2, 4, and 5.
