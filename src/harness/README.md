# The SAS Test Harness

This folder contains code and testcases for SAS certification test as defined in
"CBRS Release 2 Regulatory Impacting (RI) Test Specification for SAS as UUT" (WINNF-TS-4006).

## Code prerequisites

Use Python 3.11 or later, targeting Python 3.15. This repository has no
third-party dependencies because its only Python file is currently a testcase
placeholder. It has no standalone harness entry point.

Keep the collaborating repositories beside each other:

```text
GitRepos/
    Citizens-Broadband-Radio-Service-Device/
    Spectrum-Access-System/
    Spectrum-Access_System-Release-2/
    Common-Data/
```

For Release 1 harness development, install and configure the sibling
`Spectrum-Access-System` checkout using its own `README.md` and
`src/harness/README.md`. Its requirements, native propagation extensions,
certificate configuration, and shared-data setup apply when reusing that
harness.

After configuring `sas.cfg` and `test.cfg` in the Release 1 harness, its
operational entry point is:

```sh
cd ../Spectrum-Access-System/src/harness
uv run --no-sync python test_main.py
```

This runs the Release 1 harness against the configured SAS under test. It
does not discover or run the unimplemented Release 2 placeholder.

From the root of this repository, check Python syntax with:

```sh
uv sync --locked
uv run --no-sync python -m compileall -q src
```

CI runs this check on Python 3.11, 3.12, 3.13, 3.14, and 3.15. There are no
implemented Release 2 compliance tests to execute yet; source compilation
alone does not validate behavior or SAS interoperability. Implemented
testcases will need functional tests and explicit integration with the
Release 1 harness before that claim can be made.

## Alignment to Release 1 Testcases

To the maximum extent possible, testcase code from 
[Release 1 Spectrum Access System](https://github.com/Wireless-Innovation-Forum/Spectrum-Access-System)
will be reused in developing Release 2 testcases.

Alignment assessments for each test case are provided in:

    ./rel1_alignment

In most cases, there is not a perfect alignment between Release 1 testcases and Release 2 requirements.

## Code Location

* **./**: Test Framework
* **./rel1_alignment**: Testcase alignment assessment
* **./testcases**: Test cases, grouped by section in the test specification.
  `WINNF_PT_S_REL2_RI_GRU_testcase.py` is currently an empty placeholder.
