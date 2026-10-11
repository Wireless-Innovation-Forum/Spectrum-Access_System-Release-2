# Spectrum Access System Release 2 Testing and Interoperability Repository

This repository contains code and data for testing the compliance of Spectrum 
Access System (SAS) software. The SAS is defined by the FCC in proceeding 
12-354 as the system authorizing priority access and general access to the 
3550-3700MHz Citizens Broadband Radio Service.

This repository contains procedures, documentation, and tests for such software, 
and for the devices authorized by it.

To contribute, please first read the CONTRIBUTING file in the repository for 
instructions.

## Implementation status and Python compatibility

This repository currently contains Release 2 alignment documentation and one
testcase placeholder containing only a license header. It does not yet contain
an executable Release 2 harness or implemented compliance tests.

New Python code should support Python 3.11 through 3.15. CI checks the existing
source syntax on each of these versions; that check does not certify a SAS or
validate an implemented Release 2 test procedure. There are no third-party or
cryptographic dependencies in this repository.

Use [uv](https://docs.astral.sh/uv/) to create the pinned Python environment
and check the source from the repository root:

```sh
uv sync --locked
uv run --no-sync python -m compileall -q src
```

`.python-version` selects Python 3.15.0. To use an older supported version,
pass the same `--python` selection to both commands, for example
`uv sync --locked --python 3.11` and
`uv run --no-sync --python 3.11 python -m compileall -q src`.
This is a uv project without an installable package or runtime dependencies.

Release 2 development reuses the sibling
[Spectrum-Access-System](https://github.com/Wireless-Innovation-Forum/Spectrum-Access-System)
Release 1 harness and the shared
[Common-Data](https://github.com/Wireless-Innovation-Forum/Common-Data) repository.
See [the harness instructions](src/harness/README.md) for the workspace layout
and validation commands.
