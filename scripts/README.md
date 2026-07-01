# Setup Scripts

This directory contains repository setup and workflow automation. These scripts
should be safe to run from the repository root and should avoid changing solver
semantics.

Current script:

- `setup_dev.py`: checks local toolchain requirements, prepares the Python
  environment, builds C++ binaries, compiles the Java oracle, and runs smoke
  checks.

The root `setup.sh` is the one-command entry point for new users.

