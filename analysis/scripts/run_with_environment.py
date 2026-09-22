#!/usr/bin/env python
"""Run a script with an explicit compatible installed package directory.

This recovery launcher does not rebuild a virtual environment. Use a working
Python with the same version/architecture as the installed binary packages.
"""
from __future__ import annotations

import argparse
import importlib
import json
import runpy
import sys
import sysconfig
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--site-packages", type=Path, required=True)
    parser.add_argument("--check", action="store_true", help="Import core numerical packages and report versions")
    parser.add_argument("script", nargs="?")
    parser.add_argument("args", nargs=argparse.REMAINDER)
    args = parser.parse_args()
    directory = args.site_packages.resolve()
    if not directory.is_dir():
        parser.error(f"Missing package directory: {directory}")
    sys.path.insert(0, str(directory))
    if args.check:
        versions = {name: importlib.import_module(name).__version__
                    for name in ("numpy", "scipy", "pandas", "h5py", "anndata")}
        print(json.dumps({"python": sys.version, "executable": sys.executable,
                          "platform": sysconfig.get_platform(), "packages": versions}, indent=2))
    if args.script:
        script = Path(args.script).resolve()
        sys.path.insert(0, str(script.parent))
        sys.argv = [str(script), *args.args]
        runpy.run_path(str(script), run_name="__main__")
    elif not args.check:
        parser.error("Provide a script or --check")


if __name__ == "__main__":
    main()
