# SPDX-License-Identifier: MIT
"""Package-level tests: version string and import hygiene."""

import subprocess
import sys


def test_version():
    import aeroprops

    assert isinstance(aeroprops.__version__, str)
    assert aeroprops.__version__.count(".") >= 1


def test_import_hygiene():
    """Importing aeroprops must pull in neither RDKit nor thermo.

    Pins the layering discipline: the structure layer loads RDKit lazily
    and nothing in the package top level needs thermo's OO machinery.
    Run in a subprocess so it is independent of test ordering.
    """
    code = (
        "import sys; import aeroprops; "
        "assert 'rdkit' not in sys.modules, 'rdkit imported at package import'; "
        "assert 'thermo' not in sys.modules, 'thermo imported at package import'"
    )
    subprocess.run([sys.executable, "-c", code], check=True)
