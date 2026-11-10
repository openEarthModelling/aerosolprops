# SPDX-License-Identifier: MIT
"""Shared UManSysProp oracle harness (behavioural oracle only).

UManSysProp is run live and unmodified as the reference behaviour for
aerosolprops' re-implementations; none of its SMARTS, data files or code
are copied into aerosolprops (MIT discipline). This module centralises the
RDKit-based ``openbabel.pybel`` / ``pkg_resources`` shims from the
UManSysProp tutorial so every oracle test module (SIMPOL.1, Nannoolal)
shares one harness.

The one OpenBabel-dialect SMARTS in the UMS catalogue that RDKit
rejects is rewritten inside the shim - a compatibility fix of the
oracle harness, not part of any aerosolprops catalogue.
"""

import os
import sys
import types

#: Local UManSysProp checkout (oracle only; never a runtime dependency).
DEFAULT_USP_REPO = (
    "/home/zhangfan/Project/20260529_ChemSim/jlbox_study/"
    "umansysprop_src/UManSysProp_public-master"
)

USP_REPO = os.environ.get("AEROSOLPROPS_UMS_REPO", DEFAULT_USP_REPO)


def available():
    """True when the UManSysProp source tree is present."""
    return os.path.isdir(USP_REPO)


def install_shims():
    """Install the tutorial shims and import the UMS method modules.

    Returns a namespace with ``groups``, ``vapour_pressures`` and
    ``boiling_points``. Re-importing resets the cached UMS modules so
    repeated calls start from a clean state.
    """
    if USP_REPO in sys.path:
        sys.path.remove(USP_REPO)
    for name in (
        "umansysprop",
        "umansysprop.groups",
        "umansysprop.vapour_pressures",
        "umansysprop.boiling_points",
    ):
        sys.modules.pop(name, None)

    from rdkit import Chem as _Chem
    from rdkit import RDLogger as _RDLogger

    # The first parse attempt of the OpenBabel-dialect pattern below
    # logs a (harmless, recovered-from) RDKit error; keep test output
    # clean.
    _RDLogger.DisableLog("rdApp.error")

    rewrite = {
        # UMS SIMPOL group-2 pattern uses an OpenBabel dialect RDKit
        # rejects; drop the second alternative (UMS overwrites
        # result["2"] for amides anyway - see the SIMPOL whitelist).
        r"([#6][!$([NX3][CX3]=O),$([CX3](=O)[NX3])])": r"([#6][!$([NX3][CX3]=O)])",
    }

    class _Smarts(object):
        def __init__(self, smartsmi):
            self.source = smartsmi
            s = rewrite.get(smartsmi, smartsmi)
            patt = _Chem.MolFromSmarts(s)
            if patt is None and s.startswith("(") and s.endswith(")"):
                patt = _Chem.MolFromSmarts(s[1:-1])
            if patt is None:
                raise IOError("bad SMARTS %s" % smartsmi)
            self._patt = patt

        def findall(self, mol):
            # uniquify=True reproduces OpenBabel findall's deduplicated
            # match list, the semantics the method notes pin down.
            return mol._rdmol.GetSubstructMatches(self._patt, uniquify=True)

    class _Molecule(object):
        def __init__(self, smi):
            self._rdmol = _Chem.MolFromSmiles(smi)
            self._smiles = smi

        @property
        def molwt(self):
            from rdkit.Chem import Descriptors

            return Descriptors.MolWt(self._rdmol)

        def __str__(self):
            return self._smiles

    def _readstring(fmt, smi):
        return _Molecule(smi.decode() if isinstance(smi, bytes) else smi)

    ob = types.ModuleType("openbabel")
    pybel = types.ModuleType("openbabel.pybel")
    pybel.Smarts = _Smarts
    pybel.readstring = _readstring
    ob.pybel = pybel
    sys.modules["openbabel"] = ob
    sys.modules["openbabel.pybel"] = pybel

    pkg_resources = types.ModuleType("pkg_resources")

    def _resource_stream(package, filename):
        return open(os.path.join(USP_REPO, "umansysprop", "data", filename), "rb")

    pkg_resources.resource_stream = _resource_stream
    sys.modules["pkg_resources"] = pkg_resources

    sys.path.insert(0, USP_REPO)
    from umansysprop import boiling_points, groups, vapour_pressures

    return types.SimpleNamespace(
        groups=groups,
        vapour_pressures=vapour_pressures,
        boiling_points=boiling_points,
    )


def readstring(ums, smi):
    """Parse `smi` through the installed shim (what UMS code sees)."""
    import openbabel.pybel as _pybel

    return _pybel.readstring("smi", smi)


def example_output_smiles():
    """The 26 SMILES of the UMS Example_output.txt walkthrough."""
    path = os.path.join(USP_REPO, "umansysprop", "examples", "Example_output.txt")
    with open(path) as handle:
        lines = handle.read().splitlines()
    return [line.split("\t")[0] for line in lines[1:] if line.strip()]
