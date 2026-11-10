# SPDX-License-Identifier: MIT
"""Shared pytest fixtures for the aeroprops test suite."""

import pytest

import ums_oracle


@pytest.fixture(scope="session")
def ums():
    """Live UManSysProp modules behind the RDKit tutorial shims."""
    return ums_oracle.install_shims()
