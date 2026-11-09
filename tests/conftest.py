# SPDX-License-Identifier: MIT
"""Shared pytest fixtures and helpers for the aeroprops test suite."""


def pytest_configure(config):
    config.addinivalue_line(
        "markers", "rdkit: test requires RDKit (deselect with -m 'not rdkit')"
    )
