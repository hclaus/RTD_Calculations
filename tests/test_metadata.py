import re
from pathlib import Path

import pytest

from app.metadata.loader import load_variables

STATIC_DIR = Path(__file__).parent.parent / "static"


def test_variables_yaml_parses_and_validates():
    variables = load_variables()
    assert len(variables) > 0
    for v in variables:
        assert v.id
        assert v.category in {
            "shared", "microbiological", "solver_settings", "derived", "chemical", "ozone",
        }


def test_no_duplicate_ids():
    variables = load_variables()
    ids = [v.id for v in variables]
    assert len(ids) == len(set(ids))


def test_every_methodology_anchor_exists_in_help_html():
    """Every methodology_section_anchor referenced by a variable must
    correspond to a real id="..." in help.html, or tooltip 'Learn more'
    links will 404/no-op silently."""
    help_html = (STATIC_DIR / "help.html").read_text(encoding="utf-8")
    existing_ids = set(re.findall(r'id="([^"]+)"', help_html))

    variables = load_variables()
    anchors = {v.methodology_section_anchor for v in variables}
    missing = anchors - existing_ids
    assert not missing, f"Anchors referenced by variables.yaml but missing from help.html: {missing}"


def test_required_input_variables_present():
    variables = load_variables()
    ids = {v.id for v in variables}
    required = {"L", "a", "Q", "mu", "rho", "D", "k", "Eavg", "measured_C_C0",
                "k_min", "k_max", "eavg_min", "eavg_max"}
    missing = required - ids
    assert not missing, f"Missing required variable ids: {missing}"
