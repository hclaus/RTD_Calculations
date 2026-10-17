"""Loads and validates app/metadata/variables.yaml -- the single source of
truth for variable text, consumed by both the input tooltips and the Help
page glossary via GET /api/variables."""

from functools import lru_cache
from pathlib import Path

import yaml
from pydantic import BaseModel

_VARIABLES_PATH = Path(__file__).parent / "variables.yaml"


class VariableMeta(BaseModel):
    id: str
    symbol: str
    label: str
    units: str
    category: str
    default: float | None = None
    tooltip: str
    description: str
    methodology_section_anchor: str


@lru_cache(maxsize=1)
def load_variables() -> list[VariableMeta]:
    raw = yaml.safe_load(_VARIABLES_PATH.read_text(encoding="utf-8"))
    variables = [VariableMeta(**entry) for entry in raw]

    ids = [v.id for v in variables]
    duplicates = {i for i in ids if ids.count(i) > 1}
    if duplicates:
        raise ValueError(f"Duplicate variable id(s) in variables.yaml: {duplicates}")

    return variables
