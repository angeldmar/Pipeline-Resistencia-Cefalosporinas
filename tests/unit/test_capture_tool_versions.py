"""Pruebas unitarias de capture_tool_versions.py.

Encontrado corriendo el pipeline real (no en pruebas sinteticas): el script
original llamaba cada herramienta por nombre suelto, pero corria dentro del
ambiente generico python.yaml -- que nunca tiene instaladas fastp, spades,
quast, etc. (cada una vive en su propio ambiente Conda por regla). Por eso
"data/metadata/tool_versions.tsv" mostraba "not installed" para todo, aunque
las herramientas si estuvieran instaladas en sus propios ambientes.
"""

from capture_tool_versions import (
    CHECKM_VERSION_PATTERN,
    MissingCondaEnvironmentError,
    get_tool_version,
)


def test_checkm_version_pattern_extracts_from_help_banner():
    # checkm no soporta --version (no es un subcomando valido); la version
    # real solo aparece en el encabezado de "checkm -h".
    help_output = "\n                ...::: CheckM v1.0.18 :::...\n\n  Lineage-specific marker set:\n"
    match = CHECKM_VERSION_PATTERN.search(help_output)
    assert match.group(0) == "CheckM v1.0.18"


def test_get_tool_version_returns_not_installed_when_env_missing(monkeypatch):
    import capture_tool_versions

    def raise_missing_env(_env_yaml_relative_path):
        raise MissingCondaEnvironmentError("no existe")

    monkeypatch.setattr(capture_tool_versions, "resolve_conda_env_bin", raise_missing_env)

    assert get_tool_version("prokka", ["prokka", "--version"], "workflow/envs/prokka.yaml") == "not installed"
