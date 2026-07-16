"""Configuracion compartida de pytest para toda la suite de pruebas.

Los scripts de workflow/scripts/ son ejecutables independientes (pensados
para invocarse por linea de comandos o desde una regla de Snakemake), no un
paquete Python instalado. Para poder probar sus funciones directamente sin
duplicar logica en los tests, se agrega esa carpeta a sys.path aqui, una
sola vez, antes de que pytest recolecte ningun test.
"""

import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parent.parent
SCRIPTS_DIR = REPO_ROOT / "workflow" / "scripts"

sys.path.insert(0, str(SCRIPTS_DIR))


@pytest.fixture
def repo_root() -> Path:
    """Ruta absoluta a la raiz del repositorio, para que los tests puedan
    referenciar archivos de config/ sin depender del directorio desde el
    que se invoco pytest."""
    return REPO_ROOT
