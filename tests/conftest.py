import importlib
from pathlib import Path

import pytest

import pydecoys
from pydecoys.strategies import RAND


@pytest.fixture(autouse=True)
def reseed():
    RAND.seed(10)


@pytest.fixture
def root():
    return Path(__file__).parent


@pytest.fixture
def without_bio(missing_modules):
    warn = "Module 'Biopython' not found: .+"
    with missing_modules('Bio'), pytest.warns(UserWarning, match=warn):
        importlib.reload(pydecoys.core)  # type: ignore
        yield
    importlib.reload(pydecoys.core)  # type: ignore
