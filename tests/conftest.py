from pathlib import Path

import pytest

from pydecoys.strategies import RAND


@pytest.fixture(autouse=True)
def reseed():
    RAND.seed(10)


@pytest.fixture
def root():
    return Path(__file__).parent
