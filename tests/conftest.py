import pytest

from pydecoys.strategies import RAND


@pytest.fixture(autouse=True)
def reseed():
    RAND.seed(10)
