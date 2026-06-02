Quick start
===========

PyDecoys is available at PyPI and can be easily set up.

Dependencies
------------

- Requires Python 3.12 or superior.
- CLI and IO functionalities require Biopython. This should be automatically
  installed.

Installation
------------

If you only care about the CLI, PyDecoys can be easily installed or uninstalled
via `pipx`:

.. code-block:: sh

    # Install current version
    pipx install git+https://github.com/CPBMF/pydecoys@v0.2.0

    # Uninstall
    pipx uninstall pydecoys

`pipx` should automatically install PyDecoys and make it globally available on
PATH. You can confirm it's available by running ``pydecoys -v``. In case it
isn't, run:

.. code-block:: sh

    pipx ensurepath

This will ensure all `pipx` scripts are available. The `pipx` app itself can be
installed via `pip`:

.. code-block:: sh

    pip install --user pipx

In case `pipx` didn't automatically install [Biopython], run:

.. code-block:: sh

    pipx inject pydecoys biopython

For API usage, you can use `pip` or other package managers:

.. code-block:: sh

    # Install
    pip install git+https://github.com/CPBMF/pydecoys@v0.2.0

    # Uninstall
    pip uninstall pydecoys

Installing without Biopython
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To install PyDecoys without `Biopython`, use the ``--no-deps`` flag:

.. code-block:: sh

    pip install git+https://github.com/CPBMF/pydecoys@v0.2.0 --no-deps

If you change your mind later, simply run ``pip install biopython``.
`Biopython`'s dependency group is `biopython`.

Note that IO functions and the CLI app **aren't available without Biopython**.
