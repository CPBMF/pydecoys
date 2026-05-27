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

If you only care about the CLI, PyDecoys can be easily installed, upgraded or
uninstalled via `pipx`:

.. code-block:: sh

    # Install current version
    pipx install pydecoys

    # Upgrade to latest version
    pipx upgrade pydecoys

    # Uninstall
    pipx uninstall pydecoys

`pipx` should automatically install PyDecoys and make it globally available on
PATH. You can confirm it's available by running ``pydecoys -v``. You should see
PyDecoy's version along with its license. In case you don't, run:

.. code-block:: sh

    pipx ensurepath

This will ensure all `pipx` scripts are available. The `pipx` app itself can be
installed via `pip`:

.. code-block:: sh

    pip install --user pipx

For API usage, you can use `pip` or other package managers:

.. code-block:: sh

    # Install
    pip install pydecoys

    # Upgrade
    pip install --upgrade pydecoys

    # Uninstall
    pip uninstall pydecoys

Without Biopython
^^^^^^^^^^^^^^^^^
PyDecoys can be used without Biopython:

.. code-block:: sh

    pip install pydecoys --no-deps

If you change your mind later, simply run ``pip install biopython``.
Biopython's dependency group is ``biopython``.

Note that IO functions and the CLI app **aren't available without Biopython**.
