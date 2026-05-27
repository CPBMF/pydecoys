.. pydecoys documentation master file, created by
   sphinx-quickstart on Sun May 17 15:48:49 2026.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

PyDecoys documentation
======================

Welcome to the GitHub page of PyDecoys:
Proteomics decoy utilities for Python!

About
-----

PyDecoys is a bioinformatics Python package and application for generating
decoy proteins from target proteins.

Decoy generation is an essential step of Proteomics workflows: decoy proteins
are included in the search space as necessarily false hits; assuming decoy
matches are as likely as target false-positive matches, decoy/target match
ratios allow for the estimation of false-positive error ratios.

PyDecoys aims to facilitate this decoy generation step in a way that's highly
customizable and easily integratable to existing Proteomics workflows.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   guide/index
   cli
   api/index
   available
