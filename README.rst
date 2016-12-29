========
Overview
========

.. start-badges

.. list-table::
    :stub-columns: 1

    * - docs
      - |docs|
    * - tests
      - |
        |
    * - package
      - |version| |downloads| |wheel| |supported-versions| |supported-implementations|

.. |docs| image:: https://readthedocs.org/projects/pk_milk/badge/?style=flat
    :target: https://readthedocs.org/projects/pk_milk
    :alt: Documentation Status

.. |version| image:: https://img.shields.io/pypi/v/pk_milk.svg?style=flat
    :alt: PyPI Package latest release
    :target: https://pypi.python.org/pypi/pk_milk

.. |downloads| image:: https://img.shields.io/pypi/dm/pk_milk.svg?style=flat
    :alt: PyPI Package monthly downloads
    :target: https://pypi.python.org/pypi/pk_milk

.. |wheel| image:: https://img.shields.io/pypi/wheel/pk_milk.svg?style=flat
    :alt: PyPI Wheel
    :target: https://pypi.python.org/pypi/pk_milk

.. |supported-versions| image:: https://img.shields.io/pypi/pyversions/pk_milk.svg?style=flat
    :alt: Supported versions
    :target: https://pypi.python.org/pypi/pk_milk

.. |supported-implementations| image:: https://img.shields.io/pypi/implementation/pk_milk.svg?style=flat
    :alt: Supported implementations
    :target: https://pypi.python.org/pypi/pk_milk


.. end-badges

From Mother To Child: Generic, Generational Pharmacokinetic 1-Box Model of the Transport of Pollutants in Breastmilk

* Free software: BSD license

Installation
============

::

    pip install pk_milk

Documentation
=============

https://pk_milk.readthedocs.io/

Development
===========

To run the all tests run::

    tox

Note, to combine the coverage data from all the tox environments run:

.. list-table::
    :widths: 10 90
    :stub-columns: 1

    - - Windows
      - ::

            set PYTEST_ADDOPTS=--cov-append
            tox

    - - Other
      - ::

            PYTEST_ADDOPTS=--cov-append tox
