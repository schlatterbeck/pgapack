.. _chp:quick-start:

Quick Start
-----------

If you wish to get started by just typing a few lines and running an
example, this section is for you. We trust you know how to check out the
latest `version from github`_. Once you have a directory with the
checked-out version you can build with:

#. ``cd pgapack``

#. ``make``

The Makefile will auto-detect if you have an MPI-Implementation
installed and will build a parallel version. If no MPI is detected, the
serial version will be built. If the version is not correctly
auto-detected or you want to force a certain MPI backend, refer to the
build documentation in ``README.rst``.

Chapter :ref:`chp:examples`  (example problems) and Sections
:ref:`sec:big-picture` and
:ref:`sec:evaluation`
should be read next.

.. _`version from github`: https://github.com/schlatterbeck/pgapack
