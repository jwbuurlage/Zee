.. Zee documentation master file, created by
   sphinx-quickstart on Mon Jul  6 14:59:52 2015.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Zee; a partitioning framework for sparse matrices
=================================================

Zee is a header-only C++ library which serves as a framework for sparse matrix (SpM) operations. In particular it has been developed to ease the implementation of SpM partitioners to aid the execution speed of parallel algorithms such as the distributed SpMV algorithm and iterative solvers.

Using Zee
---------

.. toctree::
   :maxdepth: 2

   installing
   getting_started
   support

Components
----------

The library is made up of a number of components, which are introduced in detail in the pages below:

.. toctree::
   :maxdepth: 2

   distributed_matrices
   storage
   linear_algebra
   partitioning


API Overview
============

.. toctree::
   :maxdepth: 2

   DMatrixBase
   DSparseMatrix
   DSparseMatrixImage
   DSparseStorage
   Triplet

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
