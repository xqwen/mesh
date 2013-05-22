MeSH - Meta analysis with Subgroup Heterogeneity
==========================================

This directory contains MeSH, a package implementing Bayesian genetic association analysis in heterogeneous subgroups.

MeSH is free software, you can redistribute it and/or modify it under
the terms of the GNU General Public License.

The GNU General Public License does not permit this software to be
redistributed in proprietary programs.

This software is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.



Compilation and Installation
=============================================

The compilation of MeSH from the source code requires the GNU c++ compiler (g++), GNU make and GNU scientific library (GSL) pre-installed in the compiling machine. 

To compile the executable, run the following commands from the current directory:

     cd src
     make

Upon successful compilation, a binary executable "mesh" should be produced.   


Documentation and Example Data
=============================================

A detailed documentation "mesh.pdf" can be found in doc/ directory, we also include a simulated sample data set in the sample_data/ directory.  


Citation
=============================================

Please cite the following for usage of MeSH:

Wen, X and Stephens, M. 2012 "Bayesian Methods for Genetic Association Analysis with Heterogeneous Subgroups: from Meta-Analyses to Gene-Environment Interactions". arXiv pre-print: 1111.1210 