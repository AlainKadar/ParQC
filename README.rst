=====
ParQC
=====

Overview
========
A parallelized phason strain algorithm for constructing the hyperlattice
matrix from molecular simulations. The user must write the position and
neighbour list, which can be done with the provided jupyter notebook. The C++
file will write H*.txt

Use
===
To use, clone and compile from the `GitHub repository
<https://github.com/AlainKadar/ParQC>`__. You must have the OpenMP library and
a compatible OpenMP compiler

.. code:: bash

   git clone https://github.com/AlainKadar/ParQC
   g++ -std=c++11 -O0 -fopenmp -o lift lift.cpp 

You can check the executable works by testing whether running the ``Verifier.ipynb``
notebook reconstructs the gsd file from the hyperlattice matrix computed
by ``lift``. This notebook requires the freud and gsd python libraries.

The code can be run on new gsd files by first writing the neighbor and
position lists with the notebook and then running the code from either the
terminal or a bash script.

.. code:: bash

GSD_NAME='new'                                                               
Q_THREADS=1                                                                    
NL_THREADS=2                                                                   
V_THREADS=2                                                                    
ORIGIN=0                                                                       
                                                                               
./lift $GSD_NAME $Q_THREADS $NL_THREADS $V_THREADS $ORIGIN
# ParQC
