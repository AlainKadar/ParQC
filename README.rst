=====
ParQC
=====

Overview
========
A parallelized phason strain algorithm for constructing the hyperlattice matrix from molecular simulations. The user must write the position and neighbor list, which can be done with the provided jupyter notebook. The C++ file will write ``H*.txt``. This project was initially conceived for fulfilling a requirement for `EECS 587 <https://web.eecs.umich.edu/~qstout/587_Overview.pdf>`__ at the University of Michigan

Use
===
To use, clone and compile from the `GitHub repository
<https://github.com/AlainKadar/ParQC>`__. You must have the OpenMP library and a compatible OpenMP compiler

.. code:: bash

   git clone https://github.com/AlainKadar/ParQC
   g++ -std=c++11 -o lift lift.cpp -fopenmp

To use the GPU implementation, you must have the CUDA compiler, ``nvcc``:

.. code:: bash

   nvcc -std=c++11 -o lift lift.cu  -Xcompiler -fopenmp

You can check the executables work by testing whether running the ``Verifier.ipynb`` notebook reconstructs the gsd file from the hyperlattice matrix computed by ``lift``. This notebook requires the `freud <https://anaconda.org/conda-forge/freud>`__ and `gsd <https://anaconda.org/conda-forge/gsd>`__ python libraries.
The code can be run on new gsd files by first writing the neighbor and position lists with the notebook and then running the code from either the terminal or a bash script:

.. code:: bash

   GSD_NAME='new'
   Q_THREADS=1
   NL_THREADS=2
   V_THREADS=2
   ORIGIN=0
   SCALE=1                                                                
                                                                                 
   ./lift $GSD_NAME $Q_THREADS $NL_THREADS $V_THREADS $ORIGIN $SCALE

