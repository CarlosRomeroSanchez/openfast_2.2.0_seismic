OpenFAST
========

|travisci| |nbsp| |rtfd|

.. |travisci| image:: https://travis-ci.org/OpenFAST/openfast.svg?branch=dev
   :target: https://travis-ci.org/OpenFAST/openfast
   :alt: Build Status
.. |rtfd| image:: https://readthedocs.org/projects/openfast/badge/?version=dev
   :target: https://openfast.readthedocs.io/en/dev
   :alt: Documentation Status
.. |nbsp| unicode:: 0xA0
   :trim:

This repository contains a modified version of OpenFAST 2.2.0 in which the Subdyn module has been modified to include the two following features:

* Seismic foundation input motions at the base of the support structure. 
* A dynamic soil-structure interaction model that employs a Simplified Lumped Parameter Model [1] to represent the lateral, vertical, rocking and torsional impedances of the foundation.

**OpenFAST is under active development**.

 OpenFAST v2.2.0
-------------------------
In order to be able to use these features, the following changes have been introduced into the OpenFAST input files:

1. File con extensión tal incluye una nueva línea con la variable tal para indicar si hay sismo o no.
2. Ese file incluye una nueva línea tal con la dirección del fichero donde lee el sismo
3. Ese fichero tiene formato tal
4. La definición de los SLPM se indica de tal manera ..

Uniform sismic input motions are assumed (i.e., if there exists more than una support, all of them experience the same foundation input motion). On the other hand, the parameters that define the SLPM must be obtained by fitting the target impedance functions. 

The original version of OpenFAST 2.2.0 can be found here (incluir link).

[1] Carbonari S., Morici M., Dezi F., Leoni G., A lumped parameter model for timedomain inertial soil-structure interaction analysis of structures on pile foundations , Earthquake Engineering & Structural Dynamics, Vol. 47, 2147-2171.(2018)

Referenciar el paper del congreso como documentación adicional para entender esto.
