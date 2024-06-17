Changelog
=========

Release 2.0.0
-------------
- enable electrostatic embedding QM/MM with Turbomole and xtb as QM calculators
- interface to any SCINE calculator for QM
- automated construction of QM regions around multiple QM region center atoms
- identify critical failed calculations for a parametrization as debug information
- semi-automated processing of molecular input structures (ASAP)
- Python bindings for SFAM parametrization and QM region selection objects with release of the global interpreter lock for threading long-running operations
- Update address in license

Release 1.0.0
-------------
- parametrization of SFAM molecular mechanics model (large and small systems)
- GAFF molecular mechanics model
- interface through SCINE to xTB and Sparrow
- automated set-up of QM/MM calculations (including automated QM region selection)
- single-point calculations with the molecular mechanics models and with QM/MM models
- Hessian matrix available with MM models
- structure optimizations and molecular dynamics simulations with MM models and with QM/MM models
- library functions to build molecular machine learning models from quantum chemical reference data
