Automated re-parametrization of global non-covalent parameters in SFAM
======================================================================

Introduction
------------

With the script "noncovalent_opt.py" from this directory it is possible to parametrize the global non-covalent
parameters of SFAM, i.e., those for dispersion (a1, s8, and a2), the one for repulsion (beta), and the atomic
charges scaling factor. These are the parameters that can also be modified in the SFAM parameter file in the
section "! non-covalent".

Why is this necessary?
----------------------

Currently, we have only done this parametrization for the reference method PBE-D3(BJ)/def2-SVP for the programs ORCA
(with correct CM5 charges) and Turbomole (pseudo-CM5 charges with Löwdin origin). These parameters will be written
to the parameter file as well after a parametrization with a different reference method or with the Sparrow or
xTB programs. For more accurate results, we suggest to re-parametrize the global non-covalent parameters with the
script "noncovalent_opt.py" for your reference method.

How can I do this?
------------------

We provide reference data for a test set of 22 molecular trajectories as presented in the Supporting Information
of the original SFAM publication. This reference data can be found in the directory
"data_for_noncovalent_parametrization". In subdirectories, the data for each of the systems are stored.
In each of these system directories, one can add a parameter and a connectivity file (for a chosen method).
The default names for these files are "Parameters.dat" and "Connectivity.dat", however, this can be adapted in the
Python script. Hence, one must first parametrize each one (or a subset) of these test systems. If this does not work
for one of the systems, one can simply delete this system's directory and perform the parameter optimization without
considering it. The reference data was calculated with PBE0-D3(BJ)/def2-SVP and is already available in the files
"orca_trj.xyzact.dat" in each of the system directories.

Once the parameters and connectivity files are stored in their respective directories, one can run the Python script
and obtain the optimized parameters. Note that the script uses SCINE Python libraries, hence, you must make sure
that the Python bindings of Swoose are installed and that the PYTHONPATH is set correctly (see Swoose manual).
The script also creates a folder with the name "figures" (name can be adapted in the settings at the top of
the Python script), in which helpful figures are stored that document the results of the model fit to the trajectories.

Where can I get help?
---------------------

If any questions arise, do not hesitate to contact the developers of Swoose via scine@phys.chem.ethz.ch.
