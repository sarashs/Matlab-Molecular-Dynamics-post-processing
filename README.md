# Matlab-Molecular-Dynamics-post-processing
A GUI to perform simple post MD analysis on Molecular Dynamics simulations
This GUI is easy to use. You only need to structurecreator.m from inside Matlab and then the GUI will pop up.
You can then Load MD data that is previously saved in an "extended xyz format" from Ovito software. Note that the extended xyz format should be like:
AtomID AtomType Charge x y z
The GUI can calculate RDF, Neighbors list for all of the atoms
Angles, bond lengths, coordination numbers, percentage of each atom, and charge distribution along various dimensions
each of these functionalities are via a different m file that you can change for your needs.
