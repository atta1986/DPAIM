 # Dispersion from polarizability of the atom in molecule: DPAIM 
 ## Usage:
/path/to/the/dpaim.exe file.dat

Important: The files param_damp.dat (containing 4 damping parameters) and param_ckl.dat (containing 3 c_kl parameters) must be present in the current working directory

The file.dat has the following structure. The XYZ coordinates are in Bohr units.

 Number of atoms in monomer A
 X(atom 1) Y(atom 1) Z(atom 1) charge(atom 1) V3(atom 1) V5(atom 1) atomic_mass(atom 1) .... \
 X(atom 2) Y(atom 2) Z(atom 2) charge(atom 2) V3(atom 2) V5(atom 2) atomic_mass(atom 2) ....  \ monomer
 X(atom 3) Y(atom 3) Z(atom 3) charge(atom 3) V3(atom 3) V5(atom 3) atomic_mass(atom 3) ....  /    A
   ...       ...       ...          ...       .... /
 Number of atoms in monomer B
 X(atom 1) Y(atom 1) Z(atom 1) charge(atom 1) V3(atom 1) V5(atom 1) atomic_mass(atom 1) .... \
 X(atom 2) Y(atom 2) Z(atom 2) charge(atom 2) V3(atom 2) V5(atom 2) atomic_mass(atom 2) ....  \ monomer
 X(atom 3) Y(atom 3) Z(atom 3) charge(atom 3) V3(atom 3) V5(atom 3) atomic_mass(atom 3) ....  /    B
   ...       ...       ...          ...       .... /

***

The following are the files needed to compile the DPAIM code:

DPAIM.f90 -> This is the main DPAIM code file
polar_interp.f90 -> This a module file used when the polarizabilities are interpolated based on the effective atomic volume
dipole_polarizabilities_f90.f90 -> Contain the dipole polarizabilities of the reference molecules
quad_polarizabilities_f90.f90 -> Contain the quadrupole polarizabilities of the reference molecules
volume_V3 -> Contain the effective volume of the reference molecules calculted using r^3
volume_V5 -> Contain the effective volume of the reference molecules calculted using r^5
dpaim.exe -> is the DPAIM executable


***

The following command is used to compile the code:
gfortran -O2 -cpp  polar_interp.f90 quad_polarizabilities_f90.f90 dipole_polarizabilities_f90.f90 DPAIM.f90 -o dpaim.exe

***

The two main subroutine in DPAIM.f90 are:
Subroutine fitted_dispersion() -> calculates the dispersion energy
Subroutine AIM_polarizabilities() -> calculates the AIM polarizabilities

***
