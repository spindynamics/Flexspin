# Flexspin
Flexspin is an open source software simulator, time resolving mutliple Landau-Lifchitz-Gilbert equation for highly configurable macrospin simualtion. Developped at the French Alternative Energies and Atomic Energy Commission, Flexspin is a macrospin simulation tool for exploring the dynamics of magnetic tunnel junction (MTJ) based systems. Initially developed for internal use, it is now being opened to serve as a general-purpose tool for the community, some features might be missing, don't hesitate to reach to the contributor or implement them yourself via forking.

# Capabilities
Macrospin simulation of the static and dynamic magnetic properties of magnetic tunnel junctions, including description of any number of macrospin layer, in any number of separeted junction. Thermally activated switching, spin-transfer torque driven dynamics via DC or AC bias (Current, Voltage, Field).

### General simulation framework
 A parallelized sweep usable via MPI for quicker simulation
 A continuous and hysteretic sweep for each MPI node
 End state or time step printing of physical parameters
### Simulation framework for each layer
 Stochastic Landau-Lifshitz-Gilbert-Slonczewski equation
 Macrospin description
 Free or fixed layer
### Magnetic anisotropy
 Bulk uniaxial anisotropy, first and second order
 Interface anisotropy, first and second order
 Shape anisotropy from layer dimensions
 Voltage-controlled magnetic anisotropy (VCMA)
### Magnetic interactions
 Dipolar coupling between layers
 Interlayer exchange coupling (RKKY)
 Exchange bias
 Zeeman energy
### Spin torques and transport
 Spin-transfer torque, damping-like and field-like
 Bias-dependent spin-transfer torque
 Angle-dependent resistance through tunnel magnetoresistance (TMR)
 Bias dependence of TMR
### Thermal effects
 Stochastic thermal field
 Temperature dependence of magnetic parameters (Callen-Callen scaling)
 Temperature dependence of TMR
 Temperature dependence of spin-transfer torque
 Joule heating with dynamic junction temperature
### Excitations
 DC and AC current
 DC and AC voltage
 DC and AC magnetic field of arbitrary direction
 Trapezoidal pulses (delay, rise, plateau, fall)
 Bias applied to the whole system or to each junction individually
### Code features
 Multi-junction systems with a common bias
 Heun, fourth-order Runge-Kutta and symplectic solvers
 Optional adaptive time step
 JSON input file
 User-defined output columns

## Contributors
* [Arthur Courberand](https://github.com/ArthurCourb)

 
Flexspin is governed by the CeCILL license under French law and abiding by the rules of distribution of free software. You can use, modify and/ or redistribute the software under the terms of the CeCILL license as circulated by CEA, CNRS and INRIA.