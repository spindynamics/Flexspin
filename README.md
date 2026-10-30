# Flexspin

Flexspin is an MPI-parallel C++11 simulator for stochastic Landau-Lifshitz-Gilbert dynamics on macrospin spin-transfer-torque magnetic multilayer junctions.

![Flexspin](https://mark.sylphx.com/api/v1/mark/hero.svg?text=Flexspin&desc=Fast.+Stochastic.+Parallel.+High-performance+STT+macrospin+simulator&theme=dark)

## Capabilities

Flexspin simulates the static and dynamic magnetic properties of magnetic tunnel junction (MTJ) systems with configurable macrospin layers.

### General Simulation Framework
- Parallelized sweeps via MPI for faster computation
- Continuous and hysteresis sweeps per MPI node
- End state or time step output of physical parameters

### Simulation Features per Layer
- Stochastic Landau-Lifshitz-Gilbert-Slonczewski equation
- Macrospin description with free or fixed layers
- Multiple macrospin layers in any number of junctions

### Magnetic Anisotropy
- Bulk uniaxial anisotropy (first and second order)
- Interface anisotropy (first and second order)
- Shape anisotropy from layer dimensions
- Voltage-controlled magnetic anisotropy (VCMA)

### Magnetic Interactions
- Dipolar coupling between layers
- Interlayer exchange coupling (RKKY)
- Exchange bias
- Zeeman energy

### Spin Torques and Transport
- Spin-transfer torque with damping-like and field-like components
- Bias-dependent spin-transfer torque
- Angle-dependent resistance via tunnel magnetoresistance (TMR)
- Bias dependence of TMR

### Thermal Effects
- Stochastic thermal field (Langevin dynamics)
- Temperature dependence of magnetic parameters
- Temperature dependence of TMR
- Temperature dependence of spin-transfer torque
- Joule heating with dynamic junction temperature

### Excitations
- DC and AC current and voltage
- DC and AC magnetic field of arbitrary direction
- Trapezoidal pulses (delay, rise, plateau, fall)
- Bias applied to whole system or each junction individually

### Code Features
- Multi-junction systems with common bias
- Multiple solvers: Heun, fourth-order Runge-Kutta, symplectic
- Optional adaptive time step
- JSON input file
- User-defined output columns

## Installation

### Prerequisites
- C compiler (gcc, clang, or similar)
- CMake >= 3.26
- MPI library
- json-c development headers

### Build
```bash
cmake -B build .
cmake --build build
```

Configure with:
```bash
cmake -B build . \
  -DFLEXSPIN_WARNINGS=ON \
  -DCMAKE_BUILD_TYPE=Debug
```

## Running

```bash
mpirun -np 4 ./build/Flexspin input.json
```

## Usage

- Provide input via `input.json` (specify filename as command argument if different)
- input.json is required but not in the repo
- JSON keys include units: `"Ms (A/m)"`, `"Layer Size (nm)"`, `"pulse (ns)"`
- JSON schema documented in `src/read_input.cpp`
- See `AGENTS.md` for detailed usage, architecture, and conventions

## Architecture

Flexspin uses a modular design with three core headers and 16 source files:

- `objects.hpp` - All data structures and declarations
- `printing.hpp` - Output formatting (includes objects.hpp)
- `Constants.hpp` - Physics constants and M_PI fallback

The main flow: `main` (MPI init, JSON parse) → `Simulation::range` (sweep loops) → `Simulation::timeStep` (pulse → LLG integration → Joule heating).

### Available Solvers
- Heun's method
- Fourth-order Runge-Kutta (RK4)
- Symplectic integrator

## Contributors

- [Arthur Courberand](https://github.com/ArthurCourb)
- [Pascal Thibaudeau](https://github.com/pthibaud)
- [Liliana Prejbeanu](https://github.com/liliprejbeanu)

## License

Flexspin is governed by the CeCILL license under French law. See `Licence_CeCILL_V2.1-en.txt` for details. You can use, modify, and redistribute the software under the CeCILL license as circulated by CEA, CNRS, and INRIA.