# AGENTS.md

## What this is

Flexspin: an MPI-parallel C++11 simulator for the stochastic Landau-Lifshitz-Gilbert (LLG)
equation on macrospin spin-transfer-torque magnetic multilayer junctions — demag, dipolar,
RKKY, STT, Joule self-heating, Langevin thermal field. Single binary, not a library.
CeCILL v2.1 (`Licence_CeCILL_V2.1-en.txt` — leave it alone).

## Build / run

CMake >= 3.26, out of source. `Unix Makefiles` is the default generator on Linux:

```
cmake -S . -B build
cmake --build build -j
```

```
mpirun -np 4 ./build/Flexspin input.json        # or: ./build/Flexspin input.json
```

Add `-G "Unix Makefiles"` where the platform default differs (e.g. Windows).

Useful options: `-DFLEXSPIN_WARNINGS=ON` (`-Wall -Wextra`, noisy today), `-DFLEXSPIN_INSTALL=OFF`,
`-DCMAKE_BUILD_TYPE=Debug` (defaults to `Release`, i.e. `-O3`).

### json-c

Every file includes `objects.hpp` -> `<json-c/json.h>`, so the json-c *development*
headers are needed — a runtime `.so` alone is not enough. `CMakeLists.txt` resolves this in
three steps and prints which one it picked:

1. system CMake package (`json-c::json-c`; needs json-c >= 0.15 installed with CMake files)
2. system headers + library via `JSONC_INCLUDE_DIR` / `JSONC_LIBRARY` (old distro packages,
   e.g. RHEL 8's `json-c-devel`, which ship headers and `json-c.pc` but no config file)
3. **vendored**: compile json-c from `https://github.com/json-c/json-c.git` at
   `FLEXSPIN_JSONC_GIT_TAG` into the build tree and link it statically

So a bare `cmake -S . -B build` works on a machine with no json-c package at all. Step 3 needs
git and network access on the *first* configure only. To pin a different version, build offline,
or force a route:

```
cmake -S . -B build -DFLEXSPIN_JSONC_GIT_TAG=json-c-0.18-20240915   # different json-c
cmake -S . -B build -DFETCHCONTENT_SOURCE_DIR_JSON-C=/path/to/json-c # offline / local checkout
cmake -S . -B build -DFLEXSPIN_JSONC=vendored   # ignore any system json-c
cmake -S . -B build -DFLEXSPIN_JSONC=system     # fail instead of downloading
```

`input.json` is required at runtime (default name when no `argv[1]`) and is **not** in the
repo — `src/read_input.cpp` is the only schema documentation, and JSON key names carry the
units (`"Ms (A/m)"`, `"Layer Size (nm)"`, `"pulse (ns)"`); nm/ns -> SI conversions happen at
parse time.

## Architecture

Two headers only: `objects.hpp` (every struct + every declaration; include it first) and
`printing.hpp` (output formats, includes `objects.hpp`). All 16 `.cpp` live in `src/`.

Flow: `main` (MPI init, parse, per-layer RNG seed) -> `Simulation::range`
(`simulation_range.cpp`: iteration x hysteresis pass x range1 x range2 x time loop) ->
`Simulation::timeStep` (`engine.cpp`: pulse -> `CalcHeff` -> adaptive dt -> LLG step ->
Joule heating -> save row).

| File | Role |
| --- | --- |
| `main.cpp` | MPI boot, JSON parse, opens `out.evol` (rank 0 only), seeds each layer's RNG |
| `read_input.cpp` | all `::fromJson` factories + `parseLayerList`/`parseJunctionList`; logs every value read |
| `simulation_range.cpp` | `Simulation::range`: the 4-level sweep loop, MPI sharding, `.evol` file setup |
| `engine.cpp` | `Simulation::timeStep` (the hot loop) and `reset_mag`; string-dispatches the solver |
| `CalcHeff.cpp` | effective field: thermal, anisotropy, demag, dipolar, RKKY, STT |
| `LLG_{Heun1,Heun2,RK4,Symplectic}.cpp` | integrators |
| `junction_init.cpp` | `Layer::updateResistance`, `Junction::initializeRangePoint`, `initializePulse` |
| `pulse.cpp` | `getValue` (linspace) and `pulse()` waveform types 1-4 |
| `Ndip.cpp` | analytic demag/dipolar tensors (`coeffdemag`) via surface integrals |
| `Rand_dist.cpp` | Numerical Recipes `ran2`/`gasdev` RNG |
| `path_resolver.cpp` | `resolvePath`: dot-path grammar over the object graph, `[n]` index / `[x]` wildcard |
| `evol_filename.cpp` | `.evol` filename builder with SI prefixes |
| `prints.cpp` | output formats + MPI row streaming |

## Conventions

- `/////`-style banner comments separating sections; explain **why** above a function, not
  what. This codebase is heavily commented about MPI and filesystem correctness — match it.
- Physics constants (`kB`, `gamma0`, `mu0`, `h_bar`, `e_charge`) are deliberately
  re-declared per translation unit rather than centralised. Keep values identical.
- `std::`-qualified names; only `main.cpp` and `Rand_dist.cpp` use `using namespace std`.
- Indentation is mixed: tabs in `pulse.cpp`, `junction_init.cpp`, `Ndip.cpp`,
  `Rand_dist.cpp`; 4-space everywhere else. Match the file you edit.
- Input validation is exception-based (`throw std::out_of_range` / `std::runtime_error`),
  never error codes. `jsonAt` throws on a missing required key; `jsonContains` guards optional ones.
- `CalcHeff` writes `layer.Heff` in place and solvers read it back — solvers mutate `layer.m`
  and re-call `CalcHeff` between stages, so statement order is load-bearing.
- Commits: lowercase imperative one-liners ("add Licence file", "readme first draft").

## Gotchas

- **Only rank 0 may write `1-State-Diagram.txt`.** Every other rank streams rows to it over
  MPI, because append-at-EOF is not atomic across handles on a shared/network filesystem.
  Keep `finalStateRowSize` / `packFinalStateRow` / `writeFinalStateHeader` /
  `writeFinalStateRow` column order in lockstep or received rows desync.
- Adding a `double` field to any struct also requires a case in that struct's `resolve()`
  overload in `path_resolver.cpp`, or it is invisible to `saving columns` and `range.type`.
- Adding a solver: declare in `objects.hpp`, implement as `LLG_*`, add a branch in
  `engine.cpp` (~line 82).
- `systemBias.Resistance` and the thermal field's `layer.dt` are *intentionally* one step
  behind. The in-code comments explain why — do not "fix" them silently.
- `range1` is split by integer division (`range1.step / np`), so unless `np` divides `step`,
  some points get computed twice and the tail is skipped.
- `.evol` filenames encode the range1 value; with `range1.flag` false every rank derives the
  same name and they clobber each other.
- `junction` is a flat JSON array consumed 3 blocks at a time (layer, Temperature,
  junction-bias), not an array of junction objects.
- Parsed but **not** implemented: the `Symplectic` solver (branch commented out in
  `engine.cpp`), field-like STT, Callen-Callen. Don't wire them up expecting results.
- `saving columns` entries must resolve to exactly one field — no `[x]` wildcards.
- The schema is untyped and there is no example file, so adding a *required* key breaks
  every existing input. Prefer `jsonContains` for new optional keys.
- **Keep `#include <json-c/json.h>` prefixed.** That is the *installed* layout. A vendored
  json-c generates `json.h` into its build tree and includes its siblings (`arraylist.h`,
  `json_object.h`, …) by plain quoted name from its *source* tree, so neither tree alone
  satisfies the prefixed form. `CMakeLists.txt` assembles the installed layout with symlinks
  under `build/json-c-vendored/include/json-c/`. Switching to `<json.h>` to "simplify" the
  vendored build breaks every system-installed json-c.
- Vendored json-c compiles as C, so `project()` must keep `C` in `LANGUAGES` even though
  Flexspin is pure C++. `BUILD_SHARED_LIBS` / `DISABLE_WERROR` are forced in the cache for
  the subpackage only; a dependency's `-Werror` must never be able to fail our build.

## Verify

No test suite, no CI, no linter, no formatter config. Smoke test only:

```
mpirun -np 2 ./build/Flexspin my_test.json    # tiny tmax, iterations=1, final-state saving=false
```

then check exit 0, `out.evol` content, and the `.evol` files. Single-file syntax check:
`mpicxx -std=c++11 -Isrc -fsyntax-only src/<file>.cpp`.

## Don't

- Commit `out.evol`, `1-State-Diagram.txt`, or `*K.evol` — they are written to the CWD (they
  are gitignored, but check you didn't `git add -f` them). Never commit `build/`.
- Add dependencies or headers beyond MPI + json-c + libm.
- Use `file(GLOB)` for the source list: CMake won't re-configure when a new `.cpp` appears,
  so a glob silently omits it. Add new sources to `CMakeLists.txt` explicitly.
