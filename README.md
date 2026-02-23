# Context

This project is made for the course **Introduction to Molecular Simulation** (Introduction à la simulation moléculaire) from the **Paris-Saclay Master** in **High Performance Computing and Simulation** (CHPS)

# Dependencies

- Rust and cargo (https://www.rust-lang.org/tools/install)

# How to run

Run `cargo run --release` in the project root to execute the simulation, it takes in a configuration file and produces a plot of the energy evolution over time. 
Energy and temperature will be printed at each step, and you can modify the `main` function in `src/main.rs` to change the input file and output plot path.

There is also a small visualizer in `visualization`

# Documentation

The code is documented with Rust doc comments, you can generate the documentation with `cargo doc --open`.
Most of the documentation is located inside the `System` struct because Rust groups documentations by structs instead of by files/modules, so for example the documentation for the `movement` module is located in the `System` struct because it contains most of the functions of the `movement` module.

# Features

- Lennard-Jones potential
- Periodic boundary conditions
- Smooth cut-off for potential and forces
- Velocity-Verlet integration
- Berendsen thermostat
- Verlet neighbor lists for efficient force computation

# Missing features

- CLI to specify input/output paths and simulation parameters
- Parallelization
- Buffer reusage for neighbor lists and forces