# GoldGates.jl
Repository containing scripts and utilities for generating gate solutions on Gold System using the [MSSim.jl package](https://github.com/euriqa-brassboard/MSSim.jl).
Activate this package and install all the dependencies by

```julia
using Pkg
Pkg.activate(".") # or the path that points to the directory this is cloned to
Pkg.instantiate()
```
`

## Project Structure
Below is a description of all the subdirectories

### system params/
Contains JSON files that detail system parameters relevant to optimization for gates.

### notebooks/
All notebook related to gate generation for specific purposes or system parameters.
Notebooks should start with a 6 digit date (MMDDYY) and contain code generally for one system configuration.

### src/
Contains the GoldGates module which has any helper code/scripts for generating, analyzing, reading, or writing gold solutions.

### out/
Contains all the output files from solution generations.

## How to Use this Project:
