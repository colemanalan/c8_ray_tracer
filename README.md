# C8 Ray Tracer

This is the implementation of the Ray Tracer in C8 for finding solutions
for in-ice problems.

## Installation

Currently I don't think I have been good at keep up the requirements. So feel free to make an PR with updates to the pyproject.toml file
with the things that I have missed.

### As a pure dependency

```bash
pip install git+https://github.com/colemanalan/c8_ray_tracer.git
```

### To have a local copy of the repository

```bash
git clone git@github.com:colemanalan/c8_ray_tracer.git
cd c8_ray_tracer
pip install .
```

## For linking with NRMC

There is a basic sketch of an implementation of a class that overloads the NRMC ray tracer base class
located at [./python/examples/nrmc_ray_tracing_class.py](nrmc_ray_tracing_class.py). You _should_ be able
to use this class directly in NRMC to provide ray tracing solutions. However, not all features are totally complete yet.
