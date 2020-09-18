![image](render.png)

# rust-fractal
A mandelbrot fractal renderer implementing both perturbation and series approximation. A reference point is iterated at high-precision, arbitrary precision and differences from this are calculated in machine precision. This allows for a large reduction in computation required to render and image, especially at high zoom levels. This generator features:

- Perturbation based iteration with glitch detection.
- Glitch correction through automatic reference movement and recalculation.
- Series approximation calculation to skip (and approximate) large amounts of perturbation iterations.
- Probe based method to determine series approximation skip.
- Multithreading of core loops through rayon.
- Configurable location and rendering options.
- Multiple save formats including PNG, EXR and KFR.
- Utilises scaling and mantissa-exponent based extended precision to allow for arbitrary zoom, whilst maintaining good performance. Verified to be working at depths exceeding E50000. Theoretically, this is only limited by MPFR's precision.

## Usage
You need to be able to compile the 'rug' crate which requires a rust GNU toolchain. Look in the documentation for rug for more information on how to do this. Once all required dependencies have been installed, build the crate with:

```cargo build --release```

The created binary can then be used with some of the config files provided. For help with the options do:

```./target/release/main --help```

By default, an image of the set is placed in the ```./output``` folder with it's specific details. A typical call with custom settings would look like:

```./target/release/main -o default.toml locations/flake.toml```

## Acknowledgements
- claude (blog, Kalles Fraktaler 2+)
- pauldelbrot (glitch detection, nanoscope)
- knightly (superMB)



