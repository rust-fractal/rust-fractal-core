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

## Compiling
You need to be able to compile the 'rug' crate which requires a rust GNU toolchain. Look in the documentation for rug for more information on how to do this. Once all required dependencies have been installed, build the crate with:

```cargo build --release```

## Usage
Information on the flags which can be passed to the rendered can be found with the command ```rust-fractal --help```. The renderer takes .toml files as input. There are two seperate files which can be defined to render an image, the options file and the location file. Settings in these files can be changed in order to change the output of the program. By default, there are 3 options files provided, which are:

- ```low.toml```: low quality settings for fast rendering and preview.
- ```default.toml```: settings that are used by default if no config file is provided.
- ```high.toml```: higher quality settings for final rendering.

Location files contain information on the specific location to be rendered, including the location, zoom level and rotation. Some examples of these files are stored in the ```./locations``` directory. A typical call to the renderer would then look like:

- Linux: ```rust-fractal -o default.toml locations/flake.toml```
- Windows: ```rust-fractal.exe -o default.toml locations/flake.toml```

Output images are placed in the ```./output``` folder.

## Acknowledgements
- claude (blog, Kalles Fraktaler 2+)
- pauldelbrot (glitch detection, nanoscope)
- knighty (superMB)



