![image](render.png)

# rust-fractal
A mandelbrot fractal image generator featuring perturbation theory and series approximation. A high precision reference point is iterated at arbitrary precision and the differences from this are calculated in machine precision. This allows for a large speedup in image generation specifically at high zoom levels. This generator features:

- Perturbation based iteration count with glitch detection.
- Glitch fixing through automatic reference movement and recalculation.
- Series approximation calculation to skip (and approximate) large amounts of perturbation iterations.
- Multithreading of core perturbation loops through rayon.
- Multiple colouring methods including iteration, histogram and distance.
- Utilises scaling and mantissa-exponent based extended precision to allow for arbitrary zoom, whilst maintaining good performance. Verified to be working at depths exceeding E20000.

## Usage
You need to be able to compile the 'rug' crate. Look in the documentation for rug for more information on how to do this. Once all required dependencies have been installed, change the settings in the ```main.rs``` file located in ```src/```. Build and run the crate with:

```cargo run --release```
