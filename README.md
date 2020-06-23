![image](render.png)

# rust-fractal
A mandelbrot fractal image generator featuring perturbation theory and series approximation. A high precision reference point is iterated at arbitrary precision and the differences from this are calculated in machine precision. This allows for a large speedup in image generation specifically at high zoom levels. This generator features:

- Perturbation based iteration count with glitch detection.
- Glitch fixing through automatic reference movement and recalculation.
- Series approximation calculation to skip a large number of iterations.
- Multithreading of core perturbation loops through rayon.
- Multiple colouring methods including iteration, histogram and distance.
- Adaptive precision allowing for both double, and extended double (mantissa-exponent) types for deltas.
- Up to E280 zoom level for double precision, extended to E10000+.

## Usage
You need to be able to complie the 'rug' crate. Look in the documentation for rug for more information on how to do this. Once all required dependencies have been installed, change the settings in the ```main.rs``` file located in ```src/```. Build and run the crate with:

```cargo run --release```
