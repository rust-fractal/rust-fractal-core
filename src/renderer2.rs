use crate::util::image::Image;
use crate::math::perturbation::Perturbation;
use crate::util::{ComplexArbitrary, ComplexFixed};
use crate::math::series_approximation::SeriesApproximation;
use crate::math::reference_arbitrary::Reference;

fn plot(image: &mut Image, perturbation: &Perturbation, dt: f64) {
    for i in 0..perturbation.iteration.len() {
        let (r, g, b) = if !perturbation.glitched[i] && !perturbation.escaped[i] {
            (255, 255, 255)
        } else if perturbation.escaped[i] {
            let de = 2.0 * perturbation.dz[i].norm() * perturbation.dz[i].norm().ln() / perturbation.dzdc[i].norm();
            let out = (255.0 * (de / dt).tanh()) as u8;
            (out, out, out)
        } else {
            (255, 0, 0)
        };

        image.plot(perturbation.image_x[i], perturbation.image_y[i], r, g, b);
    }
}

pub fn test() {
    let width = 320;
    let height = 320;
    let order = 10;

    let center = (
        "-1.99996619445037030418434688506350579675531241540724851511761922944801584242342684381376129778868913812287046406560949864353810575744772166485672496092803920095332",
        "+0.00000000000000000000000000000000030013824367909383240724973039775924987346831190773335270174257280120474975614823581185647299288414075519224186504978181625478529");
    let zoom = 2.3620330788506154104770818136626E157;

    let precision = 1000;
    let max_iterations = 100000;

    let c = ComplexArbitrary::with_val(
        precision as u32,
        ComplexArbitrary::parse("(".to_owned() + center.0 + "," + center.1 + ")").expect("Location is not valid!"));

    let aspect = width as f64 / height as f64;
    let dt =  (-2.0 * (4.0 / height as f64 - 2.0) / zoom) / height as f64;
    let top_left_delta = ComplexFixed::new((4.0 / width as f64 - 2.0) / zoom as f64 * aspect as f64, (4.0 / height as f64 - 2.0) / zoom as f64);

    let t_max = top_left_delta.norm();

    // get the maximum difference from the reference location (will be one of the corner points)

    let mut s = SeriesApproximation::new(c, dt, t_max, order);
    s.run();

    println!("skipped: {}", s.iteration);

    let mut reference = s.get_reference();
    // let mut reference = Reference::new(c.clone(), c, 0);
    reference.iterate(max_iterations);



    let mut image_x = Vec::new();
    let mut image_y = Vec::new();
    let mut iterations = Vec::new();
    let mut dc = Vec::new();
    let mut dz = Vec::new();
    let mut dzdc = Vec::new();

    for i in 0..width {
        for j in 0..height {
            let k = j * width + i;
            image_x.push(i as i32);
            image_y.push(j as i32);
            iterations.push(reference.start_iteration);
            let element = ComplexFixed::new(i as f64 * dt + top_left_delta.re, j as f64 * dt + top_left_delta.im);
            dc.push(element);
            dz.push(s.series_dz(element));
            dzdc.push(s.series_dzdc(element));
        }
    }

    let mut perturbation = Perturbation::new(image_x, image_y, iterations, dc, dz, dzdc);
    perturbation.iterate(&reference, reference.current_iteration);

    let mut image = Image::new(width as i32, height as i32);

    plot(&mut image, &perturbation, dt);

    image.save();
}