use std::{cmp::max};
use std::time::Instant;

use std::sync::Arc;
use std::sync::atomic::{Ordering, AtomicUsize, AtomicBool};

use std::thread;
use std::time::Duration;

use parking_lot::Mutex;

use crate::util::{ComplexArbitrary, ComplexExtended, FloatArbitrary, FloatExtended, to_extended};
use crate::math::Reference;

pub struct BoxPeriod {
    pub box_center: ComplexExtended,
    pub points_z: [ComplexExtended; 4],
    pub points_c: [ComplexExtended; 4],
    pub period: usize
}

impl BoxPeriod {
    pub fn new(box_center: ComplexExtended, delta_box: [ComplexExtended; 4]) -> Self {
        BoxPeriod {
            box_center,
            points_z: delta_box.clone(),
            points_c: delta_box,
            period: 1
        }
    }

    pub fn crosses_origin(a: ComplexExtended, b: ComplexExtended) -> usize {
        if a.mantissa.im.signum() as i32 != b.mantissa.im.signum() as i32 {
            let mut d = b - a;
            d.reduce();

            let s = d.mantissa.im.signum() as i32;
            let t = (d.mantissa.im * a.mantissa.re - d.mantissa.re * a.mantissa.im).signum() as i32;

            (s == t) as usize
        } else {
            0
        }
    }

    pub fn points_surrond_origin(&self, reference_z: ComplexExtended) -> bool {
        let mut a = self.points_z[0] + reference_z;
        let mut b = self.points_z[1] + reference_z;
        let mut c = self.points_z[2] + reference_z;
        let mut d = self.points_z[3] + reference_z;

        a.reduce();
        b.reduce();
        c.reduce();
        d.reduce();

        let out = BoxPeriod::crosses_origin(a, b) 
            + BoxPeriod::crosses_origin(b, c) 
            + BoxPeriod::crosses_origin(c, d) 
            + BoxPeriod::crosses_origin(d, a);

        out & 1 == 1
    }

    pub fn find_period(&mut self, reference: &Reference) {
        while self.period < reference.current_iteration {
            if self.points_surrond_origin(reference.reference_data_extended[self.period - 1]) {
                // do some calculation here to work out a good estimate for the starting point
                // println!("{} {} {} {}", self.points_z[0], self.points_z[1], self.points_z[2], self.points_z[3]);

                // println!("box method period is: {}", self.period);
                break;
            };

            // TODO maybe add some glitch tests?
            for i in 0..4 {
                self.points_z[i] = self.points_z[i] * (reference.reference_data_extended[self.period - 1] * 2.0 + self.points_z[i]);
                self.points_z[i] += self.points_c[i];
                self.points_z[i].reduce();
            }

            self.period += 1;
        }
    }

    pub fn find_atom_domain_period(&mut self, reference: &Reference) {
        let mut z = self.box_center;

        let mut min_value = (reference.reference_data_extended[0] + z).norm();
        let mut min_value_at = 0;

        for i in 2..reference.current_iteration {
            z *= reference.reference_data_extended[i - 1] * 2.0 + self.box_center;
            z += self.box_center;
            z.reduce();

            let distance_from_origin = (reference.reference_data_extended[i] + z).norm();

            // point has escaped
            if distance_from_origin.to_float() > 1e16 {
                break;
            }

            if distance_from_origin < min_value {
                min_value = distance_from_origin;
                min_value_at = i;
            }
        }

        println!("atom domain period is: {}", min_value_at);
    }
}

pub struct BallMethod {
    pub radius: FloatExtended,
    pub radius_dz: FloatExtended,
    pub radius_z: FloatExtended,
    pub radius_radius: FloatExtended,
    pub ei: FloatExtended,
    pub point_c: ComplexExtended,
    pub point_z: ComplexExtended,
    pub point_dz: ComplexExtended,
    pub period: usize,
}

impl BallMethod {
    pub fn new(radius: FloatExtended, point_c: ComplexExtended) -> Self {
        BallMethod {
            radius,
            radius_dz: FloatExtended::new(0.0, 0),
            radius_z: FloatExtended::new(0.0, 0),
            radius_radius: FloatExtended::new(0.0, 0),
            ei: FloatExtended::new(0.0, 0),
            point_c: point_c.clone(),
            point_z: point_c,
            point_dz: ComplexExtended::new2(1.0, 0.0, 0),
            period: 1
        }
    }

    pub fn find_period(&mut self, reference: &Reference) {
        let time = Instant::now();

        self.period = 1;

        for k in 0..reference.current_iteration {
            self.radius_dz = self.point_dz.norm();
            self.radius_dz.reduce();

            self.radius_z = (self.point_z + reference.reference_data_extended[k]).norm();
            self.radius_z.reduce();

            self.radius_radius = self.radius * (self.radius_dz + self.radius * self.ei);
            self.radius_radius.reduce();

            if self.radius_radius > self.radius_z {
                self.period = k + 1;
                println!("ball method period is: {}", self.period);
                break;
            }

            self.ei = self.radius_dz * self.radius_dz + (2.0 * self.radius_z + self.radius * (2.0 * self.radius_dz + self.radius * self.ei)) * self.ei;
            self.ei.reduce();

            self.point_dz *= (self.point_z + reference.reference_data_extended[k]) * 2.0;
            self.point_dz += ComplexExtended::new2(1.0, 0.0, 0);
            self.point_dz.reduce();

            self.point_z *= reference.reference_data_extended[k] * 2.0 + self.point_z;
            self.point_z += self.point_c;
            self.point_z.reduce();
        }

        println!("period: {}ms", time.elapsed().as_millis());
    }
}

pub fn get_nucleus(mut guess_c: ComplexArbitrary, period: usize, iteration_flag: Arc<AtomicUsize>, progress_flag: Arc<AtomicUsize>, stop_flag: Arc<AtomicBool>, current_estimate: Arc<Mutex<ComplexExtended>>) -> Option<ComplexArbitrary> {
    let complex_precision = guess_c.prec();
    let precision = 3 * max(complex_precision.0, complex_precision.1);

    let start_c = guess_c.clone();

    guess_c.set_prec(precision);

    let temp = FloatArbitrary::with_val(precision, 2);
    let mut epsilon = temp.clone();

    epsilon.next_up();
    epsilon -= &temp;

    let epsilon_squared = epsilon.square();

    for _ in 0..64 {
        let mut z_current = ComplexArbitrary::new(precision);
        let mut derivative_current = ComplexArbitrary::new(precision);
        let mut h = ComplexArbitrary::with_val(precision, (1.0, 0.0));
        let mut dh = ComplexArbitrary::new(precision);
    
        for i in 1..=period {
            derivative_current *= 2.0;
            derivative_current *= &z_current;
            derivative_current += 1.0;
    
            z_current.square_mut();
            z_current += &guess_c;
    
            if i < period && period % i == 0 {
                h *= &z_current;
                dh += derivative_current.clone() / &z_current;
            }

            progress_flag.fetch_add(1, Ordering::SeqCst);

            if stop_flag.load(Ordering::SeqCst) {
                stop_flag.store(false, Ordering::SeqCst);
                return None
            };
        }
    
        dh *= &h;

        let df = ((derivative_current.clone() * &h - &z_current * &dh) / &h) / &h;

        let new_c = (-z_current.clone() / &h) / &df + &guess_c;

        let difference_from_start = new_c.clone() - &start_c;
        let difference_from_start_extended = to_extended(&difference_from_start);

        let mut current = current_estimate.lock();

        current.mantissa = difference_from_start_extended.mantissa;
        current.exponent = difference_from_start_extended.exponent;

        drop(current);

        thread::sleep(Duration::from_millis(500));

        let difference_norm = (new_c.clone() - &guess_c).norm();

        // need to set epsilon squared to floating point difference
        if difference_norm.real() <= &epsilon_squared {
            return Some(new_c);
        } else if difference_norm.real().is_infinite() || difference_norm.real().is_nan() {
            return None;
        }
        
        guess_c = new_c;
        iteration_flag.fetch_add(1, Ordering::SeqCst);
        progress_flag.store(0, Ordering::SeqCst);
    }

    None
}

pub fn get_nucleus_position(nucleus: ComplexArbitrary, period: usize) -> (FloatExtended, f64) {
    let time = Instant::now();

    let mut z = nucleus.clone();
    let mut l = ComplexExtended::new2(1.0, 0.0, 0);
    let mut b = ComplexExtended::new2(1.0, 0.0, 0);

    let temp = ComplexExtended::new2(1.0, 0.0, 0);
    let temp2 = FloatExtended::new(2.0, 0);

    l *= to_extended(&z) * 2.0;
    l.reduce();

    b += temp / l;
    b.reduce();
    

    for _ in 2..period {
        z.square_mut();
        z += &nucleus;

        l *= to_extended(&z) * 2.0;
        l.reduce();

        b += temp / l;
        b.reduce();
    }

    let mut size = temp / (b * l * l);
    size.reduce();

    let mut zoom = temp2 / size.norm();
    zoom.reduce();

    println!("nucleus size: {}ms", time.elapsed().as_millis());

    (zoom, size.mantissa.arg())
}
