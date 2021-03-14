use crate::util::ComplexExtended;
use crate::math::Reference;

pub struct BoxPeriod {
    pub points_z: Vec<ComplexExtended>,
    pub points_c: Vec<ComplexExtended>,
    pub period: usize
}

impl BoxPeriod {
    pub fn new(mut delta_top_left: ComplexExtended) -> Self {
        delta_top_left.reduce();

        let mut points = Vec::new();

        points.push(delta_top_left);
        delta_top_left.mantissa.re *= -1.0;
        points.push(delta_top_left);
        delta_top_left.mantissa.im *= -1.0;
        points.push(delta_top_left);
        delta_top_left.mantissa.re *= -1.0;
        points.push(delta_top_left);

        BoxPeriod {
            points_z: points.clone(),
            points_c: points,
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

    pub fn get_nucleus(&mut self, reference: &Reference, mut guess_c: ComplexExtended) {
        for _ in 0..16 {
            let mut z = ComplexExtended::new2(0.0, 0.0, 0);
            let mut dc = ComplexExtended::new2(0.0, 0.0, 0);

            let mut perturbation_z = ComplexExtended::new2(0.0, 0.0, 0);
    
            let mut h = ComplexExtended::new2(1.0, 0.0, 0);
            let mut dh = ComplexExtended::new2(0.0, 0.0, 0);
    
            // TODO step using perturbation
            for i in 1..=self.period {
                dc = z * dc * 2.0 + ComplexExtended::new2(1.0, 0.0, 0);
                z = z * z + guess_c;
    
                if i < self.period && self.period % i == 0 {
                    h = h * z;
                    dh = dh + dc / z;
                }
            }
    
            dh = dh * h;
            
            let g = z;
            let dg = dc;
            
            let f = g / h;
            let df = (dg * h - g * dh) / (h * h);
    
            let new_c = guess_c - f / df;
    
            let d = new_c - guess_c;
    
            if d.to_float().norm_sqr() < 1e-30 {
                break;
            }

            guess_c = new_c;
        }
    }
}



