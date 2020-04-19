use crate::util::{ComplexArbitrary, ComplexFixed};
use float_extended::complex_extended::ComplexExtended;
use crate::math::reference_arbitrary::Reference;
use rug::Float;
use float_extended::float_extended::FloatExtended;

pub struct SeriesApproximation {
    pub iteration: i32,
    pub dt: f64,
    pub t_max: f64,
    pub z: ComplexArbitrary,
    pub c: ComplexArbitrary,
    pub order: usize,
    pub a: Vec<ComplexExtended>,
    pub b: Vec<ComplexExtended>,
}

impl SeriesApproximation {
    pub fn new(c: ComplexArbitrary, dt: f64, t_max: f64, order: usize) -> Self {
        assert!(dt > 0.0);
        assert!(t_max > 0.0);
        assert!(order >= 1);

        let mut a = vec![ComplexExtended::new(0.0, 0, 0.0, 0); order as usize + 2];
        let mut b = vec![ComplexExtended::new(0.0, 0, 0.0, 0); order as usize + 1];

        a[0] = ComplexExtended::new(c.real().to_f64(), 0, c.imag().to_f64(), 0);
        a[1] = ComplexExtended::new(1.0, 0, 0.0, 0);
        b[0] = a[1];

        SeriesApproximation {
            iteration: 1,
            dt,
            t_max,
            z: c.clone(),
            c,
            order,
            a,
            b
        }
    }

    pub fn step(&mut self) -> bool {
        let z_next = self.z.clone().square() + self.c.clone();
        let mut a_next = vec![ComplexExtended::new(0.0, 0, 0.0, 0); self.order as usize + 2];
        a_next[0] = ComplexExtended::new(z_next.real().to_f64(), 0, z_next.imag().to_f64(), 0);
        a_next[1] = 2.0 * self.a[0] * self.a[1] + 1.0;

        // Calculate the new a coefficents
        for k in 2..=self.order {
            let mut sum = ComplexExtended::new(0.0, 0, 0.0, 0);

            for j in 0..=((k - 1) / 2) {
                sum += self.a[j] * self.a[k - j];
            }
            sum *= 2.0;

            if k % 2 == 0 {
                sum += self.a[k/2] * self.a[k/2];
            }

            a_next[k] = sum;
        }

        let mut b_next = vec![ComplexExtended::new(0.0, 0, 0.0, 0); self.order as usize + 1];
        b_next[0] = a_next[1];

        // Calculate the derivative coefficents (b)
        for k in 1..=self.order {
            let mut sum = ComplexExtended::new(0.0, 0, 0.0, 0);

            for j in 0..=k {
                sum += self.a[j] * self.b[k - j];
            }

            b_next[k] = 2.0 * sum;
        }

        let mut sum = FloatExtended::new(0.0, 0);
        let mut t_max_n = FloatExtended::new(1.0, 0);

        for k in (self.order + 1)..=(2 * self.order) {
            let mut sum2 = ComplexExtended::new(0.0, 0, 0.0, 0);

            for j in (k - self.order)..=((k - 1) / 2) {
                sum2 += self.a[j] * self.a[k - j];
            }
            sum2 *= 2.0;

            if k % 2 == 0 {
                sum2 += self.a[k / 2] * self.a[k / 2];
            }
            sum += sum2.norm() * t_max_n;
            t_max_n *= self.t_max;
        }
        t_max_n = FloatExtended::new(1.0, 0);
        let mut sum2 = FloatExtended::new(0.0, 0);

        for j in 0..=self.order {
            sum2 += self.a[j].norm() * t_max_n;
            t_max_n *= self.t_max;
        }

        sum += 2.0 * self.a[self.order + 1].norm() * sum2;
        sum += 2.0 * self.a[self.order + 1].norm_square() * t_max_n;

        a_next[self.order + 1] = ComplexExtended::new(sum.mantissa, sum.exponent, 0.0, 0);

        // quaz0r
        let mut a_abs = vec![FloatExtended::new(0.0, 0); self.order as usize + 2];

        for k in 0..=(self.order + 1) {
            a_abs[k] = a_next[k].norm();
        }

        let mut zsa = FloatExtended::new(0.0, 0);
        let mut t_max_n = 1.0;

        for k in 0..=self.order {
            zsa += a_abs[k] * t_max_n;
            t_max_n *= self.t_max;
        }

        let mut dsa = FloatExtended::new(0.0, 0);
        t_max_n = 1.0;

        for k in 1..=self.order {
            dsa += k as f64 * a_abs[k] * t_max_n;
            t_max_n *= self.t_max;
        }

        let mut rhs = zsa * dsa * self.dt;
        let mut lhs = a_abs[self.order + 1];

        for k in 0..self.order {
            lhs *= self.t_max;
        }

        // if lhs <= rhs {
        if self.iteration <= 10000 {
            self.iteration += 1;
            self.z = z_next;
            self.a = a_next;
            self.b = b_next;
            true
        } else {
            false
        }
    }

    pub fn run(&mut self) {
        while self.step() {

        };
    }

    pub fn get_reference(&self) -> Reference {
        Reference::new(self.z.clone(), self.c.clone(), self.iteration)
    }

    pub fn series_dz(&self, dc: ComplexFixed<f64>) -> ComplexFixed<f64> {
        let mut dc0 = ComplexExtended::new(dc.re, 0, dc.im, 0);
        let mut dcn = dc0;
        let mut sum = ComplexExtended::new(0.0, 0, 0.0, 0);

        for k in 1..=self.order {
            sum += self.a[k] * dcn;
            dcn *= dc0;
        };

        sum.to_float()
    }

    pub fn series_dzdc(&self, dc: ComplexFixed<f64>) -> ComplexFixed<f64> {
        let mut dc0 = ComplexExtended::new(1.0, 0, 0.0, 0);
        let mut dcn = dc0;
        let mut sum = ComplexExtended::new(0.0, 0, 0.0, 0);

        for k in 0..=self.order {
            sum += self.b[k] * dcn;
            dcn *= ComplexExtended::new(dc.re, 0, dc.im, 0);
        };

        sum.to_float()
    }
}