// Returns the count at which we exited the loop testing if the number is in the set
inline int perturbation(double delta_0_real,
                        double delta_0_imag,
                        int max_iters,
                        int* counts,
                        float* smooth,
                        bool* glitched,
                        const uniform double reference_real[],
                        const uniform double reference_imag[]) {
	int i = 0;

    double delta_n_real = delta_0_real;
    double delta_n_imag = delta_0_imag;

    double reference_real_tmp = reference_real[0];
    double reference_imag_tmp = reference_imag[0];

	for (; i <= max_iters; ++i) {
	    unmasked {
                double temp_real = delta_n_real;

                delta_n_real = 2.0 * (reference_real_tmp * delta_n_real - reference_imag_tmp * delta_n_imag)
                        + delta_n_real * delta_n_real - delta_n_imag * delta_n_imag
                        + delta_0_real;
                delta_n_imag = 2.0 * (reference_real_tmp * delta_n_imag + reference_imag_tmp * temp_real + delta_n_imag * temp_real)
                        + delta_0_imag;

                reference_real_tmp = reference_real[i + 1];
                reference_imag_tmp = reference_imag[i + 1];

                double temp1 = delta_n_real + reference_real_tmp;
                double temp2 = delta_n_imag + reference_imag_tmp;
        };

        double norm_sqr = temp1 * temp1 + temp2 * temp2;

        if (norm_sqr > 256.d) {
            *smooth = 1.0 - log(log((float)norm_sqr) / 1.38629436112) / 0.69314718056;
			break;
		}
	}
    *counts = i;
}

task void perturbation_group(const uniform int task_size,
                             const uniform int max_iters,
                             const uniform double delta_real[],
                             const uniform double delta_imag[],
                             uniform int counts[],
                             uniform float smooth[],
                             uniform bool glitched[],
                             const uniform double reference_real[],
                             const uniform double reference_imag[]) {
	foreach (i = 0 ... task_size) {
	    int index = taskIndex0 * task_size + i;
        perturbation(delta_real[index], delta_imag[index], max_iters, &counts[index], &smooth[index], &glitched[index], reference_real, reference_imag);
	}
}

export void mandelbrot_perturbation(const uniform int task_size,
                                    const uniform int tasks,
                                    const uniform int max_iters,
                                    const uniform double delta_real[],
                                    const uniform double delta_imag[],
                                    uniform int counts[],
                                    uniform float smooth[],
                                    uniform bool glitched[],
                                    const uniform double reference_real[],
                                    const uniform double reference_imag[]) {
	launch[tasks] perturbation_group(task_size, max_iters, delta_real, delta_imag, counts, smooth, glitched, reference_real, reference_imag);
}