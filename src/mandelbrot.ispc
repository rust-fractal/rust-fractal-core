// Returns the count at which we exited the loop testing if the number is in the set
inline int perturbation(const double delta_0_real,
                        const double delta_0_imag,
                        const uniform int max_iters,
                        int* counts,
                        float* smooth,
                        int* glitched,
                        const uniform double reference_real[],
                        const uniform double reference_imag[],
                        const uniform double tolerance_check[]) {
	int i = 0;

    double delta_n_real = delta_0_real;
    double delta_n_imag = delta_0_imag;

    double reference_real_tmp = reference_real[0];
    double reference_imag_tmp = reference_imag[0];

    double temp1 = reference_real_tmp + reference_real_tmp;
    double temp2 = delta_n_imag + reference_imag_tmp;

	for (; i <= max_iters; ++i) {
	    unmasked {
                double temp_real = delta_n_real;

                delta_n_real *= temp1 + reference_real_tmp;
                delta_n_real -= (temp2 + reference_imag_tmp) * delta_n_imag;
                delta_n_real += delta_0_real;

                delta_n_imag *= 2 * temp1;
                delta_n_imag += 2 * reference_imag_tmp * temp_real;
                delta_n_imag += delta_0_imag;

                reference_real_tmp = reference_real[i + 1];
                reference_imag_tmp = reference_imag[i + 1];

                temp1 = delta_n_real + reference_real_tmp;
                temp2 = delta_n_imag + reference_imag_tmp;

                double norm_sqr = temp1 * temp1 + temp2 * temp2;
        }

        if (norm_sqr < tolerance_check[i + 1]) {
            *glitched = 1;
            break;
        }

        if (norm_sqr > 256.d) {
            *smooth = 1.0f - log(log((float)norm_sqr) / 1.38629436112) / 0.69314718056;
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
                             uniform int glitched[],
                             const uniform double reference_real[],
                             const uniform double reference_imag[],
                             const uniform double tolerance_check[]) {

	foreach (i = 0 ... task_size) {
	    int index = taskIndex0 * task_size + i;
        perturbation(delta_real[index], delta_imag[index], max_iters, &counts[index], &smooth[index], &glitched[index], reference_real, reference_imag, tolerance_check);
	}
}

export void mandelbrot_perturbation(const uniform int task_size,
                                    const uniform int tasks,
                                    const uniform int max_iters,
                                    const uniform double delta_real[],
                                    const uniform double delta_imag[],
                                    uniform int counts[],
                                    uniform float smooth[],
                                    uniform int glitched[],
                                    const uniform double reference_real[],
                                    const uniform double reference_imag[],
                                    const uniform double tolerance_check[]) {
	launch[tasks] perturbation_group(task_size, max_iters, delta_real, delta_imag, counts, smooth, glitched, reference_real, reference_imag, tolerance_check);
}