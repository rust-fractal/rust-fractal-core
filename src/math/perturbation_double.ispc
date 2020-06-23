// Returns the count at which we exited the loop testing if the number is in the set
static inline int perturbation_double(
    int* iterations,
    const double delta_reference_re,
    const double delta_reference_im,
    double* delta_current_re,
    double* delta_current_im,
    int* glitched,
    int* escaped,
    const uniform int reference_current_iteration,
    const uniform int reference_start_iteration,
    const uniform double reference_re[],
    const uniform double reference_im[],
    const uniform double tolerance_check[]) {

	int i = 0;

    while (*iterations < reference_current_iteration) {
        double ref_re_tmp = reference_re[i];
        double ref_im_tmp = reference_im[i];


        double z_re = *delta_current_re + ref_re_tmp;
        double z_im = *delta_current_im + ref_im_tmp;

        double z_norm = z_re * z_re + z_im * z_im;

        if (z_norm < tolerance_check[i]) {
            *glitched = 1;
            *delta_current_re = z_re;
            *delta_current_im = z_im;
            break;
        }

        if (z_norm > 1e16) {
            *escaped = 1;
            *delta_current_re = z_re;
            *delta_current_im = z_im;
            break;
        }

        double temp_real = *delta_current_re;
        double temp_imag = *delta_current_im;

        *delta_current_re = 2.0 * (ref_re_tmp * temp_real - ref_im_tmp * temp_imag) + temp_real * temp_real - temp_imag * temp_imag + delta_reference_re;
        *delta_current_im = 2.0 * (ref_re_tmp * temp_imag + ref_im_tmp * temp_real + temp_imag * temp_real) + delta_reference_im;
        *iterations += 1;
        i += 1;
    }
}

task void perturbation_double_group(
    const uniform int task_size,
    uniform int iterations[],
    const uniform double delta_reference_re[],
    const uniform double delta_reference_im[],
    uniform double delta_current_re[],
    uniform double delta_current_im[],
    uniform int glitched[],
    uniform int escaped[],
    const uniform int reference_current_iteration,
    const uniform int reference_start_iteration,
    const uniform double reference_re[],
    const uniform double reference_im[],
    const uniform double tolerance_check[]) {

	foreach (i = 0 ... task_size) {
	    int index = taskIndex0 * task_size + i;

        perturbation_double(
            &iterations[index], 
            delta_reference_re[index], 
            delta_reference_im[index], 
            &delta_current_re[index],
            &delta_current_im[index],
            &glitched[index],
            &escaped[index],
            reference_current_iteration,
            reference_start_iteration,
            reference_re,
            reference_im,
            tolerance_check
        );
	}
}

export void perturbation_double_ispc(
    const uniform int task_size,
    const uniform int tasks,
    uniform int iterations[],
    const uniform double delta_reference_re[],
    const uniform double delta_reference_im[],
    uniform double delta_current_re[],
    uniform double delta_current_im[],
    uniform int glitched[],
    uniform int escaped[],
    const uniform int reference_current_iteration,
    const uniform int reference_start_iteration,
    const uniform double reference_re[],
    const uniform double reference_im[],
    const uniform double tolerance_check[]) {
	
    launch[tasks] perturbation_double_group(
        task_size, 
        iterations, 
        delta_reference_re, 
        delta_reference_im, 
        delta_current_re, 
        delta_current_im, 
        glitched, 
        escaped, 
        reference_current_iteration, 
        reference_start_iteration,
        reference_re,
        reference_im,
        tolerance_check
    );

    // foreach (i = 0 ... remaining) {
    //     int index = tasks * task_size + i;
    //     perturbation_f64(delta_real[index], delta_imag[index], max_iters, &counts[index], &smooth[index], &glitched[index], reference_real, reference_imag, tolerance_check);
    // }
}