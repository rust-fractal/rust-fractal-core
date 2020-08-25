use std::time::Instant;
use rust_fractal::renderer::FractalRenderer;

use clap::{crate_version, crate_name, crate_authors, crate_description, App, Arg};


fn main() {
    let matches = App::new(crate_name!())
        .about(crate_description!())
        .version(crate_version!())
        .author(crate_authors!())
        .arg(
            Arg::new("location")
                .short('l')
                .long("location")
                .value_name("FILE")
                .about("Sets the location file to use.")
                .takes_value(true)
                .required(true)
        )
        .arg(
            Arg::new("parameters")
                .short('p')
                .long("parameters")
                .value_name("FILE")
                .about("Sets the parameters file to use.")
                .takes_value(true)
                .required(true)
        ).get_matches();


    let mut settings = config::Config::default();

    if let Some(l) = matches.value_of("location") {
        settings.merge(config::File::with_name(l).required(true)).unwrap();
    };

    if let Some(p) = matches.value_of("parameters") {
        settings.merge(config::File::with_name(p).required(true)).unwrap();
    };

    let mut renderer = FractalRenderer::new(settings);

    let time = Instant::now();
    // renderer.render_sequence(2.0);
    // renderer.render("output/output".to_owned());
    renderer.render_sequence(2.0);
    println!("{:<14}{:>6} ms", "TOTAL", time.elapsed().as_millis());
}