use rust_fractal::renderer::FractalRenderer;
use clap::{crate_version, crate_name, crate_authors, crate_description, App, Arg};
use config::{Config, File};


fn main() {
    let matches = App::new(crate_name!())
        .about(crate_description!())
        .version(crate_version!())
        .author(crate_authors!())
        .arg(
            Arg::new("INPUT")
                .value_name("FILE")
                .about("Sets the location file to use.")
                .takes_value(true)
                .required(false)
        )
        .arg(
            Arg::new("options")
                .short('o')
                .long("options")
                .value_name("FILE")
                .about("Sets the options file to use.")
                .takes_value(true)
                .required(false)
        ).get_matches();


    let mut settings = Config::default();

    if let Some(l) = matches.value_of("INPUT") {
        settings.merge(File::with_name(l).required(true)).unwrap();
    };

    if let Some(p) = matches.value_of("options") {
        settings.merge(File::with_name(p).required(true)).unwrap();
    };

    let mut renderer = FractalRenderer::new(settings);
    // renderer.render("output/output".to_owned());
    renderer.render();
}