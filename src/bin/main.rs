use std::time::Instant;
use rust_fractal::renderer2::FractalRenderer;

fn main() {
    println!("Mandelbrot Renderer");
    // let center = (
    //     "-1.99996619445037030418434688506350579675531241540724851511761922944801584242342684381376129778868913812287046406560949864353810575744772166485672496092803920095332",
    //     "+0.00000000000000000000000000000000030013824367909383240724973039775924987346831190773335270174257280120474975614823581185647299288414075519224186504978181625478529");
    // let zoom = 2.3620330788506154104770818136626E157;

   // let center = (
   //     "0.0",
   //     "0.0");
   // let zoom = 1.0;

    let center = (
       "-0.749999987350877481864108888020013802837969258626230419972587823828734338471228477079750588709551510361714463695461745528645748607681279674273355384334270208362211787387351792878073779449767292692440",
       "0.001000038688236832013124581230049849132759425863378894883003211011278068229551274712347955044740933397589760194545872789087012331273586364914484522575986336846199522726507205442204060303594956029930");
   let zoom = 3.7E191;

   // let center = (
   //     "-1.689573325432279612075987864633746594591438093139394112928000260",
   //     "0.000000000000000000000000000000000145514706909258179374258");
   // let zoom = 1.2980742146337048E34;

    // let center = (
    //     "-1.999999999138270118722935763129859470012913240693218269085000",
    //     "-0.000000000000000000322680215517275822769282166130504217892606");
    // let zoom = 1.63e030;

    // let center = (
    //     "-1.76904090125288168240284276810365472222460651367280300612031519864776551390840630119504817227707689276670981344397551593371805167279428144061791279408701383958080000000000000000000000002e+00",
    //     "-3.1605400272558568759745268636087334444350515177571908056976983793515506960262149704114708190147368048572976192049103272074188848768835600195188510326583157513e-03");
    // let zoom = 7.5E139;



    let mut renderer = FractalRenderer::new(
        2000,
        2000,
        zoom,
        5000000,
        center.0,
        center.1,
        1000,
        0.001,
        false,
        64
    );

    let time = Instant::now();
    renderer.render();
    println!("{:<14}{:>6} ms", "TOTAL", time.elapsed().as_millis());
}