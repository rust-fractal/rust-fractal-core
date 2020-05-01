use std::time::Instant;
use rust_fractal::renderer::FractalRenderer;
use float_extended::float_extended::FloatExtended;

fn main() {
    println!("Mandelbrot Renderer");
    let center = (
        "-1.99996619445037030418434688506350579675531241540724851511761922944801584242342684381376129778868913812287046406560949864353810575744772166485672496092803920095332",
        "+0.00000000000000000000000000000000030013824367909383240724973039775924987346831190773335270174257280120474975614823581185647299288414075519224186504978181625478529");
    let zoom = "2.3620330788506154104770818136626E157";
   //
   // let center = (
   //     "-0.75",
   //     "0.0");
   // let zoom = "1E0";

   //  let center = (
   //     "-0.749999987350877481864108888020013802837969258626230419972587823828734338471228477079750588709551510361714463695461745528645748607681279674273355384334270208362211787387351792878073779449767292692440",
   //     "0.001000038688236832013124581230049849132759425863378894883003211011278068229551274712347955044740933397589760194545872789087012331273586364914484522575986336846199522726507205442204060303594956029930");
   // let zoom = 3.7E191;

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

    // let center = (
    //     "-1.7685304554715107439359975226423287998937254896566074423088827415719065562223946222490318219590860313620613389226996388",
    //     "0.0059445386042946892550190097831583153039526746987365071641624831873576346016435037438597231970921541448697901394734283301");
    // let zoom = 2.00000000000E113;

    // let center = (
    //     "-0.74977496303986055825690204249855929432518723234763348374372443797716597482713646451452176727030299366505937952743219986590574037612687173829307656886118099999999999999999999999",
    //     "-0.06638008421017090687747407166586874456952046519386815946449195256714882710925246541317934635067121915196760463342031197755105049354085337445454326026704199999999999999999999999");
    // let zoom = 4.70017300819E149;

    // let center = (
    //     "-1.99996619445037030418434688506350579675531241540724851511761922944801584242342684381376129778868913812287046406560949864353810575744772166485672496092803920095332176654846785710146671838403830455177782981641506704068474778442",
    //     "-0.000000000000000000000000000000000300138243679093832407249730397759249873468311907733352701742572801204749756148235811856472992884140755192241865049781816254785289455481579902975133989545509287024987599289178526901586954397865");
    // let zoom = 8.80307862787E167;

    // let center = (
    //     "-1.479796270901760024425279934770411525645551054432599517909807632824286254403907594526888466099962805022975196472549771681831234491695559852583955204986197762293872769474806903995564259040667568599770094383157857518790853771783314763231302599999999999999999999999998",
    //     "-0.001199443952813447746281973233374468444560314114132538362037569205657422216739564521471119107626453330996365067987088146663639996715939831819152248618042255824652268918299630897525386638029428706473919823922522752497780934312003352081931299999999999999999999999998");
    // let zoom = "3.27799999998E235";

    // let center = (
    //     "-0.74962449737876168207862620426836138684529527812364025754481571424558286479672801698750203976779135608148687747196595174858125388297577788573082753210469806739721377901297451172624762450940529814713048377873612297220313016947039287469999999999999999",
    //     "0.03427010874046016945172951545749474868051534672439111853174854738370595772935930171842999778713222794215453436628998200591914208207632674780780978496784843807690510401829301309069542198413574392425885166367231192087416338497065495509999999999999999");
    // let zoom = "7.58065474756E227";
    //
    // let center = (
    //     "-2.050559006608784126354192867479099536440881617905555543747957936759876883454787702544180601524947662497989006561199941686503035814450105179339348138642387820802369220633343749297359877883369716818327662391405150740569448193222460236250918606066946656283237234744284933274878594112062480369152009986108070633680147249036896791379281578431154790589622808147355741922535460027418509046159327830418839213357901328376216091465e-01",
    //     "6.714743990341850337919198679903175165708088772538990562494722242077350466986308963230331310791316412038624840828891415492559655884511352162724555529889388758002461559779419425873846339905834917826468747033849428932742988848171919023575874770141322083036622877010093575734754383955939380341331934522400057334366353132499917358060509093125946583754723092644426625535267499421294405389961174546102608924084746257236480182315e-01");
    // let zoom = "1.01671E401";

    // let center = (
    //     "-1.765255999389876234461152788785797644984633297655685666807481695846225912749655788093382329902715252392041460635885215850563792707279422539272781989295044035426373000613583591812021412405839940885783437366123652481019380496979520083898930887713503041323571033603267228732673881697627769631428767123079395338351895929545805513311149873254735634309800137254227353398493256006860933209587567509764016267722040584438664806755046349843252626021333596704283019326356534015635940495176966999611978220530750450843633573332991709349235725836809717722489976855541370575297961999390728282862817643005070619283157278020063873338224185426911490249877871908715519769993172252374558860975947234420276611064476299388955501971748803106758767381698583971521419007502885315003960294562814779263086726297113472688640615919829005087413164109846911484778894071443404207869354099285048615521762074034505836870947172363003569419561423084",
    //     "0.010448551737598706745099252605690625514295821842087741126875066858507024028622836332518011540214316990121188237193165533662492151873053535637312564647820075236423297779966109124237194221863354328221742652612069510430958445467568178885152483937848755343908808386977024289761709193968352919453948811678376125506326115996878653850756567342989058840713452162224873967866286764727771837281150865087243687129517443053689780343701902952426157131551709911536659251965545108966408202732131331366045475823757117158382885010926367768874015936942512124751884298663006249170055709799596487712912447394634628841851042104496855357298917956312668599435094932060198053752763342060348335338961707184554017308596827745257219064604143338288686451014411445123684128575901091265931918578117071646391893958981510330826990739118209425971517710287638838014409321765909156207492203435311198488275792102002106910744385277272767630353946649");
    // let zoom = "1.84488169976E890";

    // wiggly fungus spores
    let center = (
        "-0.2267766406738792238869135624479225519188035826708281258502443997140161010
53977243723620294561540538246226091646753716291965025945283313405812939567649164
52361982232203599593202957434690189191961247864733543762102752002017592954984861
92772278887878430351036469371829027133630840776427236400354986259421863107411213
25817524490907217561730464961939546452022764365871067724974385964779767954989349
38920671405021451930251382767250771898412213360984720364300006454269638404408517
94657206262472717467131890959883863577805840569063746152860486008421969838778780
99145894891433447761389044647748539645203345333422078197625491628967265718389221
44964818146115543749166217211082405974113708807317518936653409820885161298178304
00012590753830967255843485306892756674642256346655072217604137003717484195757936
82748557647207571432901605882657939755147101788446692307094358589890309883156893
01561385657042770968968300144987961141393421793589703004348337248018632614272227
82463739991608615887646230367881025825827344477685943533640098554924572952050641
20730377077726502722353931091060078568461148585210334548537120923902988084281147
74090692212634221582993286895040512008114816183469848055960129261029696004533444
20195490458563055860253655892486787760508428522990873741890869077442954459234571
27892977911954804221488267945432258431923698438283346329002887430137253426572309
88142781908020843388302899573862331968911117978082564655038027657803993941691718
16846767574201920635627208744412333432276302905277786787500582335288884286000987
61631770859189485519310056992658158699598853752270",
        "-1.1172396839363245916885582764842365153821675295821199843265465854803158962
50189264406156057191817849297272681286515999107473779581849650803960450730872327
82710121613416134156743550096796403977936151195888886076985793798795579476101212
03011651139903856048640702795757303764570742179468201491485788183126044330801498
73532115229542447757909367221100622837140990167041370344735702803399027683108943
14529930378897367489079610832146157903020045883604017904410952344929119156586774
37571288174459606782758145700990996030548586065791092932716537485787841131023660
99747478791966495905948007137235188983158460149463449922034330734982862417215211
94920671073109720323372489407914277193020201908834117690375864533131960533772917
23420000341442234991787129245296297322583623844610632239643805006686576274507312
50842266641375573779426661284050619522247569520340792371304202131566752454849217
96133563465638934200264488064679921718947552306971840958725886020337415034047287
44899320534879889216066670454205268778750305620255887360437005239989763262700795
86814574563207915555030173426494039202725529890039977786012777674851405459998351
37092039684336942081632361534302994306164613484897249342090617191002553681974137
20244104274254937842765584569655260372309009736741231415922832034935515699492950
50830004529526405168860372252803729345337450515287057640831305017979643140487986
77904990720304540601716185396133829287607447480986792934677435018584167609217420
35999826617596762133330744045066396240793475890186487207277589768298477420779690
53152566849560497837995206178903123387582461235780");
    let zoom = "6.36118599961E1520";

    // let center = (
    //     "-1.99999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999837485382669331911167033522316925483202642018042879214193261311019861349448517336056952538528636093441258319101518468947695052924660781278387926133612368999268678559118848357284879481437613458535637431620616343911585652641236167774248507873e+00",
    //     "1.37178052090609414145729128541350514525330839736695026708710030908880680107546729978680204559181338295039364573873196577319068006687895587728672855990826461331344843897843243379157581698185171855284714599499709801317677506729468150382879822938996916336090559197836550022966441254042166415444219105895943405807252402353950295056242178258440758e-211");
    // let zoom = "3.16254E321";

    // let center = (
    //     "-1.9415760324164384670258010207152963455836310703337636971553632153699953231748214265693439392809766889417990339495535138773600497608548184388478075120625274839384321083343796823744257245389639452038113844104471889025392801864432077245926602571922452821224354975005121169785109942357298421404262126590474492753926495821445368096060482566619420154152939020824848975839084899689147298562610145950467841704864062817746346590104375512093161193972226132292179170224058565214504713982409588034545465922379814453776711845006605474088172618318942723528198770270842997806802943245728179735224256977104485442662961648724699412079153183876067367809973539830625542616045035461492929366422046303130505120748785562498270824728380184790562313004257176283245764063853864890465963338139549448722374322325847127941993851820299660956298000794335295949327247034323965131694262466963369192793577044784245148385268311841785601289659869162188880350037606180185463830241628687029093400314547978708301711518714589761564546957594975557098179630412245677113459195619873119846884085035243183174578116395125830682200658696286847156509827144474739045902935784228299532539111233689746412869541803112073689174050805844695929308923321494334362799025940235848995",
    //     "0.0076169857177728112621246687606229814376622366239407338589772024020209499045976759825679034381527471870014537925705319003799456550890681045428875205696663254968549439704823347664655689021784716622056103356553863275910883556307935836668993039832203763244374511148663044834664066241136816718229114799629782439880378800360851777633838391773907366702059239604515317866206986243101122772063236617795555487457950265111373793657577209416496705008385380416635017788503016723893215173690307323803615044066492212811613536959807743048177509072251153187383237619353752702124457553324086682938780519928781817312106537599078988384676657029474483798795503599822706463824679809868792892065538509890720755789639529959046521787794638689936197549362968217869917857805310501253238233774705734055258130577943634688808709066929865581290058286559813519715385095276557193248610003430971112737927220580029937755075146758999539012278384936509827720561352748948755680390703939257600842162771410572206281828761630084038381535334400389424501769474806194318111015933255242113676030859575164143245275289818536248604205950967219024346107676764997921191983085684967333168366296725048419934415992920024669793857240558195607632144850125664688940387964132949595");
    // let zoom = "1.22472082766E1201";

    // let center = (
    //     "-1.47979627090176002442527993477041152564555105443259951790980763282428625440390759452688846609996280502297519647254977168183123449169555985258395520498619776229387276947480690399556425904066756859977009411263793785708219875610509744561935075696679379311446743203293779934657401446195720264889073827267581042992819255713966061960029600956488528039400110136025734674049423503781537254029386762974344877408936287464895271483783475726166964791550705530567412698465671131878693095307095060837046189521991933766461743373729335296823174958306177205874642213026688441839954998548228632916129130573679795132881474634980787481342928171810735776115910804591672359583232080214738756842020870108037672667163903668087702858242174790868294425161262311856593292749892164489981430120266269153332857313430037597424562936561035727564362469665765019460520403640792002805990427627245393548636119808720871232272469193417494166358296291456834369706631208728780590905622375604719695861753115611478829057158195173800020964886643199769833465772022595458243306922032161250434130129820404390020915953207304689209284167535233740807833771336648466892063079592833953210343248122153378823185483319111541760121993192532661867139807165430954987430491336008865394019103037240624123333041717070928774406557607194525549224437374332530755746472",
    //     "0.00119944395281344774628197323337446844456031411413253836203756920565742221673956452147111910762645333099636506798708814666363999671593983181915224861804225582465226891829963089752538663802942870647392089835850619549479809217517094705650800512589818845547651334193125696320994576851799996500166707629933364468442367780342152576726258487886944602820690396656443321807505102112291679480014093584294578224083185257286745793076460273084876732143503263079501882798843539950193014646525364332469881962607205094865925865351112950537542564548777525400698939569363893864361732830557576700300360031975362326696113051055653491167123233446246234988576885769895507893147644630600399665244899219737218184118073119994326089479053957247707087215783777715907242297581086324396128829436162952683551759071704991554446451841240100031034310522721953529643169811764354610360193619060982311619234487812680034284111950675036959896599313671848631132415414499871467394799245430624381488729998819803918236746598499649284882744230522324464162767597098718432389237393827299334316699893788725219849149531795033791868300431075317702471430360707557631331751720744479394157671912922443327484191678965009944496434405375655303850657861911827715432664919797198519668642248616949358878333185898189670994771943016313195294570183521403347534231");
    // let zoom = "2.0490799999984173E1298";

    let mut renderer = FractalRenderer::new(
        500,
        500,
        zoom,
        1500000,
        center.0,
        center.1,
        0.001,
        true,
        8
    );

    let time = Instant::now();
    renderer.render();
    println!("{:<14}{:>6} ms", "TOTAL", time.elapsed().as_millis());
}