% Plots vertically integrated LWP on the 2 x 2 grid

function [outputArg1,outputArg2] = plot_lwp(file, folder, nam, params)

% Load in the plot paramters 
output = params.output;
i = params.index;
export = params.export;

fname = strcat(folder, file); 

%% ========================================================================
% read all relevant 3-D variables

time = nam.dt * i;                                                         % Time [hours]

ph =ncread(fname,'PH' );                                                   % geopotential perturbation [m2/s2]
phb=ncread(fname,'PHB');                                                   % base geopotential [m2/s2)
p  =ncread(fname,'P'  );                                                   % pressure perturbation [Pa]
pb =ncread(fname,'PB' );                                                   % base pressure [Pa]
th =ncread(fname,'T'  )+nam.T0;                                            % Potential temperature [K]

qv=ncread(fname,'QVAPOR');                                                 % Water vapor mixing ratio [kg/kg]
qc=ncread(fname,'QCLOUD');                                                 % Cloud water mixing ratio [kg/kg]
%qr=ncread(fname,'QRAIN');                                                 % Rain water mixing ratio [kg/kg]
qi=ncread(fname,'QICE');                                                   % Ice mixing ratio [kg/kg]

%%

%% ========================================================================
% initial calculations

s=size(th);
n=s(1);
m=s(2);
l=s(3);
nm=n*m;

HS=mean(reshape(ph+phb,nm,l+1));
H=0.5*(HS(:,1:end-1)+HS(:,2:end));
Z=H./nam.g';                                                               % height at mass-levels [m]
p=p+pb;                                                                    % pressure
exn=(p/nam.P0).^(nam.R/nam.cp);                                            % exner function
qt=qv+qc+qi;                                                               % total water mixing ratio [kg/kg] (no precipitating elements)

t=th.*exn;                                                                 % temperature
tv=t.*(1+0.608*qv);                                                        % virtual temperature, bouyancy is tv - ql (eq. 1 Marcin)
rho=p./(nam.R*tv);                                                         % density

qtrho = qt.*rho;

LWP = trapz(Z',qtrho,3);                                                   % vertically integrated LWP


%% Plot Liquid Water Path

% Get the screen size
scrsz = get(0,'ScreenSize');

figure1 = figure('Position', [0 0 scrsz(3)*0.5 scrsz(3)*0.4]);
axes1 = axes('Parent',figure1);
hold(axes1,'on');

x_km = linspace(0, (n*nam.dx)/1000, n);
y_km = linspace(0, (n*nam.dy)/1000, n);

image(x_km, y_km, LWP,'CDataMapping','scaled'); hold on; grid on
ylabel('y [km]','LineWidth',1,'FontSize',14);
xlabel('x [km]','LineWidth',1,'FontSize',14);
title(num2str(time,'%.1f')+" hours")


% colormap = cmocean('ice'); % can put in a nicer colourmap later
cbar = colorbar;
xlim([min(x_km) max(x_km)])
ylim([min(y_km) max(y_km)])

box(axes1,'on');

% Set the remaining axes properties
set(axes1,'CLim',[40 60],'ColorScale','log','Colormap',...
    [0.0153116743554373 0.0225205938869953 0.0727287373590776;0.01800549591959 0.0254455160838977 0.0784187911682551;0.0209013300620317 0.0285265224507104 0.0840777157742097;0.0239981865058799 0.0317626432721827 0.0897075041635133;0.0272977574900411 0.0351490949685136 0.0953264723049982;0.0308042814605239 0.0386770816467026 0.100966267344003;0.0345098437341545 0.0422995563671478 0.106582816574358;0.0384135978411689 0.04587404239999 0.112177634064342;0.0424578961639736 0.0493996367983215 0.117796077309418;0.0464577753234479 0.0528848984910881 0.123417603628675;0.0504162176208768 0.0563378882456592 0.129021683092768;0.0543377951687331 0.0597573750328046 0.134623924667266;0.0582300025512884 0.0631358584082301 0.140265173920284;0.0620864735127598 0.0664886588959431 0.145892172402662;0.0659089082757606 0.0698175953283308 0.151505618549782;0.0697066971255817 0.0731090980407833 0.157170105577591;0.0734736080565012 0.0763781558146736 0.162830826512931;0.0772098218588232 0.0796282040661184 0.168480039405165;0.0809219057254782 0.0828485657183703 0.174169730470864;0.0846072252062268 0.086047198208808 0.179872577391429;0.0882640798355726 0.0892308723058959 0.185565254824214;0.0918974885664391 0.0923899162782522 0.191295325154019;0.095506056373538 0.0955293654823137 0.197046634472381;0.0990876019983124 0.0986573559892884 0.202788508027131;0.102646427697998 0.101763283260989 0.208572900686604;0.106180748129312 0.104853398782841 0.214377161960038;0.109688902376089 0.107935170201639 0.220172187124978;0.113174816829629 0.11099513380448 0.226023644343488;0.116635565471796 0.114044838311235 0.231883631858028;0.120070677078405 0.117088369338741 0.237737234513734;0.123483112691702 0.12010957852661 0.243664085460049;0.126869476918905 0.123127197370705 0.249580923435607;0.13023079358556 0.12613482732304 0.2555217697948;0.13356695365438 0.129129338194395 0.261506047383826;0.136876873589179 0.132122839019011 0.267479226823866;0.140160722230301 0.135100775126727 0.273516463830988;0.143417719783112 0.138076480854388 0.279556434897203;0.146647490226777 0.141046835919909 0.285618379003358;0.149848965736202 0.144008282628215 0.291724668757343;0.153022936107079 0.14697296430958 0.297817052765649;0.156165889104442 0.149923783987301 0.30398847518814;0.15928036968221 0.152879682418173 0.310144105187575;0.162362751309763 0.155828761226028 0.316351221127571;0.165414243936029 0.158779621285636 0.322568010825442;0.168432796509061 0.161729308406336 0.328814893078184;0.17141745326567 0.164679056030592 0.3350897686327;0.174368256224111 0.167631910590234 0.341379882095892;0.177281895377815 0.170584638243905 0.347708880914031;0.180160363423371 0.173543507516549 0.354044758307851;0.18299849497697 0.176503483429546 0.360422974163663;0.185799379924209 0.179471545482064 0.366806039921701;0.188557237056137 0.182443215954959 0.373227571969878;0.191274641082766 0.18542401466967 0.379657962810939;0.193947218744642 0.188411998473316 0.386115871001139;0.196574621820292 0.191409468731977 0.392592257954186;0.199156723720246 0.194418542575984 0.399078527624376;0.201687028061793 0.197437023129919 0.405597945214553;0.204173323195979 0.2004720971994 0.412103463877426;0.206598917596433 0.203516328594515 0.418661153097802;0.208976998112478 0.20657922518971 0.425203222940029;0.211296854169097 0.20965751618649 0.431764975135799;0.213558006203499 0.212753402887682 0.438337407058063;0.215765485712072 0.215870595983409 0.44489463585537;0.217902746703169 0.219005689637082 0.451482166635808;0.21998015574811 0.222164081757603 0.458056515039495;0.221997032731354 0.22534706636685 0.464612933092324;0.223935940016192 0.22855284524071 0.471190050407052;0.225809325329802 0.231785965709774 0.477744572372931;0.227615038623459 0.235047480233127 0.484274917938962;0.229341112984248 0.238337569159647 0.490800122509867;0.230990510277274 0.241658462682364 0.49730374060116;0.232565389994512 0.245011820726866 0.50377331532637;0.234063421706667 0.248398801117271 0.510206196911593;0.235473939185548 0.25182056094885 0.516613364973523;0.236801026810658 0.25527868910579 0.522978761262223;0.238045905601405 0.258774353795457 0.529292995816387;0.239206916672632 0.26230861756345 0.535551744864301;0.240282648714022 0.265882488863808 0.541750276269916;0.241269037903654 0.269497145174665 0.547887066806436;0.242164653482724 0.273153607535969 0.553956728335319;0.24297376639755 0.276852195096576 0.559947885292339;0.243696427132769 0.280593400273896 0.565854993652327;0.24433314064881 0.284377550064945 0.571672398667092;0.244884901991187 0.288204786196639 0.577394413237392;0.2453532261045 0.292075047930716 0.583015400606626;0.245740170742314 0.295988058192435 0.588529859447889;0.246048351487498 0.29994331361254 0.59393250922176;0.246280948074455 0.303940078957737 0.599218373567152;0.246441701430281 0.307977386280304 0.604382859465438;0.246534901117969 0.312054038948954 0.609421829997155;0.246565363158923 0.316168620539951 0.614331668692064;0.246538398519473 0.320319508380464 0.619109333750879;0.246459772849957 0.324504891356841 0.623752400776989;0.24633565834843 0.328722791440109 0.628259093078747;0.246172578868182 0.332971088249406 0.632628299062546;0.245977349585916 0.337247545878876 0.636859576706593;0.245757012685783 0.341549841159386 0.640953145557386;0.245518770588334 0.345875592514631 0.644909867100444;0.245269918261997 0.350222388599969 0.648731214702924;0.245017776101522 0.354587815976946 0.652419234593584;0.244769624750394 0.358969485170182 0.655976499526834;0.244532643092561 0.363365054567717 0.659406056870784;0.244313850454874 0.367772251752496 0.662711372868709;0.244120053857863 0.372188891982797 0.665896274758654;0.243957800941002 0.376612893665593 0.668964892309873;0.243833338980068 0.381042290782907 0.671921600162475;0.243752716608817 0.385475204852549 0.674770992721115;0.243721586806758 0.389909896583832 0.677517810493398;0.243744969103434 0.394344830666655 0.680166822877906;0.243827594298859 0.398778578587252 0.682722859277318;0.24397378347042 0.403209848198766 0.685190748242695;0.244187443496841 0.407637481219901 0.687575282620461;0.244472067164551 0.412060449304047 0.689881189397166;0.244830737373611 0.416477848963876 0.692113104094962;0.245266134952096 0.420888895616558 0.694275549483495;0.245780549592735 0.425292916990024 0.696372918310862;0.246375893440102 0.429689346103265 0.698409459714725;0.247053716878776 0.434077714005146 0.700389268951508;0.247815226100181 0.438457642427967 0.702316280073691;0.248661302056467 0.442828836484956 0.704194261189378;0.249592520442447 0.447191077515817 0.706026811951787;0.250609172380137 0.451544216161703 0.7078173629464;0.251711285513826 0.455888165730907 0.70956917666821;0.252898645256474 0.460222895899027 0.711285349808729;0.254170815959892 0.464548426772544 0.712968816600996;0.255527161811432 0.468864823332328 0.714622352999297;0.25696686728864 0.473172190263367 0.716248581498091;0.258488957030131 0.477470667168967 0.717849976420898;0.260092315006059 0.481760424161218 0.719428869534396;0.261775702894665 0.486041657814752 0.720987455865331;0.263537777592665 0.490314587467271 0.722527799618077;0.265377107806529 0.494579451847879 0.724051840108644;0.267292189689159 0.498836506012735 0.725561397646821;0.269281461502063 0.503086018566678 0.727058179311933;0.271343317296859 0.50732826914926 0.72854378457966;0.273476119622043 0.511563546163765 0.730019710767536;0.275678211271224 0.515792144728339 0.731487358275442;0.27794792609795 0.52001436482909 0.732948035604621;0.280283598929434 0.524230509655942 0.734402964144834;0.282683574617556 0.528440884103063 0.73585328272421;0.285146216270106 0.532645793416741 0.737300051920429;0.287669912708876 0.536845541974728 0.738744258135122;0.290253085203672 0.54104043218213 0.740186817435976;0.292894193533009 0.545230763469999 0.741628579173082;0.295591741423026 0.549416831383799 0.743070329377597;0.298344281416412 0.553598926749889 0.744512793952004;0.301150419222662 0.557777334909026 0.745956641662107;0.304008817600157 0.561952335006755 0.747402486941509;0.306918199819239 0.56612419933129 0.748850892519755;0.309877352754004 0.570293192690155 0.750302371885581;0.312885129648551 0.574459571817497 0.751757391596862;0.315940452601721 0.578623584804508 0.753216373448947;0.319042314811976 0.582785470545864 0.754679696513081;0.322189782622161 0.586945458195537 0.756147699056641;0.325381997401538 0.591103766625661 0.757620680356895;0.328618177300318 0.595260603882464 0.759098902420032;0.331897618909697 0.59941616663355 0.760582591617254;0.335219698858202 0.603570639600989 0.762071940249804;0.338583875372922 0.607724194974922 0.763567108054989;0.341989689832065 0.611876991802498 0.765068223665432;0.345436768333062 0.616029175347089 0.766575386034115;0.348924823298166 0.620180876412899 0.76808866583814;0.352453655137278 0.624332210630138 0.76960810687461;0.356023153985287 0.628483277696078 0.771133727462591;0.359633301528688 0.632634160567456 0.77266552186578;0.363284172933444 0.636784924599826 0.774203461751298;0.366975938883185 0.640935616629703 0.775747497700877;0.370708867733294 0.645086263995593 0.777297560791693;0.374483327782751 0.649236873494405 0.778853564265204;0.378299789661268 0.653387430270163 0.780415405303493;0.38215882882428 0.657537896632609 0.781982966933888;0.386061128142583 0.66168821080404 0.783556120083949;0.390007478849065 0.665838331031395 0.785134355547658;0.393998816169696 0.669988229416051 0.78671676205655;0.398036194326838 0.674137642824924 0.788304118055203;0.402120766656388 0.678286397063216 0.789896260397811;0.406253808869597 0.682434286049913 0.791493027957701;0.410436721548229 0.686581070131777 0.793094265309994;0.414671032368832 0.690726474360986 0.794699826775071;0.418958397945686 0.694870186754463 0.796309580853497;0.423300761549649 0.699011946260045 0.797922372037922;0.427700313139677 0.703151440221996 0.799536931816084;0.432158877273799 0.70728805165131 0.801155146437994;0.436678659484717 0.711421296390224 0.80277699667406;0.441261993497231 0.715550644554885 0.804402509512385;0.445911751923532 0.719675607660483 0.806030220928078;0.450631097766465 0.723795610762785 0.807658697998032;0.455422053360729 0.727909790987073 0.809291018044729;0.460287429613298 0.732017417469098 0.810927497554066;0.465231292995562 0.736117801507425 0.812565468088111;0.470256713212495 0.740210037864374 0.814205863148238;0.475365838917035 0.744293157911274 0.815851824819198;0.480562998989697 0.748366256556545 0.817501428549443;0.485851941201485 0.752428305219046 0.819154718051073;0.491234012361052 0.756478271585542 0.820816423078898;0.49671495679985 0.760515033473341 0.822482883601655;0.502295859707065 0.764537545992227 0.82415903136277;0.507979802498958 0.768544720246433 0.825845880691831;0.513770305653815 0.772535405623769 0.827543791816192;0.519667266878089 0.776508735656468 0.829258036256343;0.525674880370134 0.780463512637111 0.830987486471496;0.53179094151322 0.784399169529831 0.832738943598449;0.538018396345586 0.788314674334739 0.834512095957484;0.544353998015145 0.792209749803427 0.836313308498176;0.550797894410568 0.796083806195143 0.838144273543443;0.557346484750783 0.799936825697726 0.840009916922582;0.563996468518369 0.80376887330811 0.841913909985766;0.57074395045765 0.807580210811632 0.843859687136516;0.577582779315276 0.81137158033826 0.845851835241853;0.584508041943242 0.815143657217797 0.847893050463381;0.591512407998754 0.818897652532509 0.849987273018448;0.598589234360718 0.822634797950755 0.8521371651407;0.605731063578472 0.826356600007005 0.854345393801487;0.612931204139293 0.830064541934871 0.856613567178452;0.620181011883371 0.833760603385139 0.858944173773476;0.627475222688865 0.837446199935457 0.861337123590939;0.634805254687901 0.841123511153864 0.863794198852656;0.642166245076179 0.844794033685975 0.866314682226992;0.649551482953897 0.848459720845143 0.868898875084377;0.656955726394944 0.852122275510665 0.871546147479705;0.66437457944816 0.855783259371445 0.874255409158836;0.671802732456118 0.859444485917602 0.877026076134084;0.679237639981711 0.863107165518515 0.879856177129826;0.686674499674875 0.866773062939076 0.882744965490452;0.694111823878682 0.870443176568704 0.88569013517031;0.701545870701022 0.874119062064615 0.888690572364126;0.708975473275799 0.877801657414808 0.891744038328246;0.716398233094345 0.881492211368792 0.894848966559081;0.723812905154688 0.885191689231587 0.898003353967589;0.731218300317441 0.888901042777423 0.90120526120658;0.738613164399826 0.892621243218943 0.90445284662503;0.745996949314165 0.896353077537404 0.90774402851674;0.753368898657665 0.900097387269868 0.911076860417586;0.760728214258829 0.903855028013013 0.91444943824897;0.768075410485711 0.907626498534773 0.917859340770426;0.775408638046378 0.911412948937628 0.921305083288652;0.782729865690245 0.915214478005682 0.924783630754258;0.790036244556315 0.919032533279899 0.928293790777229;0.797330579828727 0.922866984587656 0.931832058593841;0.804609360361943 0.926719502515997 0.935397276108901;0.811876252186985 0.930589728810683 0.938985308363717;0.819126240304709 0.934479822131786 0.942595204174578;0.826364951124252 0.938388895029916 0.946221540017936;0.833585789978712 0.94231964805425 0.94986338237759;0.840792561737736 0.946271800064117 0.953515204467415;0.847982074911587 0.950247177694912 0.957173743064266;0.855150949964391 0.954247748207962 0.960835211869264;0.862302591933276 0.958273569685134 0.964492028423876;0.869429809242311 0.962327952845034 0.968140192127759;0.876528370393481 0.966413479837428 0.971773303620591;0.883593264528653 0.97053313520895 0.975383934307272;0.890617407410809 0.974690713609854 0.978964095113807;0.897591748513743 0.978890804105412 0.982505152803506;0.904501342713877 0.983139900522397 0.986000577365069;0.911330054230196 0.987444945995618 0.989444264262369;0.918059296008126 0.991813535883849 0.99283286383148],...
    'FontSize',14,'Layer','top','LineWidth',1);
% Create colorbar
colorbar(axes1,'Limits',[0 100]);

f = gcf;

if strcmp('frames',export)
    exportgraphics(f,strcat(output, 'LWP/', file, '_', nam.txt, '.png'),'Resolution',150)
else
    exportgraphics(f,strcat(output, 'LWP/', 'movie', '.gif'),'Resolution',150, 'Append',true)
end

end