clear all;

N = 10000;


h = 6.626e-34; % Planck constant
q_e = 1.60217646e-19; % Charge on electron
c = 2.99792458e8; % Speed of light in vacuum
m_e = 9.109e-31;  % Electron rest mass 
eps_0 = 8.854188e-12; % Permittivity of free space

Lc = 1.591788158834600e-05;
Lg = 0.159154938663010;
gamma_r = 100;
aw = 1;
npk_bar = 3.596654943274300e+09;



%phscons = getPhysConsts(); % physical constants
%pms = getParams();








sigmax = 5.353734158843998e-04;
sigmay = sigmax;

sigmaxd = 5.946503702917309e-04;
sigmayd = sigmaxd;

sigmag = 1e-3;

sigmat = 5.309633769487956e-13;

sigz2 = 13 - (7.5*0.4);
sigmaedge = 0.4;
sigxbar = 0.33635981;





%%%%%%%%%%%%%%%%%%%%%   Charge in center (flat-top)   %%%%%%%%%%%%%%%%%%




Ncenter = npk_bar * sigz2 * 2*pi*sigxbar^2;
Qcenter = Ncenter * q_e;





%%%%%%%%%%%%%%%%%%%%%   Charge in tails   %%%%%%%%%%%%%%%%%%

Nedge = npk_bar * (sqrt(2*pi)*sigmaedge) * 2*pi*sigxbar^2;
Qedge = Nedge * q_e;





NMcenter = round(N * Qcenter /(Qcenter + Qedge));
NMedge = N - NMcenter;



%%%%%%%%%%%%%%%%%%%   Generate and plot distributions   %%%%%%%%%%%%%

x0 = randn(1,NMcenter)*sigmax;
y0 = randn(1,NMcenter)*sigmay;
xd0 = randn(1,NMcenter)*sigmaxd;
yd0 = randn(1,NMcenter)*sigmayd;
E0 = (randn(1,NMcenter)*sigmag*gamma_r + gamma_r) * m_e ...
                        * c^2 / q_e;
s0 = rand(1,NMcenter)*sigmat;
%figure; plot(s0,E0,'.');



xe = randn(1,NMedge)*sigmax;
ye = randn(1,NMedge)*sigmay;
xde = randn(1,NMedge)*sigmaxd;
yde = randn(1,NMedge)*sigmayd;
Ee = (randn(1,NMedge)*sigmag*gamma_r + gamma_r) * m_e ...
                        * c^2 / q_e;
se = randn(1,NMedge)*sigmaedge * Lc / c;

se(se>0) = se(se>0) + sigmat;

%figure; plot(se,Ee,'.');


x0 = [x0, xe];
y0 = [y0, ye];
xd0 = [xd0, xde];
yd0 = [yd0, yde];
E0 = [E0, Ee];
s0 = [s0, se];

figure; plot(s0,E0,'.');
xlabel('time (s)');
ylabel('Energy (eV)');


%%%%%%%%%    Charge on each macroparticle is:


QMP = 4.518997285641631E-009 / N
Q0 = zeros(1,N);
Q0(:) = QMP;

%%% x0 is x coordinates in metres
%%% y0 is x coordinates in metres
%%% xd0 is dx/dz
%%% yd0 is dy/dz
%%% E0 is energy in eV
%%% s0 is temporal coordinate in seconds
%%% Q0 is charge in macroparticle in Coulombs


dlmwrite('test.txt',[x0', y0', xd0', yd0', E0', s0', Q0'], 'delimiter',',');





