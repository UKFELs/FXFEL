clear all;

N = 1000;



phscons = getPhysConsts(); % physical constants
pms = getParams();

sigmax = 5.353734158843998e-04;
sigmay = sigmax;

sigmaxd = 5.946503702917309e-04;
sigmayd = sigmaxd;

sigmag = 1e-3;

sigmat = 6.9208e-13;

sigz2 = 13 - (7.5*0.4);
sigmaedge = 0.4;
sigxbar = 0.33635981;





%%%%%%%%%%%%%%%%%%%%%   Charge in center (flat-top)   %%%%%%%%%%%%%%%%%%




Ncenter = pms.npk_bar * sigz2 * 2*pi*sigxbar^2;
Qcenter = Ncenter * phscons.q_e;





%%%%%%%%%%%%%%%%%%%%%   Charge in tails   %%%%%%%%%%%%%%%%%%

Nedge = pms.npk_bar * (sqrt(2*pi)*sigmaedge) * 2*pi*sigxbar^2;
Qedge = Nedge * phscons.q_e;





NMcenter = round(N * Qcenter /(Qcenter + Qedge));
NMedge = N - NMcenter;



%%%%%%%%%%%%%%%%%%%   Generate and plot distributions   %%%%%%%%%%%%%

x0 = randn(1,NMcenter)*sigmax;
y0 = randn(1,NMcenter)*sigmay;
xd0 = randn(1,NMcenter)*sigmaxd;
yd0 = randn(1,NMcenter)*sigmayd;
E0 = (randn(1,NMcenter)*sigmag*pms.gamma_r + pms.gamma_r) * phscons.m_e ...
                        * phscons.c^2 / phscons.q_e;
s0 = rand(1,NMcenter)*sigmat;
figure; plot(s0,E0,'.');



xe = randn(1,NMedge)*sigmax;
ye = randn(1,NMedge)*sigmay;
xde = randn(1,NMedge)*sigmaxd;
yde = randn(1,NMedge)*sigmayd;
Ee = (randn(1,NMedge)*sigmag*pms.gamma_r + pms.gamma_r) * phscons.m_e ...
                        * phscons.c^2 / phscons.q_e;
se = randn(1,NMedge)*sigmaedge * pms.Lc / phscons.c;

se(se>0) = se(se>0) + sigmat;

figure; plot(se,Ee,'.');


x0 = [x0, xe];
y0 = [y0, ye];
xd0 = [xd0, xde];
yd0 = [yd0, yde];
E0 = [E0, Ee];
s0 = [s0, se];

figure; plot(s0,E0,'.');


%%%%%%%%%    Charge on each macroparticle is:


QMP = 4.518997285641631E-009 / N



