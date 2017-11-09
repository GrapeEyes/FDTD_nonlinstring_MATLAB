%nlinstring_3DNonPhys

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Carlos de la Vega Martin
% MSc Acoustics and Music Technology at the University of Edinburgh.

% Notes: developed on MATLAB R2016b for the masters thesis.
% needs two external functions, gridInterpDirichlet.m,
% etafun.m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description
% This script models nonlinear string oscillation, as a
% system of 3 coupled pdes, for the longitudinal and
% the 2 orthogonal transversal oscillations respectively. Includes
% stiffnes and loss, x coordinate and displacements are
% normalized. Parameters can be set differently for the 2
% transversal polarizations.Grid spacing is different
% between the longintudinal and transversal polarizations,
% Interpolation is used between them

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Boundary conditions, stability & initialization
% Dirichlet boundary conditions at x=0, x=L. Clamped or
% simply supported boundary conditions can be chosen.
% We use the linear limit for the stability condition.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all

t1 = tic; % start timer

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% flags
nlin_on = 1; %linear or nonlinear
syssolv = 0; % method used for solving the linear system backlash = 0, Jacobi = 1, 2=both to compare
egy_on = 1; % calculate energy
ddtest_on = 0; % test that the system matrix is diagonally dominant
loss_fmt = 0; %0=lossless(or set the sigmas), 1=set the T60s
plot_on = 1; % plot spectra, and if calculated, energy
saveout_on = 0; % save output as .wav, variables as .mat, and if produced, plots as .fig and .eps

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parameters

% simpulation parameters
SR = 44100; % sampling rate
TF = 1; % duration of the simulation
eps = 1e-15; % convergence condition for the Jacobi method

% string parameters

ffrt_ver = 100; % fundamental frequency of transverse polarization
cents = 0;
ffrt_hor  = ffrt_ver*10^(cents/1200);
% ffrt_hor = ; % fundamental frequency of transverse polarization
ffrl = 2504; % fundamental frequency of longitudinal polarization
B_ver = 0; % inharmonicity
B_hor = B_ver;

% boundary conditions

%'clamped' , 'ssup' , 'clamped-ssup' , 'ssup-clamped'
bcond = 'clamped'; % for the Dxxxx operator
bconditp = 'revlin'; % for the interpolating matrix

% loss parameters

switch loss_fmt
    case 0 % loss coefficients directly
        
        sigmal0 = 0.0; % 1st order loss coeff, longitudinal pol
        sigmal1 = 0.0; % 3rd order loss coeff, longitudinal pol
        sigmat0_ver = 0.0; % 1st order loss coeff, transverse pol
        sigmat1_ver = 0.0; % 3rd order loss coeff, transverse pol
        sigmat0_hor = 0.0; % 1st order loss coeff, transverse pol
        sigmat1_hor = 0.0; % 3rd order loss coeff, transverse pol
        
    case 1 % loss fitted to T60 at two freq f1<f2, T601>T602
        fq1 = 100; % lower frequency for loss profile, Hz 
        fq2 = 2000; % higher frequency for loss profile, Hz
        
        T60l_fq1 = 50; % T60 longitudinal pol at freq 1, s
        T60l_fq2 = 30; % T60 longitudinal pol at freq 2, s
        T60t_ver_fq1 = 50; % T60 transverse pol at freq 1, s
        T60t_ver_fq2 = 30; % T60 tranverse pol at freq 2, s
        T60t_hor_fq1 = 50; % T60 transverse pol at freq 1, s
        T60t_hor_fq2 = 30; % T60 tranverse pol at freq 2, s
        
end


% I/O parameters

input_type = 'trinit'; % strike, sine, pluck
% for 0,1 an initial bandpassed noise can be added to any of the polarizations

switch input_type
    
    case 'strike' % raised cosine strike
        
        xh = 0.7; % excitation position
        Th = 0.003; % excitation time (s)
        fh = 50; % excitation strength (N/(Kg/m))
        angle = 0; % incidence angle, respect to the vertical (rad)
        
        u0t1 = 0.0; % maximum initial displacement, vert transv pol, rel to max disp in 2
        u0t2 = 0.0; % maximum initial displacement, horiz transv pol, rel to max disp in 1
        fcutoff_init = 1000; % cutoff frequency for bandpassing the noise term
        
    case 'sine' % sine input
        
        xh = 0.7; % excitation position
        Th = 0.05; % excitation ramp time
        fh = 5.0; % excitation strength (N/(Kg/m))
        inp_freq = 100; % frequency of the input sinusoid
        angle = 0; % incidence angle, respect to the vertical (rad)
        
        u0t1 = 0.0; % maximum initial displacement, vert transv pol, rel to max disp in 2
        u0t2 = 1; % maximum initial displacement, horiz transv pol, rel to max disp in 1
        fcutoff_init = 100; % cutoff frequency for bandpassing the noise term
        
    case 'pluck'
               
        xh = 0.7; % excitation position
        Th = 0.003; % excitation time (s)
        fh = 5.0; % excitation strength (N/(Kg/m))
        angle = 0; % incidence angle, respect to the vertical (rad)
        
        u0t1 = 0.0; % maximum initial displacement, vert transv pol, rel to max disp in 2
        u0t2 = 0.0; % maximum initial displacement, horiz transv pol, rel to max disp in 1
        fcutoff_init = 500; % cutoff frequency for bandpassing the noise term
        
    case 'rcinit'
        ctr = 0.6;
        wid = 0.2;
        u0v = 0.0001;
        
    case 'trinit'
        ctr = 0.5;
        u0v = 0.005;
        
    case 'randin'
        u0v = 1e-4;
        u0w = 0;
        
        
        
end


xo = 0.7; % output position (normalized, 0-1)

% plot parameters
maxfreq = 20000; % upper frequency limit for the spectrum plots
trans_time = 0.0; % for transverse movement plot, defines center of the transition phase, normalized respect to TF

% save parameters
SimID = 'E3_egyfix_Stiffterm96k_egy4'; % name of the folder, use an identifiable tag, like 'intrument_note'


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% derived parameters

% general simulation parameters
NF = floor(SR*TF); %number of samples of the simulation
k = 1/SR; % time step

% string parameters


gamma_ver = 2*ffrt_ver; % 2*fundamental of transverse mode
gamma_hor = 2*ffrt_hor; % 2*fundamental of transverse mode
alpha = 2*ffrl/gamma_ver; % 2*fundamental/gamma of longitudinal mode
kappa_ver = gamma_ver/pi*sqrt(B_ver); % stiffness coefficient
kappa_hor = gamma_hor/pi*sqrt(B_hor); % stiffness coefficient



% loss

if loss_fmt == 1
    eta1_ver = etafun(2*pi*fq1,gamma_ver,kappa_ver);
    eta2_ver = etafun(2*pi*fq2,gamma_ver,kappa_ver);
    eta1_hor = etafun(2*pi*fq1,gamma_hor,kappa_hor);
    eta2_hor = etafun(2*pi*fq2,gamma_hor,kappa_hor);
    
    sigmal0 = 6*log(10)/(eta2_ver-eta1_ver)*(eta2_ver/T60l_fq1 - eta1_ver/T60l_fq2); % 1st order loss coeff, longitudinal pol
    sigmal1 = 6*log(10)/(eta2_ver-eta1_ver)*(-1/T60l_fq1 + 1/T60l_fq2); % 3rd order loss coeff, longitudinal pol
    sigmat0_ver = 6*log(10)/(eta2_ver-eta1_ver)*(eta2_ver/T60t_ver_fq1 - eta1_ver/T60t_ver_fq2); % 1st order loss coeff, transverse pol
    sigmat1_ver = 6*log(10)/(eta2_ver-eta1_ver)*(-1/T60t_ver_fq1 + 1/T60t_ver_fq2);  % 3rd order loss coeff, transverse pol
    sigmat0_hor = 6*log(10)/(eta2_hor-eta1_hor)*(eta2_hor/T60t_hor_fq1 - eta1_ver/T60t_hor_fq2); % 1st order loss coeff, transverse pol
    sigmat1_hor = 6*log(10)/(eta2_hor-eta1_hor)*(-1/T60t_hor_fq1 + 1/T60t_hor_fq2);  % 3rd order loss coeff, transverse pol


end


% grid

hl = 2*ffrl*k; % minimum longitudinal grid spacing
ht_ver = sqrt(0.5*(gamma_ver^2*k^2+sqrt(gamma_ver^4*k^4+16*kappa_ver^2*k^2))); % minimum transversal grid spacing
ht_hor = sqrt(0.5*(gamma_hor^2*k^2+sqrt(gamma_hor^4*k^4+16*kappa_hor^2*k^2))); % minimum transversal grid spacing

ht = max([ht_ver,ht_hor]);

ht_ver = ht;
ht_hor = ht;

Nl = floor(1/hl); % number of segments, longitudinal grid
Nt_ver = floor(1/ht_ver); % number of segments, transversal grid
Nt_hor = floor(1/ht_hor); % number of segments, transversal grid

hl = 1/Nl; % longitudinal grid spacing
ht_ver = 1/Nt_ver; % transversal grid spacing
ht_hor = 1/Nt_hor; % transversal grid spacing

N = Nl+Nt_ver+Nt_hor-3; % size of the system matrix


% output position
lol = floor(xo*(Nl-2))+2; % grid index of output (longitudinal oscillation)
lot_ver = floor(xo*(Nt_ver-2))+2; % grid index of output (tranversal oscillation)
lot_hor = floor(xo*(Nt_hor-2))+2; % grid index of output (tranversal oscillation)



% finite difference matrix operators

% longitudinal grid
el = ones(Nl,1);
Dxlforw = spdiags([-el el], 0:1,Nl-1,Nl)/hl; %d/dx, forward difference
Dxlbackw = -Dxlforw'; % d/dx, backward difference
Dxxl = Dxlforw*Dxlbackw; % d2/dx2, centered difference

% tranversal grid operators, vertical
et_ver = ones(Nt_ver,1);
Dxtforw_ver = spdiags([-et_ver et_ver], 0:1,Nt_ver-1,Nt_ver)/ht_ver; % d/dx, forward difference
Dxtbackw_ver = -Dxtforw_ver'; % d/dx, backward difference
Dxxt_ver = Dxtforw_ver*Dxtbackw_ver; % d2/dx2, centered difference
Dxxxxt_ver = Dxxt_ver*Dxxt_ver; % d4/dx4, centered difference

% change values for u1, uN-1 depending on the boundary conditions
switch bcond
    case 'clamped'
        Dxxxxt_ver(1) = 7/ht_ver^4;
        Dxxxxt_ver(end) = 7/ht_ver^4;
    case 'clamped-ssup'
        Dxxxxt_ver(1) = 7/ht_ver^4;
    case 'ssup-clamped'
        Dxxxxt_ver(end) = 7/ht_ver^4;
end


% tranversal grid operators, horizontal
et_hor = ones(Nt_hor,1);
Dxtforw_hor = spdiags([-et_hor et_hor], 0:1,Nt_hor-1,Nt_hor)/ht_hor; % d/dx, forward difference
Dxtbackw_hor = -Dxtforw_hor'; % d/dx, backward difference
Dxxt_hor = Dxtforw_hor*Dxtbackw_hor; % d2/dx2, centered difference
Dxxxxt_hor = Dxxt_hor*Dxxt_hor; % d4/dx4, centered difference

% change values for u1, uN-1 depending on the boundary conditions
switch bcond
    case 'clamped'
        Dxxxxt_hor(1) = 7/ht_hor^4;
        Dxxxxt_hor(end) = 7/ht_hor^4;

    case 'clamped-ssup'
        Dxxxxt_hor(1) = 7/ht_hor^4;
    case 'ssup-clamped'
        Dxxxxt_hor(end) = 7/ht_hor^4;
end


% interpolation
Inpltv = gridInterpDirichlet(Nl,Nt_ver,'cubic', bconditp); % interpolation matrices
Inplth = Inpltv;
Inptvl = ht_ver/hl*Inpltv';
Inpthl = Inptvl;
[Inptvth,Inpthtv] = gridInterpDirichlet(Nt_ver,Nt_hor,'cubic', bconditp); % interpolation matrices


% create identity and zero matrices of required dimensions
It_ver = speye(Nt_ver-1); % identity matrix transversal grid
It_hor = speye(Nt_hor-1); % identity matrix transversal grid
Il = speye(Nl-1); % identity matrix longitudinal grid

Zv = sparse(Nl-1,Nt_ver-1); % zero matrix crossterms 1
Zvtr = Zv'; % zero matrix crossterms 2
Zh = sparse(Nl-1,Nt_hor-1); % zero matrix crossterms 1
Zhtr = Zh'; % zero matrix crossterms 2
Zvh = sparse(Nt_ver-1,Nt_hor-1); % zero matrix crossterms 1
Zvhtr = Zvh'; % zero matrix crossterms 2

Zt_ver = sparse(Nt_ver-1,Nt_ver-1); % zero matrix transversal grid
Zt_hor = sparse(Nt_hor-1,Nt_hor-1); % zero matrix transversal grid
Zl = sparse(Nl-1,Nl-1); % zero matrix longitudinal grid


% initialize diagonal matrices for the first derivative of
% the transverse polarization
Lambdav = sparse(Nt_ver,Nt_ver);
Lambdah = sparse(Nt_hor,Nt_hor);
Lambdav2 = Lambdav;
Lambdah2 = Lambdah;

% initialize system matrices off-diagonal blocks
Mnlin12 = Zv; Mnlin21 = Zvtr; Mnlin13 = Zh; Mnlin31 = Zhtr;
Mnlin22 = Zt_ver; Mnlin33 = Zt_hor;

% initialize outputs
out1 = zeros(NF,1); out2 = out1; out3 = out1;
H = out1;  itno = out1; T = out1; VLE = out1; VLS1=out1; VLS = out1; VNLV=out1; VNLH=out1; VNL = out1;
outv1 = out1; outv2 = out1; outv3 = out1;

% initialize system matrices
A = sparse(N,N);
B = A; C = A; D=A; Dinv = A; R=A; dd = zeros(N,1);

% vectors for linear indexing of the diagonal components
vecdv = [0:Nt_ver-1]*Nt_ver+[1:Nt_ver];
vecdh = [0:Nt_hor-1]*Nt_hor+[1:Nt_hor];
vecdfull = [0:(N-1)]*N+[1:N];

% input & initialization

fin = zeros(NF,1); % initialize vector of input (as long as output)
Jin = zeros(N,1); % initialize vector for force to displacement conversion


switch input_type
    
    case 'strike'
        
        fh1 = fh*cos(angle); % vertical maximum force
        fh2 = fh*sin(angle); % horizontal maximum force
        
        fdur_int = floor(Th*SR); %excitation duration in samples
        fin(1:fdur_int) = 0.5*(1-cos([0:fdur_int-1]'*2*pi/fdur_int)); % fill nonzero values of input
        
        xhintv = floor(xh*Nt_ver); % input position in samples
        xhinth = floor(xh*Nt_hor); % input position in samples
        Jin(Nl-1+xhintv) = fh1*k^2/ht_ver; % factor to convert from force to displacement, applied to only one grid point, vertical pol
        Jin(Nl+Nt_ver-2+xhinth) = fh2*k^2/ht_hor; % factor to convert from force to displacement, applied to only one grid point, horiz pol
        
    case 'sine'
        
        fh1 = fh*cos(angle); % vertical maximum force
        fh2 = fh*sin(angle); % horizontal maximum force
        
        fdur_int = floor(Th*SR); %excitation duration in samples
        fin(1:fdur_int) = 0.5*(1-cos([0:fdur_int-1]'*2*pi/fdur_int)); % fill nonzero values of input
        fin = fin.*sin(2*pi*inp_freq/SR*[1:NF]');
        
        xhintv = floor(xh*Nt_ver); % input position in samples
        xhinth = floor(xh*Nt_hor); % input position in samples
        Jin(Nl-1+xhintv) = fh1*k^2/ht_ver; % factor to convert from force to displacement, applied to only one grid point, vertical pol
        Jin(Nl+Nt_ver-2+xhinth) = fh2*k^2/ht_hor; % factor to convert from force to displacement, applied to only one grid point, horiz pol
        
    case 'pluck'
        
        fh1 = fh*cos(angle); % vertical maximum force
        fh2 = fh*sin(angle); % horizontal maximum force
        
        fdur_int = floor(Th*SR); %excitation duration in samples
        fin(1:fdur_int) = 0.5*(1-cos([0:fdur_int-1]'*pi/fdur_int)); % fill nonzero values of input
        
        xhintv = floor(xh*Nt_ver); % input position in samples
        xhinth = floor(xh*Nt_hor); % input position in samples
        Jin(Nl-1+xhintv) = fh1*k^2/ht_ver; % factor to convert from force to displacement, applied to only one grid point, vertical pol
        Jin(Nl+Nt_ver-2+xhinth) = fh2*k^2/ht_hor; % factor to convert from force to displacement, applied to only one grid point, horiz pol
end



% initialize string
m = zeros(N,1);
m_jac = m;
mdif = zeros(NF,1);
switch input_type
    case 'strike'
        bhi = fir1(406,2/SR*fcutoff_init,'low');
        rrv = (1-2*rand((Nt_ver-1),1));
        rrh = (1-2*rand((Nt_hor-1),1));
        rrfilv = filter(bhi,1,rrv);
        rrfilv = rrfilv/max(abs(rrfilv));
        rrfilh = filter(bhi,1,rrh);
        rrfilh = rrfilh/max(abs(rrfilh));
        m1 = [zeros((Nl-1),1);u0t1*fh2*k^2/ht_ver*rrfilv;u0t2*fh1*k^2/ht_hor*rrfilh];
        m2 = m1;
    case 'sine'
        bhi = fir1(406,2/SR*fcutoff_init,'low');
        rrv = (1-2*rand((Nt_ver-1),1));
        rrh = (1-2*rand((Nt_hor-1),1));
        rrfilv = filter(bhi,1,rrv);
        rrfilv = rrfilv/max(abs(rrfilv));
        rrfilh = filter(bhi,1,rrh);
        rrfilh = rrfilh/max(abs(rrfilh));
        m1 = [zeros((Nl-1),1);u0t1*fh2*k^2/ht_ver*rrfilv;u0t2*fh1*k^2/ht_hor*rrfilh];
        m2 = m1;
    case 'pluck'
        bhi = fir1(406,2/SR*fcutoff_init,'low');
        rrv = (1-2*rand((Nt_ver-1),1));
        rrh = (1-2*rand((Nt_hor-1),1));
        rrfilv = filter(bhi,1,rrv);
        rrfilv = rrfilv/max(abs(rrfilv));
        rrfilh = filter(bhi,1,rrh);
        rrfilh = rrfilh/max(abs(rrfilh));
        m1 = [zeros((Nl-1),1);u0t1*fh2*k^2/ht_ver*rrfilv;u0t2*fh1*k^2/ht_hor*rrfilh];
        m2 = m1;
        
        
    case 'rcinit'
        m1 = m;
        xax = [1:Nt_ver-1]'*ht_ver;
        ind = sign(max(-(xax-ctr-wid/2).*(xax-ctr+wid/2),0));
        rc = 0.5*ind.*(1+cos(2*pi*(xax-ctr)/wid));
        m1(Nl:Nl+Nt_ver-2) = u0v*rc;
        m2  = m1;
    case 'trinit'
        m1 = m;
        xax = [1:Nt_ver-1]'*ht_ver;
        tri = min(xax/ctr-1,0)+1+min((1-xax)/(1-ctr)-1,0);
        m1(Nl:Nl+Nt_ver-2) = u0v*tri;
        m2  = m1;
        
    case 'randin'
        m2 = m;
        m2(Nl:Nl+Nt_ver-2) = u0v*sin(pi*[1:Nt_ver-1]/(Nt_ver+1));
        m2(Nl+Nt_ver-1:end) = u0w*(1-2*rand(1,Nt_ver-1)).*sin(pi*[1:Nt_hor-1]/(Nt_hor+1));
        m1 = m2;
        
end
%m1 = m; m2 = m;

% linear components of the system
A00 = sparse([(1+k*sigmal0)*Il-k*sigmal1*Dxxl, Zv, Zh; Zvtr, (1+k*sigmat0_ver)*It_ver-k*sigmat1_ver*Dxxt_ver, Zvh; Zhtr, Zvhtr, (1+k*sigmat0_hor)*It_hor-k*sigmat1_hor*Dxxt_hor]);
B00 = sparse([-2*Il-(k*2*ffrl)^2*Dxxl, Zv, Zh; Zvtr, -2*It_ver-(k*gamma_ver)^2*Dxxt_ver+(k*kappa_ver)^2*Dxxxxt_ver, Zvh; Zhtr, Zvhtr, -2*It_hor-(k*gamma_hor)^2*Dxxt_hor+(k*kappa_hor)^2*Dxxxxt_hor]);
C00 = sparse([(1-k*sigmal0)*Il+k*sigmal1*Dxxl, Zv, Zh; Zvtr, (1-k*sigmat0_ver)*It_ver+k*sigmat1_ver*Dxxt_ver, Zvh; Zhtr, Zvhtr, (1-k*sigmat0_hor)*It_hor+k*sigmat1_hor*Dxxt_hor]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
phisq_ver = k^2*(ffrl^2-ffrt_ver^2); % coupling coefficient
phisq_hor = k^2*(ffrl^2-ffrt_hor^2); % coupling coefficient


%auxiliary matrices
auxmat1_ver = phisq_ver*Inptvl*Dxtforw_ver;
auxmat1_hor = phisq_hor*Inpthl*Dxtforw_hor;
auxmat2_ver = phisq_ver*Dxtforw_ver;
auxmat2_hor = phisq_hor*Dxtforw_hor;
auxmat3_ver = Dxtbackw_ver*Inpltv;
auxmat3_hor = Dxtbackw_hor*Inplth;
auxmat4_ver = Inpthtv*auxmat2_hor;
auxmat4_hor = Inptvth*auxmat2_ver;
auxmat5_ver = Dxtbackw_hor*Inptvth;
auxmat5_hor = Dxtbackw_ver*Inpthtv;

% linear

if nlin_on==0
    A = A00;
    B = B00;
    C = C00;
end


for n = 1:NF
    
    if nlin_on == 1
        % first derivative vector transverse polarizations
        lambdavecv = Dxtbackw_ver*m1(Nl:Nl+Nt_ver-2);
        lambdavech = Dxtbackw_hor*m1(Nl+Nt_ver-1:end);
        Lambdav(vecdv) = lambdavecv;
        Lambdah(vecdh) = lambdavech;
        Lambdav2(vecdv) = lambdavecv.^2;
        Lambdah2(vecdh) = lambdavech.^2;
        
        % non-linear matrix terms
        Mnlin12 = auxmat1_ver*Lambdav*Dxtbackw_ver;
        Mnlin21 = auxmat2_ver*Lambdav*auxmat3_ver;
        Mnlin13 = auxmat1_hor*Lambdah*Dxtbackw_hor;
        Mnlin31 = auxmat2_hor*Lambdah*auxmat3_hor;
        Mnlin22 = auxmat2_ver*(Lambdav2)*Dxtbackw_ver + auxmat4_ver*(Lambdah2)*auxmat5_ver;
        Mnlin33 = auxmat2_hor*(Lambdah2)*Dxtbackw_hor + auxmat4_hor*(Lambdav2)*auxmat5_hor;
        M00 = sparse([Zl, -Mnlin12, -Mnlin13; -Mnlin21, -Mnlin22, Zvh; -Mnlin31, Zvhtr, -Mnlin33]);
        
        % build full system matrices
        A = A00 + M00;
        B = B00 + sparse([Zl, Zv, Zh; -2*Mnlin21, Zt_ver, Zvh; -2*Mnlin31, Zvhtr, Zt_hor]);
        C = C00 + M00;
    end
    
    
    if ddtest_on == 1 % convergence test for Jacobi, sufficient but not necessary condition
        dd = A(vecdfull); % vector  of diagonal components
        D(vecdfull) = dd; % diagonal terms matrix
        R = A-D; % non-diagonal terms matrix
        if sum([abs(dd)]'>sum(abs(R),2))/N < 1 % check A is diagonally dominant
            sprintf('system matrix is not diagonally dominant, convergence not guaranteed')
            sprintf('simulation stopped at time = \n %1e', k*n)
            return
        end
    end
    
    % solve system fo current step
    if syssolv == 0 % solve using MATLAB backlash operation
        m = A\(-B*m1-C*m2+Jin*fin(n));
        
    elseif syssolv == 1 % Solve using Jacobi method
        count = 0;
        mk1 = m1; % last step state vector
        mk = ones(Nl+Nt_ver+Nt_hor-3,1);
        b = -B*m1-C*m2+Jin*fin(n);
        dd = A(vecdfull); % extract diagonal terms
        D(vecdfull) = dd; % diagonal terms matrix
        Dinv(vecdfull) = 1./dd; % inverse
        R = A-D; % non-diagonal terms matrix
        while sqrt(sum((A*mk-b).^2)/(sum(mk.^2)))>eps
            mk = Dinv*(b-R*mk1); % calculate new state vector
            mk1 = mk; % cycle state vector
            count = count+1;
            if count > 10000
                message = sprintf('too many iterations\n simulation stopped at time = \n %1e', k*n);
                
                return
            end
        end
        
        itno(n) = count; % iterations counter vector
        m = mk; %
        
    elseif syssolv == 2 % use both, keep backlash, calculate difference
        count = 0;
        mk1 = m1; % last step state vector
        mk = ones(Nl+Nt_ver+Nt_hor-3,1);
        b = -B*m1-C*m2+Jin*fin(n);
        dd = A(vecdfull); % extract diagonal terms
        D(vecdfull) = dd; % diagonal terms matrix
        Dinv(vecdfull) = 1./dd; % inverse
        R = A-D; % non-diagonal terms matrix
        while sqrt(sum((A*mk-b).^2)/(sum(mk.^2)))>eps
            mk = Dinv*(b-R*mk1); % calculate new state vector
            mk1 = mk; % cycle state vector
            count = count+1;
            if count > 10000
                message = sprintf('too many iterations\n simulation stopped at time = \n %1e', k*n);
                
                return
            end
        end
        
        itno(n) = count; % iterations counter vector
        m_jac = mk; %
        
        m_bl = A\(-B*m1-C*m2+Jin*fin(n));
        mdif(n) = sqrt(sum((m_jac-m_bl).^2));
        m = m_bl;
    end
    
    % write ouput positions
    out1(n) = m(lol);
    out2(n) = m(Nl-1 +lot_ver);
    out3(n) = m(Nl+Nt_ver-2 +lot_hor);
    
%     if mod(n,3) == 0
%         plot(k*[1:n],out3(1:n),'k',k*[1:n],out2(1:n),'g');
%         
%         drawnow
%     end
    

    
% write output velocities
    outv1(n) = (m(lol)-m1(lol))./k;
    outv2(n) = (m(Nl-1 +lot_ver)-m1(Nl-1 +lot_ver))./k;
    outv3(n) = (m(Nl+Nt_ver-2 +lot_hor)-m1(Nl+Nt_ver-2 +lot_hor))./k;
    
    
    
    
    
    if(egy_on==1) % energy calculation
        switch bcond
            case 'clamped'
                bc0 = 0;
                bcN = 0;
            case 'ssup'
                bc0 = 1;
                bcN = 1;
            case 'clamped-ssup'
                bc0 = 0;
                bcN = 1;
            case 'ssup-clamped'
                bc0 = 1;
                bcN = 0;
        end
        u = m(1:Nl-1);
        u1 = m1(1:Nl-1);
        v = m(Nl:Nl+Nt_ver-2);
        v1 = m1(Nl:Nl+Nt_ver-2);
        v2 = m2(Nl:Nl+Nt_ver-2);
        w = m(Nl+Nt_ver-1:end);
        w1 = m1(Nl+Nt_ver-1:end);
        
        % derivatives
        p = Dxtbackw_ver*Inpltv*m(1:Nl-1); % du/dx, current step
        p1 = Dxtbackw_ver*Inpltv*m1(1:Nl-1); % du/dx, previous step
        q = Dxtbackw_ver*m(Nl:Nl+Nt_ver-2); % dv/dx, current step
        q1 = Dxtbackw_ver*m1(Nl:Nl+Nt_ver-2); % dv/dx, previous step
        q2 = Dxtbackw_ver*m2(Nl:Nl+Nt_ver-2); % dv/dx, previous previous step
        qq = Dxxt_ver*v;
        qq1 = Dxxt_ver*v1;
        r = Dxtbackw_ver*Inpthtv*m(Nl+Nt_ver-1:end); % dw/dx, current step
        r1 = Dxtbackw_ver*Inpthtv*m1(Nl+Nt_ver-1:end); % dw/dx, previous step
        
        T(n) = 0.5/k^2*([hl*el(1:end-1)',ht_ver*et_ver(1:end-1)',ht_hor*et_hor(1:end-1)']*(m-m1).^2); % kinetic energy
        VLE(n) = 0.5*gamma_ver^2*ht_ver*(sum(q.*q1)) + 0.5*gamma_hor^2*ht_ver*sum(r.*r1) + 2*ffrl^2*ht_ver*sum(p.*p1); % potential energy due to elastic terms
        VLS(n) = 0.5*kappa_ver^2*ht_ver*(sum(qq.*qq1) + 1/ht_ver^4*(2*(1-bc0)*v(1)*v1(1)+2*(1-bcN)*v(end)*v1(end)))...
            +0.5*kappa_hor^2*ht_hor*(sum((Dxxt_hor*w).*(Dxxt_hor*w1)) + 1/ht_hor^4*(2*(1-bc0)*w(1)*w1(1)+2*(1-bcN)*w(end)*w1(end))); % potential energy due to stiffness terms
        VLS1(n) = 0.5*kappa_ver^2*ht_ver*(sum(qq.*qq1))+0.5*kappa_hor^2*ht_hor*(sum((Dxxt_hor*w).*(Dxxt_hor*w1)));% potential energy due to stiffness terms
        
        if nlin_on == 1
            VNL(n) = 0.5*(4*ffrl^2-gamma_ver^2)*ht_ver*(0.5*sum((p+p1).*(q.*q1+r.*r1))+0.25*sum((q.*q1+r.*r1).^2)+0.25*sum((q.*r1-r.*q1).^2));
        end
        H(n) = T(n)+VLE(n)+VLS(n)+VNL(n);
        
        Hnorm = (H-H(1))/max(H);
        if mod(n,10) == 0
            figure(1)
            plot(ht_ver*[1:Nt_ver-1]',v)
            title('transversal displacement')
            drawnow
            figure(2)
            plot(k*[1:n],Hnorm(1:n),'r')
            title('relative error in the energy')
            %k*[1:n],VLS(1:n),'r',,k*[1:n],T(1:n)+VLE(1:n),'g'
            drawnow
        end

    end
    
    % cycle stored string states
    m2 = m1;
    m1 = m;
    
    
end

runtime = toc(t1) % get scheme runtime
H = T+VLE+VLS+VNL;

out_pos = out1+out2+out3; % combined position output
out_vel = outv1+outv2+outv3; % combined velocity output

if egy_on==1
    s = {'strike','pluck','sine'};
    if sum(strcmp(input_type,s))>0
        Hvar = max(H(2*floor(Th*SR):end))-min(H(2*floor(Th*SR):end));
        Hvarnorm = Hvar/max(H(2*floor(Th*SR):end));
        Hinloss = max(H)-max(H(2*floor(Th*SR):end));
    else
        Hvar = max(H)-min(H);
        Hvarnorm = Hvar/max(H);
        Hnorm = (H-H(1))/max(H);
    end
end


if plot_on == 1
    
    % spectrum parameters
    nfft = 2^nextpow2(NF);
    df = SR/nfft;
    freqv = 0:df:maxfreq;
    nbins  = numel(freqv);
    win_han = [(sin(pi*[0:NF-1]/(NF-1))).^2]';
    %win_han = ones(NF,1);
    magf1 = abs(fft(win_han.*outv1,nfft));
    magf2 = abs(fft(win_han.*outv2,nfft));
    magf3 = abs(fft(win_han.*outv3,nfft));
    maxf1 = max(magf1);
    maxf2 = max(magf2);
    maxf3 = max(magf3);
    fftmax = max([maxf1,maxf2,maxf3]);
    
    fig1 = figure(1); % spectrum plots for the three polarizations, separated
    
    subplot(3,1,1)
    vfft1 = 20*log10(magf1/maxf1);
    plot(freqv,vfft1(1:nbins),'k');
    axis tight
    xlabel('frequency (Hz)');
    ylabel('amplitude (dB)');
    title('longitudinal polarization spectrum')
    
    subplot(3,1,2)
    vfft2 = 20*log10(magf2/maxf2);
    plot(freqv,vfft2(1:nbins),'k');
    axis tight
    xlabel('frequency (Hz)');
    ylabel('amplitude (dB)');
    title('vertical polarization spectrum');
    
    subplot(3,1,3)
    vfft3 = 20*log10(magf3/maxf3);
    plot(freqv,vfft3(1:nbins),'k');
    axis tight
    xlabel('frequency (Hz)');
    ylabel('amplitude (dB)');
    title('horizontal polarization spectrum');
    
    
    fig2 = figure(2); % spectrum plot of the
    vfft1_all = 20*log10(magf1/fftmax);
    vfft2_all = 20*log10(magf2/fftmax);
    vfft3_all = 20*log10(magf3/fftmax);
    plot(freqv,vfft1_all(1:nbins),'k',freqv,vfft2_all(1:nbins),'b',freqv,vfft3_all(1:nbins),'g');
    axis tight
    xlabel('frequency (Hz)');
    ylabel('amplitude (dB)');
    ylim([-100, 0]);
    title('output spectra of the three polarizations');
    legend('longitudinal','vertical','horizontal');
    
    
    fig3 = figure(3);
    ram1 = linspace(0,1,NF/2);
    ram2 = flip(ram1);
    vr = [zeros(1,NF/2),ram1]';
    vg = [ram1,ram2]';
    vb = [ram2,zeros(1,NF/2)]';
    c = [vr,vg,vb];
    scatter(out3,out2,1,c)
    axis tight
    xlabel('horizontal position (m)')
    ylabel('vertical position (m)')
    title('string transverse polarization')
    
    
    if egy_on ==1
        fig4 = figure(4);
        plot(k*[1:NF],Hnorm,'k');
        axis tight
        xlabel('time (seconds)');
        ylabel('energy (J/(Kg/m))');
        title('discrete energy');

    end
    
end % plotting

% save data
if saveout_on == 1
    
    %dirname = strcat(cd,'\\data\\',SimID,sprintf('_fh-%1$2.4g_ang-%2$2.4g_Th-%3$2.4g_',fh,angle,Th),input_type,'_',bcond);
    dirname = strcat(cd,'\\data\\',SimID,'\\',input_type,'_',bcond,'_',sprintf('B_ver-%1$2.4g_u0v-%2$2.4g',B_ver,u0v));
    mkdir(dirname);
    if plot_on == 1
        savefig(fig1,strcat(dirname,'\\spectra.fig'));
        saveas(fig1,strcat(dirname,'\\spectra.eps'),'epsc');
        
        savefig(fig2,strcat(dirname,'\\spectra1pl.fig'));
        saveas(fig2,strcat(dirname,'\\spectra1pl.eps'),'epsc');
        
        savefig(fig3,strcat(dirname,'\\trans_pos.fig'));
        saveas(fig3,strcat(dirname,'\\trans_pos.eps'),'epsc');
        
        if egy_on == 1
            savefig(fig4,strcat(dirname,'\\egyplot.fig'))
            saveas(fig4,strcat(dirname,'\\egyplot.eps'),'epsc')
            savefig(fig5,strcat(dirname,'\\egyFFTdifplot.fig'))
            saveas(fig5,strcat(dirname,'\\egyFFTdifplot.eps'),'epsc')
        end
    end
    out_arr = [out1,out2,out3];
    outv_arr = [outv1,outv2,outv3];
    maxtot = max(max(outv_arr));
    save(strcat(dirname,'\\wksp.mat'));
    audiowrite(strcat(dirname,'\\audio_ver.wav'),0.8*outv2/maxtot,SR);
    audiowrite(strcat(dirname,'\\audio_hor.wav'),0.8*outv3/maxtot,SR);
    audiowrite(strcat(dirname,'\\audio_lon.wav'),0.8*outv1/maxtot,SR);
    audiowrite(strcat(dirname,'\\audio_norm.wav'),0.8*out_vel/max(abs(out_vel)),SR);
    
    
end