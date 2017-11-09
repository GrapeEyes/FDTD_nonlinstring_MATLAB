% nlinstring_3DPhysical

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Carlos de la Vega Martin, 19/08/2017
% MSc Acoustics and Music Technology at the University of Edinburgh.

% Notes: developed on MATLAB R2016b for the masters thesis.
% needs one external function etafun.m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Description
% This script models nonlinear string oscillation, as a
% system of 3 coupled pdes, for the longitudinal and
% the 2 orthogonal transversal oscillations respectively. Includes
% stiffnes and loss, x coordinate and displacements are
% normalized. Loss can be set different for each
% polarization.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Boundary conditions, stability & initialization
% Dirichlet boundary conditions at x=0, x=L. Clamped or
% simply supported boundary conditions can be chosen.
% We use the linear limit for the stability condition.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%clear all
close all

t1 = tic; % start timer

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% flags

syssolv = 0; % method used for solving the linear system backlash = 0, Jacobi = 1
egy_on = 0; % calculate energy
ddtest_on = 0; % test that the system matrix is diagonally dominant
plot_on = 1; % plot spectra, and if calculated, energy
saveout_on = 0; % save output as .wav, variables as .mat, and if produced, plots as .fig and .eps
str_param_fmt = 1; % physical=0, fundamentals+inharm=1  eq coeff=2
loss_fmt = 1; % sigmas=0, T60+freq=1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% parameters

% simpulation parameters
SR = 192000; % sampling rate
TF = 7; % duration of the simulation
eps = 1e-15; % convergence condition for the Jacobi method

% string parameters

switch str_param_fmt
    case 0 % physical parameters
        
        T0 = 100; % string tension (N)
        E = 2e11; % Young modulus (Pa/m)
        rho = 7850; % density (kg/m^3)
        r0 = 0.4e-3; % string radius (m)
        reff = 0.2e-3; % string effective radius (m), for the stiffness component
        L = 1.3; % length (m)
        
    case 1 % polarization fundamentals and inharmonicity
        
        ffrt = 164.814; % fundamental frequency of transverse polarization
        ffrl = 2430; % fundamental frequency of longitudinal polarization
        Binh = 1.51e-4; % inharmonicity
        
    case 2 % scheme parameters directly
        
        gamma = 200; % 2*fundamental of transverse mode
        alpha = 15; % 2*fundamental/gamma of longitudinal mode
        kappa = 1.5; % stiffness coefficient
end

% boundary conditions

%'clamped' , 'ssup' , 'clamped-ssup' , 'ssup-clamped'
bcond = 'ssup'; % for the Dxxxx operator

% loss parameters

switch loss_fmt
    case 0 % loss coefficients directly
        
        sigmal0 = 0.0; % 1st order loss coeff, longitudinal pol
        sigmal1 = 0.0; % 3rd order loss coeff, longitudinal pol
        sigmatv0 = 0.0; % 1st order loss coeff, transverse pol, vert
        sigmatv1 = 0.0; % 3rd order loss coeff, transverse pol, vert
        sigmath0 = 0.0; % 1st order loss coeff, transverse pol, hor
        sigmath1 = 0.0; % 3rd order loss coeff, transverse pol, hor
        
    case 1 % loss fitted to T60 at two freq f1<f2, T601>T602
        fq1 = 100; % lower frequency for loss profile, Hz 
        fq2 = 2000; % higher frequency for loss profile, Hz
        
        T60l_fq1 = 13; % T60 longitudinal pol at freq 1, s
        T60l_fq2 = 10; % T60 longitudinal pol at freq 2, s
        T60tv_fq1 = 11; % T60 transverse pol v at freq 1, s
        T60tv_fq2 = 8; % T60 tranverse pol v at freq 2, s
        T60th_fq1 = 11; % T60 transverse pol h at freq 1, s
        T60th_fq2 = 8; % T60 tranverse pol h at freq 2, s
        
end


% I/O parameters

input_type = 'strike'; % strike, sine, pluck
% for 0,1 an initial bandpassed noise can be added to any of the polarizations

switch input_type
    
    case 'strike' % raised cosine strike
        
        xh = 0.7; % excitation position
        Th = 0.015; % excitation time (s)
        fh = 100.0; % excitation strength (N/(Kg/m))
        angle = pi/12; % incidence angle, respect to the vertical (rad)
        
        u0t1 = 1e-8; % maximum initial displacement, vert transv pol, rel to max disp in 2
        u0t2 = 1e-8; % maximum initial displacement, horiz transv pol, rel to max disp in 1
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
        fh = 10.0; % excitation strength (N/(Kg/m))
        angle = 0.0; % incidence angle, respect to the vertical (rad)
        
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
        u0tot = 1e-2;
        
        fh = 0;
        xh = 0;
        Th = 0;
        angle = 0;
        u0v = 1e-2;
        u0w = 1e-10;
%         u0v = u0tot*cos(angle);
%         u0w = u0tot*sin(angle);
end


xo = 0.4; % output position (normalized, 0-1)

% plot parameters
maxfreq = 15000; % upper frequency limit for the spectrum plots
trans_time = 0.0; % for transverse movement plot, defines center of the transition phase, normalized respect to TF

% save parameters
SimID = 'E3_noitp2'; % name of the folder, use an identifiable tag, like 'intrument_note'


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% derived parameters

% general simulation parameters
NF = floor(SR*TF); %number of samples of the simulation
k = 1/SR; % time step

% string parameters

if str_param_fmt == 0 % physical parameters given

    A0 = pi*r0^2; % string crosssection (m^2)
    I0 = 0.25*pi*reff^4; % area moment of inertia
    gamma = sqrt(T0/(rho*A0))/L; % elastic coefficient
    alpha = sqrt(E*A0/T0); % adimensional coefficient
    kappa = sqrt(E*I0/rho/A0)/(L^2); % stiffness coefficient

elseif str_param_fmt == 1 % fundamentals and inharmonicity given

    gamma = 2*ffrt; % 2*fundamental of transverse mode
    alpha = 2*ffrl/gamma; % 2*fundamental/gamma of longitudinal mode
    kappa = gamma/pi*sqrt(Binh); % stiffness coefficient
    
end
    

% loss

if loss_fmt == 1
    eta1 = etafun(2*pi*fq1,gamma,kappa);
    eta2 = etafun(2*pi*fq2,gamma,kappa);
    
    sigmal0 = 6*log(10)/(eta2-eta1)*(eta2/T60l_fq1 - eta1/T60l_fq2); % 1st order loss coeff, longitudinal pol
    sigmal1 = 6*log(10)/(eta2-eta1)*(-1/T60l_fq1 + 1/T60l_fq2); % 3rd order loss coeff, longitudinal pol
    sigmatv0 = 6*log(10)/(eta2-eta1)*(eta2/T60tv_fq1 - eta1/T60tv_fq2); % 1st order loss coeff, transverse pol v
    sigmatv1 = 6*log(10)/(eta2-eta1)*(-1/T60tv_fq1 + 1/T60tv_fq2);  % 3rd order loss coeff, transverse pol v
    sigmath0 = 6*log(10)/(eta2-eta1)*(eta2/T60th_fq1 - eta1/T60th_fq2); % 1st order loss coeff, transverse pol h
    sigmath1 = 6*log(10)/(eta2-eta1)*(-1/T60th_fq1 + 1/T60th_fq2);  % 3rd order loss coeff, transverse pol h

end


    % grid

    ht = sqrt(0.5*(gamma^2*k^2+sqrt(gamma^4*k^4+16*kappa^2*k^2))); % minimum transversal grid spacing
    hl = gamma*alpha*k; % minimum longitudinal grid spacing
    
    h = max([hl,ht]);
    
    N = floor(1/h); % number of segments
    h = 1/N; % grid spacing
    Ntot = 3*N-3; % size of the system matrix


    % output position
    lo = floor(xo*(N-2))+2; % grid index of output (tranversal oscillation)


    % finite difference matrix operators

    % tranversal grid operators
    e1 = ones(N,1);
    Dxf = spdiags([-e1 e1], 0:1,N-1,N)/h; % d/dx, forward difference
    Dxb = -Dxf'; % d/dx, backward difference
    Dxx = Dxf*Dxb; % d2/dx2, centered difference
    Dxxxx = Dxx*Dxx; % d4/dx4, centered difference

    % change values for u1, uN-1 depending on the boundary conditions
    switch bcond
        case 'clamped'
            Dxxxx(1) = 7/h^4;
            Dxxxx(end) = 7/h^4;
        case 'clamped-ssup'
            Dxxxx(1) = 7/h^4;
        case 'ssup-clamped'
            Dxxxx(end) = 7/h^4;
    end


    
    % create identity and zero matrices of required dimensions
    I = speye(N-1); % identity matrix
    Z = sparse(N-1,N-1); % zero matrix


    % initialize diagonal matrices for the first derivative of
    % the transverse polarization
    Lambdav = sparse(N,N); Lambdah = Lambdav; 


    % initialize outputs
    out1 = zeros(NF,1); out2 = out1; out3 = out1;
    H = out1;  itno = out1;
    outv1 = out1; outv2 = out1; outv3 = out1;

    % initialize system matrices
    A = sparse(Ntot,Ntot);
    B = A; C = A; D=A; Dinv = A; R=A; dd = zeros(Ntot,1);

    % vectors for linear indexing of the diagonal components
    vecd = [0:N-1]*N+[1:N];
    vecdfull = [0:(Ntot-1)]*Ntot+[1:Ntot];



    % input & initialization
    fh1 = fh*cos(angle); % vertical maximum force
    fh2 = fh*sin(angle); % horizontal maximum force
    fin = zeros(NF,1); % initialize vector of input (as long as output)


switch input_type
    case 'strike'
        fdur_int = floor(Th*SR); %excitation duration in samples
        fin(1:fdur_int) = 0.5*(1-cos([0:fdur_int-1]'*2*pi/fdur_int)); % fill nonzero values of input
    case 'sine'
        fdur_int = floor(Th*SR); %excitation duration in samples
        fin(1:fdur_int) = 0.5*(1-cos([0:fdur_int-1]'*2*pi/fdur_int)); % fill nonzero values of input
        fin = fin.*sin(2*pi*inp_freq/SR*[1:NF]');
    case 'pluck'
        fdur_int = floor(Th*SR); %excitation duration in samples
        fin(1:fdur_int) = 0.5*(1-cos([0:fdur_int-1]'*pi/fdur_int)); % fill nonzero values of input
end


Jin = zeros(Ntot,1); % initialize vector for force to displacement conversion
lh = floor(xh*(N-2))+2; % input position in samples
Jin(N-1+lh) = fh1*k^2/h; % factor to convert from force to displacement, applied to only one grid point, vertical pol
Jin(2*N-2+lh) = fh2*k^2/h; % factor to convert from force to displacement, applied to only one grid point, horiz pol

% initialize string
m = zeros(Ntot,1);
switch input_type
    case 'strike'
        lhi = fir1(106,2/SR*[10 fcutoff_init]);
        rr = (1-2*rand((N-1),1));
        rrfil = filter(lhi,1,rr);
        rrfil = rrfil/max(rrfil);
        m1 = [zeros((N-1),1);u0t1*fh2*k^2/h*rrfil;u0t2*fh1*k^2/h*rrfil];
        m2 = m1;
    case 'sine'
        lhi = fir1(106,2/SR*[10 fcutoff_init]);
        rr = (1-2*rand((N-1),1));
        rrfil = filter(lhi,1,rr);
        rrfil = rrfil/max(rrfil);       
        m1 = [zeros((N-1),1);u0t1*fh2*k^2/h*rrfil;u0t2*fh1*k^2/h*rrfil];
        m2 = m1;
    case 'pluck'
        lhi = fir1(106,2/SR*[10 fcutoff_init]);
        rr = (1-2*rand((N-1),1));
        rrfil = filter(lhi,1,rr);
        rrfil = rrfil/max(rrfil);       
        m1 = [zeros((N-1),1);u0t1*fh2*k^2/h*rrfil;u0t2*fh1*k^2/h*rrfil];
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
        m2(N:2*N-2) = u0v*sin(pi*[1:N-1]/(N+1));
        m2(2*N-1:end) = u0w*(1-2*rand(1,N-1)).*sin(pi*[1:N-1]/(N+1));
        m1 = m2;
end

%% linear components of the system
A00 = sparse([(1+k*sigmal0)*I-k*sigmal1*Dxx, Z, Z; Z, (1+k*sigmatv0)*I-k*sigmatv1*Dxx, Z; Z, Z, (1+k*sigmath0)*I-k*sigmath1*Dxx]);
B00 = sparse([-2*I-(k*gamma*alpha)^2*Dxx, Z, Z; Z, -2*I-(k*gamma)^2*Dxx+(k*kappa)^2*Dxxxx, Z; Z, Z, -2*I-(k*gamma)^2*Dxx+(k*kappa)^2*Dxxxx]);
C00 = sparse([(1-k*sigmal0)*I+k*sigmal1*Dxx, Z, Z; Z, (1-k*sigmatv0)*I+k*sigmatv1*Dxx, Z; Z, Z, (1-k*sigmath0)*I+k*sigmath1*Dxx]);


phisq = (gamma*k)^2*(alpha^2-1)/4; % coupling coefficient


%auxiliary matrices
auxmat1 = phisq*Dxf;



for n = 1:NF
    
    % first derivative vector transverse polarizations
    lambdavec1 = Dxb*m1(N:2*N-2);
    lambdavec2 = Dxb*m1(2*N-1:end);
    Lambdav(vecd) = lambdavec1;
    Lambdah(vecd) = lambdavec2;
    
    % non-linear matrix terms
    Mnlinvl = auxmat1*Lambdav*Dxb;
    Mnlinhl = auxmat1*Lambdah*Dxb;
    Mnlinvh = auxmat1*(Lambdav.^2 + Lambdah.^2)*Dxb;
    
    M00 = sparse([Z, -Mnlinvl, -Mnlinhl; -Mnlinvl, -Mnlinvh, Z; -Mnlinhl, Z, -Mnlinvh]);
    
    % build full system matrices
    A = A00 + M00;
    B = B00 + sparse([Z, Z, Z; -2*Mnlinvl, Z, Z; -2*Mnlinhl, Z, Z]);
    C = C00 + M00;
    

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
        m=A\(-B*m1-C*m2+Jin*fin(n));
        %[m,~,~,~,~] = pcg(A,(-B*m1-C*m2+Jin*fin(n)),1e-15);
    
    elseif syssolv == 1 % Solve using Jacobi method
        count = 0;
        mk1 = m1; % last step state vector
        mk = ones(Ntot,1);
        b = -B*m1-C*m2+Jin*fin(n);
        dd = A(vecdfull); % extract diagonal terms
        D(vecdfull) = dd; % diagonal terms matrix
        Dinv(vecdfull) = 1./dd; % inverse
        R = A-D; % non-diagonal terms matrix
        while sqrt(sum((A*mk-b).^2)/(sum(mk.^2)))>eps
            mk = Dinv*(b-R*mk1); % calculate new state vector
                mk1 = mk; % cycle state vector
            count = count+1;
        end
        itno(n) = count; % iterations counter vector
        m = mk; % 
    end
    
    % write ouput positions
    out1(n) = m(lo);
    out2(n) = m(N-1 +lo);
    out3(n) = m(2*N-2 +lo);
    
    % write output velocities
    outv1(n) = (m(lo)-m1(lo))./k;
    outv2(n) = (m(N-1 +lo)-m1(N-1 +lo))./k;
    outv3(n) = (m(2*N-2 +lo)-m1(2*N-2 +lo))./k;
    
    
    
    
    
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
        
        % string state
        v = m(N:2*N-2);
        v1 = m1(N:2*N-2);
        w  = m(2*N-1:end);
        w1 = m1(2*N-1:end);
        
        % derivatives
        p = Dxb*m(1:N-1); % du/dx, current step
        p1 = Dxb*m1(1:N-1); % du/dx, previous step
        q = Dxb*m(N:N+N-2); % dv/dx, current step
        q1 = Dxb*m1(N:N+N-2); % dv/dx, previous step
        r = Dxb*m(N+N-1:end); % dw/dx, current step
        r1 = Dxb*m1(N+N-1:end); % dw/dx, previous step

        T = 0.5*h*sum((m-m1).^2)/k^2; % kinetic energy
        VLE = 0.5*gamma^2*h*sum(q.*q1) + 0.5*gamma^2*h*sum(r.*r1) + 0.5*gamma^2*alpha^2*h*sum(p.*p1); % potential energy due to elastic terms
        VLS = 0.5*kappa^2*h*(sum((Dxx*v).*(Dxx*v1)) + sum((Dxx*w).*(Dxx*w1)) + ... 
            1/h^4*(2*(1-bc0)*v(1)*v1(1)+2*(1-bcN)*v(end)*v1(end)+2*(1-bc0)*w(1)*w1(1)+2*(1-bcN)*w(end)*w1(end))); % potential energy due to stiffness terms
        VNL = 0.5*gamma^2*(alpha^2-1)*h*(0.5*sum((p+p1).*(q.*q1+r.*r1))+0.25*sum((q.*q1+r.*r1).^2)+0.25*sum((q.*r1-r.*q1).^2));
        H(n) = T+VLE+VLS+VNL;
    end
    
    % cycle stored string states
    m2 = m1;
    m1 = m;
    
    
end

runtime = toc(t1) % get scheme runtime

out_pos = out1+out2+out3; % combined position output
out_vel = outv1+outv2+outv3; % combined velocity output

if plot_on == 1
    
    % spectrum parameters
    nfft = 2^nextpow2(NF);
    df = SR/nfft;
    freqv = 0:df:maxfreq;
    nbins  = numel(freqv);
    win_han = [(sin(pi*[0:NF-1]/(NF-1))).^2]';
    
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
    ylim([-100, 0]);
    xlabel('frequency (Hz)');
    ylabel('amplitude (dB)');
    title('longitudinal polarization spectrum')
    
    subplot(3,1,2)
    vfft2 = 20*log10(magf2/maxf2);
    plot(freqv,vfft2(1:nbins),'k');
    axis tight
    ylim([-100, 0]);
    xlabel('frequency (Hz)');
    ylabel('amplitude (dB)');
    title('vertical polarization spectrum');
    
    subplot(3,1,3)
    vfft3 = 20*log10(magf3/maxf3);
    plot(freqv,vfft3(1:nbins),'k');
    axis tight
    ylim([-100, 0]);
    xlabel('frequency (Hz)');
    ylabel('amplitude (dB)');
    title('horizontal polarization spectrum');
    
    
    fig2 = figure(2); % spectrum plot of the 
    vfft1_all = 20*log10(magf1/fftmax);
    vfft2_all = 20*log10(magf2/fftmax);
    vfft3_all = 20*log10(magf3/fftmax);
    plot(freqv,vfft1_all(1:nbins),'k',freqv,vfft2_all(1:nbins),'b',freqv,vfft3_all(1:nbins),'g');
    axis tight
    %ylim([-100, 0]);
    xlabel('frequency (Hz)');
    ylabel('amplitude (dB)');
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
        plot(k*[1:NF],(H(1)-H)/max(H),'k');
        axis tight
        xlabel('time (seconds)');
        ylabel('energy (adim))');
        title('normalized energy variation');
    end
    
        fig5 = figure(5); %trans displ
    
    subplot(2,1,1)
    plot(k*[1:NF],out2,'k');
    axis tight

    xlabel('time (s)');
    ylabel('vertical displacement (adim)');

    
    subplot(2,1,2)
    plot(k*[1:NF],out3,'k');
    axis tight
    xlabel('time (s)');
    ylabel('horizontal displacement (adim)');
    
end % plotting

% save data
if saveout_on == 1

    dirname = strcat(cd,'\\data\\',SimID,sprintf('_u0v-%1$2.4g_u0w-%2$2.4g_',u0v,u0w),input_type,'_',bcond);
    mkdir(dirname);
    if plot_on == 1
        savefig(fig1,strcat(dirname,'\\spectra.fig'));
        saveas(fig1,strcat(dirname,'\\spectra.eps'),'epsc');
        
        savefig(fig2,strcat(dirname,'\\spectra1pl.fig'));
        saveas(fig2,strcat(dirname,'\\spectra1pl.eps'),'epsc');
        
        savefig(fig3,strcat(dirname,'\\trans_pos.fig'));
        saveas(fig3,strcat(dirname,'\\trans_pos.eps'),'epsc');
        
        savefig(fig5,strcat(dirname,'\\trans_disp.fig'));
        saveas(fig5,strcat(dirname,'\\trans_disp.eps'),'epsc');
        
        if egy_on == 1
            savefig(fig4,strcat(dirname,'\\egyplot.fig'))
            saveas(fig4,strcat(dirname,'\\egyplot.eps'),'epsc')
        end
    end
    out_arr = [out1,out2,out3];
    outv_arr = [outv1,outv2,outv3];
    save(strcat(dirname,'\\wksp.mat'))
    maxtot = max(max(abs(outv_arr)));
    audiowrite(strcat(dirname,'\\audio_ver.wav'),0.8*outv2/maxtot,SR);
    audiowrite(strcat(dirname,'\\audio_hor.wav'),0.8*outv3/maxtot,SR);
    audiowrite(strcat(dirname,'\\audio_lon.wav'),0.8*outv1/maxtot,SR);
    audiowrite(strcat(dirname,'\\audio_norm.wav'),out_vel/max(abs(out_vel)),SR);
    
end