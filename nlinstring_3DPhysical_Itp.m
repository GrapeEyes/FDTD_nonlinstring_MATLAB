% nlinstring_3DPhysical_Itp

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Carlos de la Vega Martin, 22/08/2017
% MSc Acoustics and Music Technology at the University of Edinburgh.

% Notes: developed on MATLAB R2016b for the masters thesis.
% needs two external functions, gridInterpDirichlet.m,
% etafun.m



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Boundary conditions, stability & initialization
% Dirichlet boundary conditions at x=0, x=L. Clamped or
% simply supported boundary conditions can be chosen.
% We use the linear limit for the stability condition.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all

t1 = tic; % start timer

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% flags

syssolv = 0; % method used for solving the linear system backlash = 0, Jacobi = 1
egy_on = 1; % calculate energy
ddtest_on = 0; % test that the system matrix is diagonally dominant
plot_on = 1; % plot spectra, and if calculated, energy
saveout_on = 1; % save output as .wav, variables as .mat, and if produced, plots as .fig and .eps
str_param_fmt = 1; % physical=0, fundamentals+inharm=1  eq coeff=2
loss_fmt = 1; % sigmas=0, T60+freq=1
norm_out_on = 1; % set to 1 normalizes the out vel before saving it to .wav, 2 normalizes each comp before adding

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% parameters

% simpulation parameters
SR = 48000; % sampling rate
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
        B = 1.51e-4; % inharmonicity
        
    case 2 % scheme parameters directly
        
        gamma = 200; % 2*fundamental of transverse mode
        alpha = 15; % 2*fundamental/gamma of longitudinal mode
        kappa = 1.5; % stiffness coefficient
end

% boundary conditions

%'clamped' , 'ssup' , 'clamped-ssup' , 'ssup-clamped'
bcond = 'ssup'; % for the Dxxxx operator
bconditp = bcond; % for the interpolating matrix

% loss parameters

switch loss_fmt
    case 0 % loss coefficients directly
        
        sigmal0 = 0.0; % 1st order loss coeff, longitudinal pol
        sigmal1 = 0.0; % 3rd order loss coeff, longitudinal pol
        sigmat0 = 0.0; % 1st order loss coeff, transverse pol
        sigmat1 = 0.0; % 3rd order loss coeff, transverse pol
        
    case 1 % loss fitted to T60 at two freq f1<f2, T601>T602
        fq1 = 100; % lower frequency for loss profile, Hz 
        fq2 = 2000; % higher frequency for loss profile, Hz
        
        T60l_fq1 = 6; % T60 longitudinal pol at freq 1, s
        T60l_fq2 = 4; % T60 longitudinal pol at freq 2, s
        T60t_fq1 = 6; % T60 transverse pol at freq 1, s
        T60t_fq2 = 4; % T60 tranverse pol at freq 2, s
        
end


% I/O parameters

input_type = 'randin'; % strike, sine, pluck
% for 0,1 an initial bandpassed noise can be added to any of the polarizations

switch input_type
    
    case 'strike' % raised cosine strike
        
        xh = 0.7; % excitation position
        Th = 0.003; % excitation time (s)
        fh = 1000.0; % excitation strength (N/(Kg/m))
        angle = pi/2; % incidence angle, respect to the vertical (rad)
        
        u0t1 = 0.0; % maximum initial displacement, vert transv pol, rel to max disp in 2
        u0t2 = 0.0; % maximum initial displacement, horiz transv pol, rel to max disp in 1
        fcutoff_init = 500; % cutoff frequency for bandpassing the noise term
        
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
maxfreq = 3000; % upper frequency limit for the spectrum plots
trans_time = 0.0; % for transverse movement plot, defines center of the transition phase, normalized respect to TF

% save parameters
SimID = 'E3exampleitp'; % name of the folder, use an identifiable tag, like 'instrument_note_parameters'


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
    kappa = gamma/pi*sqrt(B); % stiffness coefficient
    
end
    

% loss

if loss_fmt == 1
    eta1 = etafun(2*pi*fq1,gamma,kappa);
    eta2 = etafun(2*pi*fq2,gamma,kappa);
    
    sigmal0 = 6*log(10)/(eta2-eta1)*(eta2/T60l_fq1 - eta1/T60l_fq2); % 1st order loss coeff, longitudinal pol
    sigmal1 = 6*log(10)/(eta2-eta1)*(-1/T60l_fq1 + 1/T60l_fq2); % 3rd order loss coeff, longitudinal pol
    sigmat0 = 6*log(10)/(eta2-eta1)*(eta2/T60t_fq1 - eta1/T60t_fq2); % 1st order loss coeff, transverse pol
    sigmat1 = 6*log(10)/(eta2-eta1)*(-1/T60t_fq1 + 1/T60t_fq2);  % 3rd order loss coeff, transverse pol

end


    % grid

    ht = sqrt(0.5*(gamma^2*k^2+sqrt(gamma^4*k^4+16*kappa^2*k^2))); % minimum transversal grid spacing
    hl = gamma*alpha*k; % minimum longitudinal grid spacing
    
    Nl = floor(1/hl); % number of segments, longitudinal grid
    Nt = floor(1/ht); % number of segments, transversal grid
    hl = 1/Nl; % longitudinal grid spacing
    ht = 1/Nt; % transversal grid spacing
    N = Nl+2*Nt-3; % size of the system matrix


    % output position
    lot = floor(xo*(Nt-2))+2; % grid index of output (tranversal oscillation)
    lol = floor(xo*(Nl-2))+2; % grid index of output (longitudinal oscillation)


    %% finite difference matrix operators

    % tranversal grid operators
    et = ones(Nt,1);
    Dxtforw = spdiags([-et et], 0:1,Nt-1,Nt)/ht; % d/dx, forward difference
    Dxtbackw = -Dxtforw'; % d/dx, backward difference
    Dxxt = Dxtforw*Dxtbackw; % d2/dx2, centered difference
    Dxxxxt = Dxxt*Dxxt; % d4/dx4, centered difference

    % change values for u1, uN-1 depending on the boundary conditions
    switch bcond
        case 'clamped'
            Dxxxxt(1) = 7/ht^4;
            Dxxxxt(Nt^2+1-2*Nt) = 7/ht^4;
        case 'clamped-ssup'
            Dxxxxt(1) = 7/ht^4;
        case 'ssup-clamped'
            Dxxxxt(Nt^2+1-2*Nt) = 7/ht^4;
    end


    % longitudinal grid
    el = ones(Nl,1);
    Dxlforw = spdiags([-el el], 0:1,Nl-1,Nl)/hl; %d/dx, forward difference
    Dxlbackw = -Dxlforw'; % d/dx, backward difference
    Dxxl = Dxlforw*Dxlbackw; % d2/dx2, centered difference

    % interpolation
    Inplt = gridInterpDirichlet(Nl,Nt,'cubic', bconditp); % interpolation matrices
    Inptl = ht/hl*Inplt';

    % create identity and zero matrices of required dimensions
    It = speye(Nt-1); % identity matrix transversal grid
    Il = speye(Nl-1); % identity matrix longitudinal grid
    Z = sparse(Nl-1,Nt-1); % zero matrix crossterms 1
    Ztr = Z'; % zero matrix crossterms 2
    Zt = sparse(Nt-1,Nt-1); % zero matrix transversal grid
    Zl = sparse(Nl-1,Nl-1); % zero matrix longitudinal grid


    % initialize diagonal matrices for the first derivative of
    % the transverse polarization
    Lambdav = sparse(Nt,Nt); Lambdah = Lambdav; 

    % initialize system matrices off-diagonal blocks
    Mnlin12 = Z; Mnlin21 = Z'; Mnlin13 = Mnlin12; Mnlin31 = Mnlin21;
    Mnlin22 = sparse(Nt-1,Nt-1); Mnlin33 = Mnlin22;

    % initialize outputs
    out1 = zeros(NF,1); out2 = out1; out3 = out1;
    H = out1;  itno = out1;
    outv1 = out1; outv2 = out1; outv3 = out1;

    % initialize system matrices
    A = sparse(N,N);
    B = A; C = A; D=A; Dinv = A; R=A; dd = zeros(N,1);

    % vectors for linear indexing of the diagonal components
    vecd = [0:Nt-1]*Nt+[1:Nt];
    vecdfull = [0:(N-1)]*N+[1:N];



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


Jin = zeros(N,1); % initialize vector for force to displacement conversion
xhint = floor(xh*Nt); % input position in samples
Jin(Nl-1+xhint) = fh1*k^2/ht; % factor to convert from force to displacement, applied to only one grid point, vertical pol
Jin(Nl+Nt-2+xhint) = fh2*k^2/ht; % factor to convert from force to displacement, applied to only one grid point, horiz pol

% initialize string
m = zeros(N,1);
switch input_type
    case 'strike'
        lhi = fir1(106,2/SR*[10 fcutoff_init]);
        rr = (1-2*rand((Nt-1),1));
        rrfil = filter(lhi,1,rr);
        rrfil = rrfil/max(rrfil);
        m1 = [zeros((Nl-1),1);u0t1*fh2*k^2/ht*rrfil;u0t2*fh1*k^2/ht*rrfil];
        m2 = m1;
    case 'sine'
        lhi = fir1(106,2/SR*[10 fcutoff_init]);
        rr = (1-2*rand((Nt-1),1));
        rrfil = filter(lhi,1,rr);
        rrfil = rrfil/max(rrfil);       
        m1 = [zeros((Nl-1),1);u0t1*fh2*k^2/ht*rrfil;u0t2*fh1*k^2/ht*rrfil];
        m2 = m1;
    case 'pluck'
        lhi = fir1(106,2/SR*[10 fcutoff_init]);
        rr = (1-2*rand((Nt-1),1));
        rrfil = filter(lhi,1,rr);
        rrfil = rrfil/max(rrfil);       
        m1 = [zeros((Nl-1),1);u0t1*fh2*k^2/ht*rrfil;u0t2*fh1*k^2/ht*rrfil];
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
        m2(Nl:Nl+Nt-2) = u0v*sin(pi*[1:Nt-1]/(Nt+1));
        
        m2(Nl+Nt-1:end) = u0w*(1-2*rand(1,Nt-1)).*sin(pi*[1:Nt-1]/(Nt+1));
        %(1-2*rand(1,Nt-1))
        m1 = m2;
end

%% linear components of the system
A00 = sparse([(1+k*sigmal0)*Il-k*sigmal1*Dxxl, Z, Z; Ztr, (1+k*sigmat0)*It-k*sigmat1*Dxxt, Zt; Ztr, Zt, (1+k*sigmat0)*It-k*sigmat1*Dxxt]);
B00 = sparse([-2*Il-(k*gamma*alpha)^2*Dxxl, Z, Z; Ztr, -2*It-(k*gamma)^2*Dxxt+(k*kappa)^2*Dxxxxt, Zt; Ztr, Zt, -2*It-(k*gamma)^2*Dxxt+(k*kappa)^2*Dxxxxt]);
C00 = sparse([(1-k*sigmal0)*Il+k*sigmal1*Dxxl, Z, Z; Ztr, (1-k*sigmat0)*It+k*sigmat1*Dxxt, Zt; Ztr, Zt, (1-k*sigmat0)*It+k*sigmat1*Dxxt]);


phisq = (gamma*k)^2*(alpha^2-1)/4; % coupling coefficient


%auxiliary matrices
auxmat1 = phisq*Inptl*Dxtforw;
auxmat2 = phisq*Dxtforw;
auxmat3 = Dxtbackw*Inplt;


for n = 1:NF
    
    % first derivative vector transverse polarizations
    lambdavec1 = Dxtbackw*m1(Nl:Nl+Nt-2);
    lambdavec2 = Dxtbackw*m1(Nl+Nt-1:end);
    Lambdav(vecd) = lambdavec1;
    Lambdah(vecd) = lambdavec2;
    
    % non-linear matrix terms
    Mnlin12 = auxmat1*Lambdav*Dxtbackw;
    Mnlin21 = auxmat2*Lambdav*auxmat3;
    Mnlin13 = auxmat1*Lambdah*Dxtbackw;
    Mnlin31 = auxmat2*Lambdah*auxmat3;
    Mnlin22 = auxmat2*(Lambdav.^2 + Lambdah.^2)*Dxtbackw;
    Mnlin33 = auxmat2*(Lambdav.^2 + Lambdah.^2)*Dxtbackw;
    M00 = sparse([Zl, -Mnlin12, -Mnlin13; -Mnlin21, -Mnlin22, Zt; -Mnlin31, Zt, -Mnlin33]);
    
    % build full system matrices
    A = A00 + M00;
    B = B00 + sparse([Zl, Z, Z; -2*Mnlin21, Zt, Zt; -2*Mnlin31, Zt, Zt]);
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
        m = A\(-B*m1-C*m2+Jin*fin(n));
    
    elseif syssolv == 1 % Solve using Jacobi method
        count = 0;
        mk1 = m1; % last step state vector
        mk = ones(Nl+2*Nt-3,1);
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
    out1(n) = m(lol);
    out2(n) = m(Nl-1 +lot);
    out3(n) = m(Nl+Nt-2 +lot);
    
    % write output velocities
    outv1(n) = (m(lol)-m1(lol))./k;
    outv2(n) = (m(Nl-1 +lot)-m1(Nl-1 +lot))./k;
    outv3(n) = (m(Nl+Nt-2 +lot)-m1(Nl+Nt-2 +lot))./k;
    
    
    
    
    
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
        v = m(Nl:Nl+Nt-2);
        v1 = m1(Nl:Nl+Nt-2);
        w  = m(Nl+Nt-1:end);
        w1 = m1(Nl+Nt-1:end);
        
        % derivatives
        p = Dxtbackw*Inplt*m(1:Nl-1); % du/dx, current step
        p1 = Dxtbackw*Inplt*m1(1:Nl-1); % du/dx, previous step
        q = Dxtbackw*m(Nl:Nl+Nt-2); % dv/dx, current step
        q1 = Dxtbackw*m1(Nl:Nl+Nt-2); % dv/dx, previous step
        r = Dxtbackw*m(Nl+Nt-1:end); % dw/dx, current step
        r1 = Dxtbackw*m1(Nl+Nt-1:end); % dw/dx, previous step

        T = 0.5*([hl*el(1:end-1)',ht*et(1:end-1)',ht*et(1:end-1)']*(m-m1).^2)/k^2; % kinetic energy
        VLE = 0.5*gamma^2*ht*sum(q.*q1) + 0.5*gamma^2*ht*sum(r.*r1) + 0.5*gamma^2*alpha^2*hl*sum((Dxlbackw*m(1:Nl-1)).*(Dxlbackw*m1(1:Nl-1))); % potential energy due to elastic terms
        VLS = 0.5*kappa^2*ht*(sum((Dxxt*v).*(Dxxt*v1)) + sum((Dxxt*w).*(Dxxt*w1)) + ... 
            1/ht^4*(2*(1-bc0)*v(1)*v1(1)+2*(1-bcN)*v(end)*v1(end)+2*(1-bc0)*w(1)*w1(1)+2*(1-bcN)*w(end)*w1(end))); % potential energy due to stiffness terms
        VNL = 0.5*gamma^2*(alpha^2-1)*ht*(0.5*sum((p+p1).*(q.*q1+r.*r1))+0.25*sum((q.*q1+r.*r1).^2)+0.25*sum((q.*r1-r.*q1).^2));
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
    ylim([-100,0])
    xlabel('frequency (Hz)');
    ylabel('amplitude (dB)');
    title('longitudinal polarization spectrum')
    
    subplot(3,1,2)
    vfft2 = 20*log10(magf2/maxf2);
    plot(freqv,vfft2(1:nbins),'k');
    axis tight
    ylim([-100,0])
    xlabel('frequency (Hz)');
    ylabel('amplitude (dB)');
    title('vertical polarization spectrum');
    
    subplot(3,1,3)
    vfft3 = 20*log10(magf3/maxf3);
    plot(freqv,vfft3(1:nbins),'k');
    axis tight
    ylim([-100,0])
    xlabel('frequency (Hz)');
    ylabel('amplitude (dB)');
    title('horizontal polarization spectrum');
    
    
    fig2 = figure(2); % spectrum plot of the 
    vfft1_all = 20*log10(magf1/fftmax);
    vfft2_all = 20*log10(magf2/fftmax);
    vfft3_all = 20*log10(magf3/fftmax);
    plot(freqv,vfft1_all(1:nbins),'k',freqv,vfft2_all(1:nbins),'b',freqv,vfft3_all(1:nbins),'g');
    axis tight
    ylim([-100,0])
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
        ylabel('energy (adim)');
        title('normalized energy variation');
    end
    
    fig5 = figure(5); % spectrum plots for the three polarizations, separated
    
    subplot(2,1,1)
    plot(k*[1:NF],out3,'k');
    axis tight

    xlabel('time (s)');
    ylabel('vertical displacement (adim)');

    
    subplot(2,1,2)
    plot(k*[1:NF],out3,'k');
    axis tight
    xlabel('time (s)');
    ylabel('horizontal displacement (adim)');
    
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

    dirname = strcat(cd,'\\data\\',SimID,sprintf('\\u0v-%1$2.4g_u0w-%1$2.4g_',u0v,u0w),bcond);
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
    maxtot = max(max(outv_arr));
    save(strcat(dirname,'\\wksp.mat'))
    audiowrite(strcat(dirname,'\\audio_ver.wav'),0.8*outv2/maxtot,SR);
    audiowrite(strcat(dirname,'\\audio_hor.wav'),0.8*outv3/maxtot,SR);
    audiowrite(strcat(dirname,'\\audio_lon.wav'),0.8*outv1/maxtot,SR);
    audiowrite(strcat(dirname,'\\audio_norm.wav'),0.8*out_vel/max(abs(out_vel)),SR);
    
end