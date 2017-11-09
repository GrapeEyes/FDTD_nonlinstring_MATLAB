% nlinstring_basic_nonorm

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Theoretical description
% This script models nonlinear string oscillation, as a
% system of 2 coupled pdes, for the longitudinal and
% transversal oscillations respectively. Doesn't include
% stiffnes or loss, and neither the coordinate x or the
% displacements are normalized. The grids are equal for both
% polarizations.

% Eq1.: rho*Ar*vtt = E*Ar*vxx-(E*Ar-T0)*(dphi/dp)x
% Eq2.: rho*Ar*utt = E*Ar*uxx-(E*Ar-T0)*(dphi/dq)x
% v : longitudinal displacement
% u : transversal displacement
% p = vx; q = ux;
% phi = sqrt((1-p)^2+q^2)-1-p,
% approx as = 1/2*q^2 - 1/2*p*q^2 - 1/8*q^4
% dphi/dp = -1/2*q^2
% dphi/dq = -1/2*(q^3 + 2*p*q + [-2q]) [] term is omitted 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Discretized scheme
% w = [v;u]
% A*w[n+1]+B*w[n]+C*w[n-1] = 0, with A, B, C 2x2 block matrices

% A11=C11 = I; A12=C12 = -beta*Dxforw*Lambdat*Dxbackw;
% B11 = -2I-alpha*Dxx
% A22=C22 = I - beta*Dxforw*(Lamdat)^2*Dxbackw - beta*Dxforw*Lambdal*Dxback
% B22  = -2*I - alpha*Dxx - 2*beta*Dxforw*Lambdal*Dxbackw
% with
% alpha = k^2*E/rho, beta = 0.25*k^2*(E*Ar-T0)/(rho*Ar)
% I = identity, Lambdat = diag(Dxback*u[n]), Lambdal = diag(Dxback*v[n])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Boundary conditions, stability & initialization
% Dirichlet boundary conditions at x=0, x=L
% We use the linear limit for the stability condition,
% defining the courant number as lambda = ct*k/h, and
% lambda<=1, hmin = ct*k, with ct = T0/(rho*Ar) the
% transversal wave velocity

% Three types of input can be used, either initialization 
% with a raised cosine o triangular distribution, or a short raised cosine
% in time impulse for a given grid point   


















