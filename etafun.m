function  eta = etafun( w, gamma, kappa )
% function for loss coeffficients
eta = 0.5*(-gamma^2+sqrt(gamma^4+4*kappa^2*w^2))/kappa^2;

end

