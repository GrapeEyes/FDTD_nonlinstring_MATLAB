function [Itpmat12,varargout] = gridInterpDirichlet( N1, N2, varargin )
% This function creates the interpolation matrix such that 
% Itpmat12*u1 = u2. Interpolates from N1+1 to N2+1 grid points. We use dirichlet
% boundary conditions, so Itpmat12 is (N2-1)x(N1-1). The inverse
% operator matrix can be also output. Output is sparse.

% linear ('linear') and cubic ('cubic') interpolation are
% implemented. For cubic interpolation, the behaviour near the
% boundaries can be chosen between reverting to linear
% interpolation ('revlin'), clamped boundary conditions
% ('clamped') or simply supported boundary conditions
% ('ssup'), as well as mixed boundary conditions 
% ('clamped-ssup'),('ssup-clamped')

% the order of the arguments is N1, N2, interpolation type, boundary conditions 

% function arguments check 
if nargin < 2
    error('not enough arguments');
elseif nargin == 2
    itpor = 'linear'; % interpolation polynomial order
elseif nargin == 3
    itpor = varargin{1};
    bcond = 'revlin';
elseif nargin == 4
    itpor = varargin{1};
    bcond = varargin{2};
else
    error('too many arguments');
end

flipflag = 0; % see below

N = [N1,N2]; % see below

if N2>N1 % this changes the order to make N(1) to N(2) always the downsampling operation
    N = flip(N);
    flipflag = 1;
elseif N2==N1
    Itpmat12 = speye(N1-1,N1-1); % downsampling matrix
    Itpmat21 = Itpmat12; % upsampling matrix
    varargout{1} = Itpmat21; % interpolation matrix N2 to N1
    return
    
end

h = 1./N; % grid spacings vector
Imat1 = sparse(N(1)-1,N(2)-1); % initialize upsampling matrix
Imat2 = Imat1'; % initialize downsampling matrix

hquot = h(1)/h(2); % ratio of the grid spacings
ind1 = [1:N(1)-1]; % indexes grid 1
l1 = floor(hquot*ind1); % truncated index in grid 2
alpha1 = hquot*ind1-l1; % truncation remain of l1
ind2 = [1:N(2)-1]; % indexes grid 2
l2 = floor(ind2/hquot); % truncated index in grid 1
alpha2 = ind2/hquot-l2; % truncation remain of l2



switch itpor % order of interpolation to be used

    % Linear interpolation
    case 'linear'  % linear interpolation
        
        % upsampling matrix
        for n = ind1
            if l1(n) == 0
                Imat1(n,1) = alpha1(n);
            elseif l1(n) == N(2)-1
                Imat1(n,N(2)-1) = 1-alpha1(n);
            else
                Imat1(n,l1(n)+1) = alpha1(n);
                Imat1(n,l1(n)) = 1-alpha1(n);
            end
        end
        
        % downsampling matrix
        ind2 = [1:N(2)-1];
        l2 = floor(ind2/hquot);
        alpha2 = ind2/hquot-l2;
        for n = ind2
            Imat2(n,l2(n)+1) = alpha2(n);
            Imat2(n,l2(n)) = 1-alpha2(n);
        end
        
    % Cubic interpolation
    case 'cubic' % cubic interpolation
        switch bcond % boundary condition determines interpolation at the ends of the string for the upsampling matrix
            
            case 'revlin' % reverts to linear interpolation for points between l=0,l=1 and l=N-1,l=N
                
                % upsampling matrix
                for n = ind1
                    
                    if l1(n) == 0
                        Imat1(n,1) = alpha1(n);
                    
                    elseif l1(n) == N(2)-1
                        Imat1(n,N(2)-1) = 1-alpha1(n);
                    
                    elseif l1(n) == 1
                        Imat1(n,1) = (alpha1(n)-1)*(alpha1(n)+1)*(alpha1(n)-2)/2;
                        Imat1(n,2) = -alpha1(n)*(alpha1(n)+1)*(alpha1(n)-2)/2;
                        Imat1(n,3) = alpha1(n)*(alpha1(n)+1)*(alpha1(n)-1)/6;
                    
                    elseif l1(n) == N(2)-2
                        Imat1(n,N(2)-3) = -alpha1(n)*(alpha1(n)-1)*(alpha1(n)-2)/6;
                        Imat1(n,N(2)-2) = (alpha1(n)-1)*(alpha1(n)+1)*(alpha1(n)-2)/2;
                        Imat1(n,N(2)-1) = -alpha1(n)*(alpha1(n)+1)*(alpha1(n)-2)/2;
                        
                    else
                        Imat1(n,l1(n)-1) = -alpha1(n)*(alpha1(n)-1)*(alpha1(n)-2)/6;
                        Imat1(n,l1(n)) = (alpha1(n)-1)*(alpha1(n)+1)*(alpha1(n)-2)/2;
                        Imat1(n,l1(n)+1) = -alpha1(n)*(alpha1(n)+1)*(alpha1(n)-2)/2;
                        Imat1(n,l1(n)+2) = alpha1(n)*(alpha1(n)+1)*(alpha1(n)-1)/6;
                    end
                end %for loop upsampling matrix with linear interpolation at the ends
                
            
            case 'clamped' % use clamped boundary conditions
                % upsampling matrix
                for n = ind1
                    
                    if l1(n) == 0
                        Imat1(n,1) = -alpha1(n)*(alpha1(n)+1)*(alpha1(n)-2)/2 - alpha1(n)*(alpha1(n)-1)*(alpha1(n)-2)/6; 
                        Imat1(n,2) = alpha1(n)*(alpha1(n)+1)*(alpha1(n)-1)/6;
                    
                    elseif l1(n) == N(2)-1
                        Imat1(n,N(2)-2) = -alpha1(n)*(alpha1(n)-1)*(alpha1(n)-2)/6;
                        Imat1(n,N(2)-1) = (alpha1(n)-1)*(alpha1(n)+1)*(alpha1(n)-2)/2 + alpha1(n)*(alpha1(n)+1)*(alpha1(n)-1)/6;; 
                        
                    elseif l1(n) == 1
                        Imat1(n,1) = (alpha1(n)-1)*(alpha1(n)+1)*(alpha1(n)-2)/2;
                        Imat1(n,2) = -alpha1(n)*(alpha1(n)+1)*(alpha1(n)-2)/2;
                        Imat1(n,3) = alpha1(n)*(alpha1(n)+1)*(alpha1(n)-1)/6;
                    
                    elseif l1(n) == N(2)-2
                        Imat1(n,N(2)-3) = -alpha1(n)*(alpha1(n)-1)*(alpha1(n)-2)/6;
                        Imat1(n,N(2)-2) = (alpha1(n)-1)*(alpha1(n)+1)*(alpha1(n)-2)/2;
                        Imat1(n,N(2)-1) = -alpha1(n)*(alpha1(n)+1)*(alpha1(n)-2)/2;
                        
                    else
                        Imat1(n,l1(n)-1) = -alpha1(n)*(alpha1(n)-1)*(alpha1(n)-2)/6;
                        Imat1(n,l1(n)) = (alpha1(n)-1)*(alpha1(n)+1)*(alpha1(n)-2)/2;
                        Imat1(n,l1(n)+1) = -alpha1(n)*(alpha1(n)+1)*(alpha1(n)-2)/2;
                        Imat1(n,l1(n)+2) = alpha1(n)*(alpha1(n)+1)*(alpha1(n)-1)/6;
                    end
                end %for loop upsampling matrix with clampled boundary conditions
                
            case 'ssup' % use simply supported boundary conditions
                
                % upsampling matrix
                for n = ind1
                    
                    if l1(n) == 0
                        Imat1(n,1) = -alpha1(n)*(alpha1(n)+1)*(alpha1(n)-2)/2 + alpha1(n)*(alpha1(n)-1)*(alpha1(n)-2)/6;
                        Imat1(n,2) = alpha1(n)*(alpha1(n)+1)*(alpha1(n)-1)/6;
                    
                    elseif l1(n) == N(2)-1
                        Imat1(n,N(2)-2) = -alpha1(n)*(alpha1(n)-1)*(alpha1(n)-2)/6;
                        Imat1(n,N(2)-1) = (alpha1(n)-1)*(alpha1(n)+1)*(alpha1(n)-2)/2 - alpha1(n)*(alpha1(n)+1)*(alpha1(n)-1)/6;
                        
                    elseif l1(n) == 1
                        Imat1(n,1) = (alpha1(n)-1)*(alpha1(n)+1)*(alpha1(n)-2)/2;
                        Imat1(n,2) = -alpha1(n)*(alpha1(n)+1)*(alpha1(n)-2)/2;
                        Imat1(n,3) = alpha1(n)*(alpha1(n)+1)*(alpha1(n)-1)/6;
                    
                    elseif l1(n) == N(2)-2
                        Imat1(n,N(2)-3) = -alpha1(n)*(alpha1(n)-1)*(alpha1(n)-2)/6;
                        Imat1(n,N(2)-2) = (alpha1(n)-1)*(alpha1(n)+1)*(alpha1(n)-2)/2;
                        Imat1(n,N(2)-1) = -alpha1(n)*(alpha1(n)+1)*(alpha1(n)-2)/2;
                        
                    else
                        Imat1(n,l1(n)-1) = -alpha1(n)*(alpha1(n)-1)*(alpha1(n)-2)/6;
                        Imat1(n,l1(n)) = (alpha1(n)-1)*(alpha1(n)+1)*(alpha1(n)-2)/2;
                        Imat1(n,l1(n)+1) = -alpha1(n)*(alpha1(n)+1)*(alpha1(n)-2)/2;
                        Imat1(n,l1(n)+2) = alpha1(n)*(alpha1(n)+1)*(alpha1(n)-1)/6;
                    end
                end %for loop upsampling matrix with simply supported boundary conditions
                
                case 'clamped-ssup' % use clamped boundary condition on 0, ssup on L
                % upsampling matrix
                for n = ind1
                    
                    if l1(n) == 0
                        Imat1(n,1) = -alpha1(n)*(alpha1(n)+1)*(alpha1(n)-2)/2 - alpha1(n)*(alpha1(n)-1)*(alpha1(n)-2)/6; 
                        Imat1(n,2) = alpha1(n)*(alpha1(n)+1)*(alpha1(n)-1)/6;
                    
                    elseif l1(n) == N(2)-1
                        Imat1(n,N(2)-2) = -alpha1(n)*(alpha1(n)-1)*(alpha1(n)-2)/6;
                        Imat1(n,N(2)-1) = (alpha1(n)-1)*(alpha1(n)+1)*(alpha1(n)-2)/2 - alpha1(n)*(alpha1(n)+1)*(alpha1(n)-1)/6;; 
                        
                    elseif l1(n) == 1
                        Imat1(n,1) = (alpha1(n)-1)*(alpha1(n)+1)*(alpha1(n)-2)/2;
                        Imat1(n,2) = -alpha1(n)*(alpha1(n)+1)*(alpha1(n)-2)/2;
                        Imat1(n,3) = alpha1(n)*(alpha1(n)+1)*(alpha1(n)-1)/6;
                    
                    elseif l1(n) == N(2)-2
                        Imat1(n,N(2)-3) = -alpha1(n)*(alpha1(n)-1)*(alpha1(n)-2)/6;
                        Imat1(n,N(2)-2) = (alpha1(n)-1)*(alpha1(n)+1)*(alpha1(n)-2)/2;
                        Imat1(n,N(2)-1) = -alpha1(n)*(alpha1(n)+1)*(alpha1(n)-2)/2;
                        
                    else
                        Imat1(n,l1(n)-1) = -alpha1(n)*(alpha1(n)-1)*(alpha1(n)-2)/6;
                        Imat1(n,l1(n)) = (alpha1(n)-1)*(alpha1(n)+1)*(alpha1(n)-2)/2;
                        Imat1(n,l1(n)+1) = -alpha1(n)*(alpha1(n)+1)*(alpha1(n)-2)/2;
                        Imat1(n,l1(n)+2) = alpha1(n)*(alpha1(n)+1)*(alpha1(n)-1)/6;
                    end
                end %for loop upsampling matrix with clampled boundary conditions
                
                case 'ssup-clamped' % use clamped boundary conditions
                % upsampling matrix
                for n = ind1
                    
                    if l1(n) == 0
                        Imat1(n,1) = -alpha1(n)*(alpha1(n)+1)*(alpha1(n)-2)/2 + alpha1(n)*(alpha1(n)-1)*(alpha1(n)-2)/6; 
                        Imat1(n,2) = alpha1(n)*(alpha1(n)+1)*(alpha1(n)-1)/6;
                    
                    elseif l1(n) == N(2)-1
                        Imat1(n,N(2)-2) = -alpha1(n)*(alpha1(n)-1)*(alpha1(n)-2)/6;
                        Imat1(n,N(2)-1) = (alpha1(n)-1)*(alpha1(n)+1)*(alpha1(n)-2)/2 + alpha1(n)*(alpha1(n)+1)*(alpha1(n)-1)/6;; 
                        
                    elseif l1(n) == 1
                        Imat1(n,1) = (alpha1(n)-1)*(alpha1(n)+1)*(alpha1(n)-2)/2;
                        Imat1(n,2) = -alpha1(n)*(alpha1(n)+1)*(alpha1(n)-2)/2;
                        Imat1(n,3) = alpha1(n)*(alpha1(n)+1)*(alpha1(n)-1)/6;
                    
                    elseif l1(n) == N(2)-2
                        Imat1(n,N(2)-3) = -alpha1(n)*(alpha1(n)-1)*(alpha1(n)-2)/6;
                        Imat1(n,N(2)-2) = (alpha1(n)-1)*(alpha1(n)+1)*(alpha1(n)-2)/2;
                        Imat1(n,N(2)-1) = -alpha1(n)*(alpha1(n)+1)*(alpha1(n)-2)/2;
                        
                    else
                        Imat1(n,l1(n)-1) = -alpha1(n)*(alpha1(n)-1)*(alpha1(n)-2)/6;
                        Imat1(n,l1(n)) = (alpha1(n)-1)*(alpha1(n)+1)*(alpha1(n)-2)/2;
                        Imat1(n,l1(n)+1) = -alpha1(n)*(alpha1(n)+1)*(alpha1(n)-2)/2;
                        Imat1(n,l1(n)+2) = alpha1(n)*(alpha1(n)+1)*(alpha1(n)-1)/6;
                    end
                end %for loop upsampling matrix with clampled boundary conditions
                
        end % cubic interpolation upsampling matrix (depends on boundary conditions)
        
        % cubic interpolation downsampling matrix (doesn't depend on boundary conditions)        
        for n = ind2
            if l2(n) == 1
                Imat2(n,1) = (alpha2(n)-1)*(alpha2(n)+1)*(alpha2(n)-2)/2;
                Imat2(n,2) = -alpha2(n)*(alpha2(n)+1)*(alpha2(n)-2)/2;
                Imat2(n,3) = alpha2(n)*(alpha2(n)+1)*(alpha2(n)-1)/6;
                
            elseif l2(n) == N(1)-2
                Imat2(n,N(1)-3) = -alpha2(n)*(alpha2(n)-1)*(alpha2(n)-2)/6;
                Imat2(n,N(1)-2) = (alpha2(n)-1)*(alpha2(n)+1)*(alpha2(n)-2)/2;
                Imat2(n,N(1)-1) = -alpha2(n)*(alpha2(n)+1)*(alpha2(n)-2)/2;
                
            else
                Imat2(n,l2(n)-1) = -alpha2(n)*(alpha2(n)-1)*(alpha2(n)-2)/6;
                Imat2(n,l2(n)) = (alpha2(n)-1)*(alpha2(n)+1)*(alpha2(n)-2)/2;
                Imat2(n,l2(n)+1) = -alpha2(n)*(alpha2(n)+1)*(alpha2(n)-2)/2;
                Imat2(n,l2(n)+2) = alpha2(n)*(alpha2(n)+1)*(alpha2(n)-1)/6;
            end
        end %for loop that creates cubic downsampling matrix
        
        
    % Otherwise      
    otherwise % if the interp order given is not valid, use linear interpolation 
        warning('%s interpolation not available, defaulting to linear interpolation',itpor);
        
        % upsampling matrix
        for n = ind1
            if l1(n) == 0
                Imat1(n,1) = alpha1(n);
            elseif l1(n) == N(2)-1
                Imat1(n,N(2)-1) = 1-alpha1(n);
            else
                Imat1(n,l1(n)+1) = alpha1(n);
                Imat1(n,l1(n)) = 1-alpha1(n);
            end
        end
        
        % downsampling matrix
        for n = ind2
            Imat2(n,l2(n)+1) = alpha2(n);
            Imat2(n,l2(n)) = 1-alpha2(n);
        end
      
end


% assign depending on the order of inputs
if flipflag == 0 % N1>N2
    Itpmat12 = Imat2; % downsampling matrix
    Itpmat21 = Imat1; % upsampling matrix
else % N2>N1
    Itpmat12 = Imat1; % upsampling matrix
    Itpmat21 = Imat2; % downsampling matrix
end


varargout{1} = Itpmat21; % interpolation matrix N2 to N1

end %function