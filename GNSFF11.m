function [y, nmodes] = GNSFF11(k, eps_out, eps_in, rc, r1_vec, r2_vec, nmin, nmax, ...
    i, j, tol, kzimax, direction)
%%% - 09.10.2017
%%% Calculates the total scattered Green's tensor, including all possible
%%%     types of modes
%%% k - k-vector value, 2pi/\lambda_0;
%%% eps_out - outside; eps_in - inside; rc - radius;
%%% r1_vec - reciever position; r2_vec - source position; 1x3 vectors [rho, phi, z]
%%% nmin - minimal mode number in the sum; nmax - maximal mode number in
%%%     the sun;
%%% i, j - tensor indicies over cylindrical coordinates \rho, \phi, z;
%%% tol - relative tolerance for the sum (how many modes 'n' to consider?); |G^{N} - G^{N-1}|/G^{N}; 
%%% kzimax - upper integration limit in kz, usually several k_in;
%%% direction - diraction of propagation - +/- direction;

    r1 = r1_vec(1, 1); r2 = r2_vec(1, 1);
    p1 = r1_vec(1, 2); p2 = r2_vec(1, 2);
    z1 = r1_vec(1, 3); z2 = r2_vec(1, 3);
    %%% CASE OF DIELECTRIC
%     k1 = sqrt(eps_out).*k;   
    k2 = sqrt(eps_in).*k; 
    dkz = -0.000000001;

    kzReMax = kzimax*sqrt(eps_in);        
    kzImMax = 0.0001;   % choose this quantity to be smaller for 
    % larger interatomic distances dz, since exp(1i k_z delta_z) 
    
    % these two tolerances are for the quadgk/integral function
    retol = 1e-9;                             
    abtol = 1e-11;
    
    %%% all modes, both directions
    c_both_start_all = -kzReMax - 1i*dkz;
    c_both_end_all = -c_both_start_all;
    c_both_all = [-1.1*k2 - 1i*dkz, -1.1*k2 + 1i*kzImMax,...
        1.1*k2 - 1i*kzImMax, 1.1*k2 + 1i*dkz];
    
    %%% all modes, forward directions
    c_forw_start_all = 0.0 + 1i*0.0;
    c_forw_end_all = kzReMax + 1i*dkz;
    c_forw_all = [1.1*k2 - 1i*kzImMax, 1.1*k2 + 1i*dkz];
    
    %%% all modes, backward directions
    c_backw_start_all = - c_forw_end_all;
    c_backw_end_all = - c_forw_start_all;
    c_backw_all = [-1.1*k2 - 1i*dkz, -1.1*k2 + 1i*kzImMax];
    
    GNS11mat = 0; % Green's tensor component
    nmodes = 0; % number of considered 'n' modes
    Gnprev = 0; % previous, it sums from nmin to nmax-1; when nmax = 0 Gnprev = 0;
    for num = nmin:1:nmax
% ------------%% REGULAR CASE: ANY STRUCTURE, ALL MODES
               if (direction == 0)
                   GNS11mat = GNS11mat + quadgk(@(x)iGNSFF11(k,eps_out, eps_in, rc, num,x,r1,r2,p1,p2,z1,z2,i,j),...
                       c_both_start_all,c_both_end_all,'RelTol',retol,'AbsTol',abtol,'Waypoints',...
                       c_both_all,'MaxIntervalCount',100000000); 
               elseif (direction == + 1)
                   GNS11mat = GNS11mat + quadgk(@(x)iGNSFF11(k,eps_out, eps_in, rc, num,x,r1,r2,p1,p2,z1,z2,i,j),...
                       c_forw_start_all,c_forw_end_all,'RelTol',retol,'AbsTol',abtol,'Waypoints',...
                       c_forw_all,'MaxIntervalCount',100000000);
               elseif (direction == -1)
                   GNS11mat = GNS11mat + quadgk(@(x)iGNSFF11(k,eps_out, eps_in, rc, num,x,r1,r2,p1,p2,z1,z2,i,j),...
                       c_backw_start_all,c_backw_end_all,'RelTol',retol,'AbsTol',abtol,'Waypoints',...
                       c_backw_all,'MaxIntervalCount',100000000); 
               end
           rel = abs(abs(GNS11mat) - abs(Gnprev))/abs(GNS11mat); % relative contribution of the current mode 'n'
        if (rel < tol) % condition for cutting the 'n'
            break;
        end
        Gnprev = GNS11mat;
        nmodes = nmodes + 1; 
    end
    y = GNS11mat;
end