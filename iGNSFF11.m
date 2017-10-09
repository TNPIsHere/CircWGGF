function y = iGNSFF11(k, eps_out, eps_in, rc, n, x, rr, rs, pr, ps, zr, zs, i, j)
%%% - 09.10.2017 - 
%%% calculates the integrand of the scattered part of the Green's tensor of the fiber;
%%% k - k-vector value, 2pi/\lambda_0;
%%% eps_out, eps_in - eps outside and inside
%%% rc - radius
%%% n - mode order
%%% x - k_z, integration variable
%%% rr, pr, zr - reciever coordinates
%%% rs, ps, zs - source coordinates
%%% i, j - rho, phi, z tensor indeces

    k1 = sqrt(eps_out).*k;     k2 = sqrt(eps_in).*k;
 	a = sqrt(eps_out.*k.^2 - x.^2); b = sqrt(eps_in.*k.^2 - x.^2);

    Hn1r = besselh(n,a.*rr);    Hn1s = besselh(n,a.*rs);
    
    DHn1r = 0.50*(besselh(n - 1, a*rr) - besselh(n + 1, a*rr));
    DHn1s = 0.50*(besselh(n - 1, a*rs) - besselh(n + 1, a*rs));

    %%%%-----------------------------------------------------------------------
    DJnb = 0.50*(besselj(n - 1, b*rc) - besselj(n + 1, b*rc));
    Jnb  =  besselj(n, b*rc);

    DJna = 0.50*(besselj(n - 1, a*rc) - besselj(n + 1, a*rc));
    Jna = besselj(n,a.*rc);

    DHna = 0.50*(besselh(n - 1, a*rc) - besselh(n + 1, a*rc));
    Hna = besselh(n, a*rc);

    Det = rc^2*(k2^2*DJnb./(b.*Jnb) - k1^2*DHna./(a.*Hna)).*...
        (DJnb./(b.*Jnb) - DHna./(a.*Hna)) - (n.^2).*(x.^2).*((b.^(-2)-a.^(-2)).^2);
    % y = ( ( (k2^2*DJn(n,b*rc))./(b.*besselj(n,b*rc)) - (k1^2*DHn(n,a*rc))./...
    %     (a.*besselh(n,a*rc)) ).*( (DJn(n,b*rc))./(b.*besselj(n,b*rc)) - (DHn(n,a*rc))./...
    %     (a.*besselh(n,a*rc)) ) )*rc^2 - n^2*(k1^2-a.^2).*(( b.^(-2) - a.^(-2) ).^2);
    %%%-----------------------------------------------------------------------

%%% -------------------- FRESNEL COEFFICIENTS -----------------------------
    Rn11mm = (-1)*Jna./Hna./Det.*((k2^2*DJnb./(b.*Jnb) - k1^2*DHna./(a.*Hna)).*...
    (DJnb./(b.*Jnb) - DJna./(a.*Jna)).*rc^2 - n.^2.*(x.^2).*((a.^(-2)-b.^(-2)).^2));
    Rn11mn = k1.*n.*rc.*x.*Jna./a./Hna./Det.*(a.^(-2)-b.^(-2)).*(DJna./Jna - DHna./Hna);
    Rn11nm = k1.*n.*rc.*x.*Jna./a./Hna./Det.*(a.^(-2)-b.^(-2)).*(DJna./Jna - DHna./Hna);
    Rn11nn = Jna./Hna.*((b.^(-2) - a.^(-2)).^2.*n.^2.*x.^2 - ...
    (DJnb./(Jnb.*b) - DHna./(Hna.*a)).*(DJnb.*k2.^2./(Jnb.*b) - DJna.*k1.^2./(Jna.*a)).*rc.^2)...
    ./(-(b.^(-2) - a.^(-2)).^2.*n.^2.*x.^2 + (DJnb./(Jnb.*b) - DHna./(Hna.*a)).*(DJnb.*k2.^2./(Jnb.*b) - DHna.*k1.^2./(Hna.*a)).*rc.^2);
%%% -----------------------------------------------------------------------

    %% rr component
    if((i == 1) && (j == 1))

          iGNrr11mm = Hn1r.*Hn1s.*n.^2.*Rn11mm./(rr.*rs.*a.^2);
          iGNrr11nm = DHn1r.*Hn1s.*n.*Rn11nm.*x./(k1.*rs.*a);
          iGNrr11mn = DHn1s.*Hn1r.*n.*Rn11mn.*x./(k1.*rr.*a);
          iGNrr11nn = DHn1r.*DHn1s.*Rn11nn.*x.^2./(k1.^2);
          
          y = (2 - double(0==n)).*1i.*cos(n.*(pr-ps)).*(iGNrr11mm + iGNrr11nm + iGNrr11mn + ...
              iGNrr11nn).*exp(1i.*x.*(zr-zs))./(8.*pi);
    end
    
    %% rp component
    if((i == 1) && (j == 2))

        iGNrp11mm = n.*Hn1r.*DHn1s.*Rn11mm./(rr.*a);
        iGNrp11nm = x.*DHn1r.*DHn1s.*Rn11nm./k1;
        iGNrp11mn = Hn1r.*Hn1s.*n.^2.*x.*Rn11mn./k1./rr./rs./a.^2;
        iGNrp11nn = DHn1r.*Hn1s.*n.*Rn11nn.*x.^2./(k1.^2.*rs.*a);

         y = (2-double(0==n))*1i*(iGNrp11mm+iGNrp11nm+...
            iGNrp11mn+iGNrp11nn).*sin(n.*(pr-ps)).*exp(1i*x.*(zr-zs))./...
            (8*pi);
    end
    
    %% rz component
    if((i == 1) && (j == 3))

        iGNrz11mm = 0.0;
        iGNrz11nm = 0.0;
        iGNrz11mn = 1i.*Hn1r.*Hn1s.*n.*Rn11mn./k1./rr;
        iGNrz11nn = 1i.*a.*DHn1r.*Hn1s.*Rn11nn.*x./k1./k1;

         y = (2-double(0==n))*1i*(iGNrz11mm+iGNrz11nm+...
            iGNrz11mn+iGNrz11nn).*cos(n.*(pr-ps)).*exp(1i*x.*(zr-zs))./...
            (8*pi);
    end

    %% pr component
    if((i == 2) && (j == 1))
        
        iGNpr11mm = -DHn1r.*Hn1s.*n.*Rn11mm./(rs.*a);
        iGNpr11nm = -Hn1r.*Hn1s.*n.^2.*Rn11nm.*x./(k1.*rr.*rs.*a.^2);
        iGNpr11mn = -DHn1r.*DHn1s.*Rn11mn.*x./k1;
        iGNpr11nn = -DHn1s.*Hn1r.*n.*Rn11nn.*(x.^2)./(k1.^2.*rr.*a);

          y = (2-double(0==n))*1i*(iGNpr11mm+iGNpr11nm+...
          iGNpr11mn+iGNpr11nn).*sin(n.*(pr-ps)).*exp(1i*x.*(zr-zs))./...
          (8*pi);

    end
    
    %% pp component
    if((i == 2) && (j == 2))

        iGNpp11mm = Rn11mm.*DHn1r.*DHn1s;
        iGNpp11nm = n.*x.*Rn11nm.*Hn1r.*DHn1s./(k1.*rr.*a);
        iGNpp11mn = n.*Rn11mn.*DHn1r.*Hn1s.*x./(k1.*rs.*a);
        iGNpp11nn = (n.^2).*(x.^2).*Rn11nn.*Hn1r.*Hn1s./(k1.^2.*rr.*rs.*a.^2);

          y = (2-double(0==n))*1i*(iGNpp11mm+iGNpp11nm+...
          iGNpp11mn+iGNpp11nn).*cos(n.*(pr-ps)).*exp(1i*x.*(zr-zs))./...
          (8*pi);

    end

    %% pz component
    if((i == 2) && (j == 3))

        iGNpz11mm = 0.0;
        iGNpz11nm = 0.0;
        iGNpz11mn = -1i*a.*DHn1r.*Hn1s.*Rn11mn./k1;
        iGNpz11nn = -1i*Hn1r.*Hn1s.*n.*Rn11nn.*x./(k1.^2.*rr);

          y = (2-double(0==n))*1i*(iGNpz11mm+iGNpz11nm+...
          iGNpz11mn+iGNpz11nn).*sin(n.*(pr-ps)).*exp(1i*x.*(zr-zs))./...
          (8*pi);

    end
    
    %% zr component
    if((i == 3) && (j == 1))

        iGNzr11mm = 0.0;
        iGNzr11nm = -1i.*Hn1r.*Hn1s.*n.*Rn11nm./(k1.*rs);
        iGNzr11mn = 0.0;
        iGNzr11nn = -1i.*a.*DHn1s.*Hn1r.*Rn11nn.*x./(k1.^2);

          y = (2-double(0==n))*1i*(iGNzr11mm+iGNzr11nm+...
          iGNzr11mn+iGNzr11nn).*cos(n.*(pr-ps)).*exp(1i*x.*(zr-zs))./...
          (8*pi);

    end
    
    %% zp component
    if((i == 3) && (j == 2))

        iGNzp11mm = 0.0;
        iGNzp11nm = -1i*a.*DHn1s.*Hn1r.*Rn11nm./k1;
        iGNzp11mn = 0.0;
        iGNzp11nn = -1i*Hn1r.*Hn1s.*n.*Rn11nn.*x./(k1.^2.*rs);

          y = (2-double(0==n))*1i*(iGNzp11mm+iGNzp11nm+...
          iGNzp11mn+iGNzp11nn).*sin(n.*(pr-ps)).*exp(1i*x.*(zr-zs))./...
          (8*pi);

    end

    %% zz component
    if((i == 3) && (j == 3))

        y = (2-double(0==n))*1i*Rn11nn.*Hn1r.*...
        ((Hn1s).*(a.^2)).*cos(n*(pr-ps)).*exp(1i*x.*(zr-zs))...
            ./(k1.^2*8*pi);
    end

end