function [I,Ez] = function_analytic_BMP(zo, zend, z_mesh, xo, xend, x_mesh, Lambda, wo, no)
    
    x = xo:x_mesh:xend-x_mesh;
    z = zo:z_mesh:zend-z_mesh;
    E = 1;
    k = 2* pi/Lambda*no;
    zR = pi*wo^2/Lambda;
    Nx = length(x);
    Nz = length(z);
    Ez = zeros([Nx,Nz]);
    epsilon = 8.854187817e-12;
    c = 299792458;
    
    
    for i=1:Nz
        w = wo*sqrt(1+(z(i)/zR).^2);
        r = z(i)+(zR^2/z(i));
        gouy = atan(z(i)/zR);

        % Amplitude Normalized
        Ez(:,i) = E * (wo./w).^0.5.*exp(-(x.^2)./w^2).*exp(1j*(k*z(i)+k*(x.^2)./(2*r)-gouy));
   
    end
    
    I = epsilon*no*c*0.5*abs(Ez).^2;
end