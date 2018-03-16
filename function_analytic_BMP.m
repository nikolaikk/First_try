function [I,Ez] = function_analytic_BMP(zo, zend, z_mesh, xo, xend, x_mesh, Lambda, wo, no, type)
    
    x = xo:x_mesh:xend-x_mesh;
    z = zo:z_mesh:zend-z_mesh;
    
    E = 1;
    k = 2* pi/Lambda/no;
    zR = pi*wo^2/Lambda;
    Nx = length(x);
    Nz = length(z);
    I = zeros([Nx,Nz]);
    Ez = zeros([Nx,Nz]);
    
    for i=1:Nz
        Ez(:,i) = E* wo./width(z(i),zR,wo).*exp(-(x.^2)./width(z(i),zR,wo).^2).*exp(1j*(k*z(i)+k*(x.^2)./(2*R(z(i),zR))-psi(z(i))));
%         Ez(:,i) = E* wo./width(z(i),zR,wo).*exp(-(x.^2)./width(z(i),zR,wo).^2).*exp(1j*(k*z(i)+k*(x.^2)./(2*R(z(i),zR))));
        I(:,i) =((wo./width(z(i),zR,wo)).^2).*exp(-2*(x.^2)/(width(z(i),zR,wo).^2));
%         I(:,i) = exp(-2*(x.^2)/(width(z(i),zR,wo).^2));       % peak intensity is the same as beam propagates
    end
    
%     I = real(conj(Ez).* Ez);
%     I = I1-I;
        if type == 'imag'
            [Z, X] = meshgrid(z,x);
            mesh(Z,X,abs(real(Ez)));
            zlabel ('Imaginary Part Electric Field');
            title('Im{Ez} Analytic Gaussian Beam');

            view([0,90]), axis tight, shading interp, xlabel ('z (\mum)');
            ylabel ('x (\mum)'), colormap jet,rotate3d on, colorbar;

        elseif type == 'real'
            [Z, X] = meshgrid(z,x);
            mesh(Z,X,abs(imag(Ez)));
            zlabel ('Real Part Electric Field');
            title('Re{Ez} Analytic Gaussian Beam');
            
            view([0,90]), axis tight, shading interp, xlabel ('z (\mum)');
            ylabel ('x (\mum)'), colormap jet,rotate3d on, colorbar;

        elseif type == 'I'
            [Z, X] = meshgrid(z,x);
            mesh(Z,X,I);
            zlabel Intensity;
            title('Intensity Analytic Gaussian Beam');
            
            view([0,90]), axis tight, shading interp, xlabel ('z (\mum)');
            ylabel ('x (\mum)'), colormap jet,rotate3d on, colorbar;

        else
        end
      
    function w = width(z,zR,wo)
        w = wo*sqrt(1+(z/zR).^2);
    end
    function r = R(z,zR)
        r = z*(1+(zR/z)^2);
    end
    function gouy = psi(z)
        gouy = atan(z/zR);
    end

end