function [I,Ez] = function_analytic_BMP(zo, zend, z_mesh, xo, xend, x_mesh, Lambda, wo, no, type)
    
    x = xo:x_mesh:xend-x_mesh;
    z = zo:z_mesh:zend-z_mesh;
    
    E = 1;
    Io = 1;
    k = 2* pi/Lambda*no;
    zR = pi*wo^2/Lambda;
    Nx = length(x);
    Nz = length(z);
    I = zeros([Nx,Nz]);
    Ez = zeros([Nx,Nz]);
    
    for i=1:Nz
        
        % Intensity Normalized
%         Ez(:,i) = E* (2/pi/width(z(i),zR,wo).^2).^0.25.*exp(-(x.^2)./width(z(i),zR,wo).^2).*exp(-1j*(k*z(i)+k*(x.^2)./(2*R(z(i),zR))-psi(z(i),zR)/2));
        
        % Amplitude Normalized
        Ez(:,i) = E * (wo./width(z(i),zR,wo)).^0.5.*exp(-(x.^2)./width(z(i),zR,wo).^2).*exp(1j*(k*z(i)+k*(x.^2)./(2*R(z(i),zR))-psi(z(i),zR)/2));
   
    end
    
    I = abs(Ez).^2;
    
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
%       plot(z,arr)
        
        
    function w = width(z,zR,wo)
        w = wo*sqrt(1+(z/zR).^2);
    end
    function r = R(z,zR)
        r = z+(zR^2/z);
    end
    function gouy = psi(z,zR)
        gouy = atan(z/zR);
    end

end