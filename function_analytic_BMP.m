function I = function_analytic_BMP(zo, zend, z_mesh, xo, xend, x_mesh, Lambda, wo, no, type,bla)
    
    x = xo:x_mesh:xend-x_mesh;
    z = zo:z_mesh:zend-z_mesh;
    
    if nargin < 11
        bla = false;
    end    
    E = 1;
    k = 2* pi/Lambda/no;
    zR = pi*wo^2/Lambda;
    Nx = length(x);
    Nz = length(z);
    I = zeros([Nx,Nz]);
    Ez = zeros([Nx,Nz]);
    
    for i=1:Nz
%         Ez(:,i) = E* wo./width(z(i),zR,wo).*exp(-(x.^2+z(i)^2)./width(z(i),zR,wo).^2).*exp(1j*(k*z(i)+k*(x.^2+z(i)^2)./(2*R(z(i),zR))-psi(z(i))));
        Ez(:,i) = E* wo./width(z(i),zR,wo).*exp(-(x.^2)./width(z(i),zR,wo).^2).*exp(1j*(k*z(i)+k*(x.^2)./(2*R(z(i),zR))-psi(z(i))));
        
%         I(:,i) = ((wo./width(z(i),zR,wo)).^2).*exp(-(x.^2+z(i)^2)./(width(z(i),zR,wo).^2));
        I(:,i) = ((wo./width(z(i),zR,wo)).^2).*exp(-(x.^2)/(width(z(i),zR,wo).^2));
    end
    
%     I = real(conj(Ez).* Ez);
%     I = I1-I;
        
    if bla == true  
        if type == 'real'
            [Z, X] = meshgrid(z,x);
            mesh(Z,X,abs(real(Ez)));
            zlabel ('Imaginary Part Electric Field');
        elseif type == 'imag'
            [Z, X] = meshgrid(z,x);
            mesh(Z,X,abs(imag(Ez)));
             zlabel ('Real Part Electric Field');
        else
            [Z, X] = meshgrid(z,x);
            mesh(Z,X,I);
             zlabel Intensity;
        end
        view([0,90]);
        title('Analytic Gaussian Beam'),
        axis tight;
        shading interp;
        xlabel ('z (\mum)');
        ylabel ('x (\mum)');
        colormap jet;
        rotate3d on
    end
        
    function w = width(z,zR,wo)
        w = sqrt(1+(z/zR).^2)*wo;
    end
    function r = R(z,zR)
        r = z.*(1+(zR./z).^2);
    end
    function gouy = psi(z)
        gouy = atan(z/zR);
    end

end