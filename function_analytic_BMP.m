function I = function_analytic_BMP(Lambda,Nz,wo,z_length,ratio,bla)

    if nargin < 6
        bla = false;
    end    
    
    k = 2* pi/Lambda;
    x = linspace(-0.001,0.001,Nz) * ratio;
    z = linspace(0,z_length,length(x)) * ratio;

    zR = pi*wo^2/Lambda;

    I = zeros([length(x),length(z)]);

    for i=1:length(z)
        I(:,i) = (wo./width(z(i),zR,wo)).*exp(-x.^2/width(z(i),zR,wo).^2);
    end
    
    if bla == true
        [Z, X] = meshgrid(z,x);
        figure,mesh(Z,X,(I));
        axis vis3d;
        shading interp;
        xlabel ('z (m)');
        ylabel ('x (m)');
        zlabel I;
        rotate3d on
    end
    
    function t = width(z,zR,wo)
        t = sqrt(1+(z/zR)^2)*wo;
    end
    function r = R(z,zR)
        r = z.*(1+(zR./z).^2);
    end
    function gouy = psi(z)
        gouy = atan(z/zR);
    end

end