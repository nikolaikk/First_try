function I = function_analytic_BMP(Lambda,wo,x,z,bla)

    if nargin < 5
        bla = false;
    end    
    
    k = 2* pi/Lambda;
    zR = pi*wo^2/Lambda;
    Nx = length(x);
    Nz = length(z);
    I = zeros([Nx,Nz]);

    for i=1:Nz
        I(:,i) = ((wo./width(z(i),zR,wo)).^2).*exp(-2*x.^2/(width(z(i),zR,wo).^2));
    end
    
    if bla == true  
        [Z, X] = meshgrid(z,x);
        figure('Name','Analytic Gaussian Beam','NumberTitle','off'),mesh(Z,X,(I));
        axis tight;
        shading interp;
        xlabel ('z (\mum)');
        ylabel ('x (\mum)');
        zlabel Intensity;
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