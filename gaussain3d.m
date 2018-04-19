clear;
zo = 0;                     % micrometer
zend = 400;                % micrometer
z_mesh = 1;                 % micrometer
ratio = 1;
xo = -200;                  % micrometer
xend = -xo;                 % micrometer
x_mesh = z_mesh/ratio;      % micrometer
Lambda = 10;                % micrometer
wo = 30;                    % micrometer
no = 1;
type = 'I';

x = xo:x_mesh:xend-x_mesh;
z = zo:z_mesh:zend-z_mesh;

E = 1;
Io = 1;
k = 2* pi/Lambda*no;
zRR = pi*wo^2/Lambda;
Nx = length(x);
Nz = length(z);
I = zeros([Nx,Nx,Nz]);
Ez = zeros([Nx,Nx,Nz]);

for i=1:Nz
    for y=1:Nx
        Ez(:,y,i) = E * wo./width(z(i),zRR,wo).*exp(-(x.^2+x(y).^2)./width(z(i),zRR,wo).^2).*exp(1j*(k*z(i)+k*(x.^2+x(y).^2)./(2*R(z(i),zRR))-psi(z(i),zRR)));
        I(:,y,i) =Io*((wo./width(z(i),zRR,wo)).^2).*exp(-2*(x.^2+x(y).^2)/(width(z(i),zRR,wo).^2));
    end
end

[Z, X] = meshgrid(z,x);
if type == 'imag'
    mesh(Z,X,abs(real(reshape(Ez(:,round(Nx/2),:),[Nx,Nz]))));
    zlabel ('Imaginary Part Electric Field');
    title('Im{Ez} Analytic Gaussian Beam');

    view([0,90]), axis tight, shading interp, xlabel ('z (\mum)');
    ylabel ('x (\mum)'), colormap jet,rotate3d on, colorbar;

elseif type == 'real'
    mesh(Z,X,abs(imag(reshape(Ez(:,round(Nx/2),:),[Nx,Nz]))));
    zlabel ('Real Part Electric Field');
    title('Re{Ez} Analytic Gaussian Beam');

    view([0,90]), axis tight, shading interp, xlabel ('z (\mum)');
    ylabel ('x (\mum)'), colormap jet,rotate3d on, colorbar;

elseif type == 'I'
    mesh(Z,X,reshape(I(:,round(Nx/2),:),[Nx,Nz]));
    zlabel Intensity;
    title('Intensity Analytic Gaussian Beam');

    view([0,90]), axis tight, shading interp, xlabel ('z (\mum)');
    ylabel ('x (\mum)'), colormap jet,rotate3d on, colorbar;

else
end
%       plot(z,arr)


% save to VTK 
saveVTK('TestFile.vtk',I)

% check conservation of power
sum(sum(I(:,:,1)))
sum(sum(I(:,:,50)))
sum(sum(I(:,:,400)))



function w = width(z,zR,wo)
    w = wo*sqrt(1+(z/zR).^2);
end
function r = R(z,zR)
    r = z*(1+(zR/z)^2);
end
function gouy = psi(z,zR)
    gouy = atan(z/zR);
end