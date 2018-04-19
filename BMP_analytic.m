zo = 0;                     % micrometer
zend = 400;                % micrometer
z_mesh = 1;                 % micrometer
ratio = 1;
xo = -200;                  % micrometer
xend = -xo;                 % micrometer
x_mesh = z_mesh/ratio;      % micrometer
Lambda = 10;              % micrometer
wo = 30;                    % micrometer
no = 1;

x = xo:x_mesh:xend-x_mesh;
z = zo:z_mesh:zend-z_mesh;
Nz = length(z);
Nx = length(x);
ko = 2*pi/Lambda*no;
n = 1;
k = 2*pi/Lambda*n*ones(1,length(x));
zR = pi*wo^2/Lambda;
I = zeros([length(x),length(z)]);
E = I;

for i=1:length(z)
    qo_z = R(z(i),zR) * pi*width(z(i),zR,wo)/(pi*width(z(i),zR,wo)-1j*Lambda*R(z(i),zR));
    E(:,i) = (1/qo_z) * exp(-1j*k.*(z(i)+x.^2./qo_z));
%     E(:,i) = E* wo./width(z(i),zR,wo).*exp(-(x.^2)./width(z(i),zR,wo).^2).*exp(1j*(k*z(i)+k*(x.^2)./(2*R(z(i),zR))-psi(z(i),zR)));
    I(:,i) = (wo./width(z(i),zR,wo)).*exp(-x.^2/width(z(i),zR,wo).^2);
end

[Z, X] = meshgrid(z,x);
% I = abs(E).^2;
I = conj(E).*E;
figure,mesh(Z,X,real(I));
axis vis3d;
shading interp;
xlabel ('z (m)');
ylabel ('x (m)');
zlabel I;
rotate3d on


function t = width(z,zR,wo)
    t = sqrt(1+(z/zR)^2)*wo;
end
function r = R(z,zR)
    r = z.*(1+(zR./z).^2);
end
function gouy = psi(z)
    gouy = atan(z/zR);
end