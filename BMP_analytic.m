ratio = 1;
Lambda = 600e-9;   %meter
wo=0.0001;
k = 2* pi/Lambda;
Nz = 20
x = linspace(-0.001,0.001,Nz) * ratio;
z = linspace(0,0.3,length(x)) * ratio;

zR = pi*wo^2/Lambda;



% x = width(3,zR,wo)
% R(1:5,zR)

I = zeros([length(x),length(z)]);

for i=1:length(z)
    I(:,i) = (wo./width(z(i),zR,wo)).*exp(-x.^2/width(z(i),zR,wo).^2);
end
[Z, X] = meshgrid(z,x);
figure,mesh(Z,X,(I));
axis vis3d;
shading interp;
xlabel ('z (m)');
ylabel ('x (m)');
zlabel I;
rotate3d on
toc



function t = width(z,zR,wo)
    t = sqrt(1+(z/zR)^2)*wo;
end
function r = R(z,zR)
    r = z.*(1+(zR./z).^2);
end
function gouy = psi(z)
    gouy = atan(z/zR);
end