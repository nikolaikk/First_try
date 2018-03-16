tic

zo = 0;
zend = 0.08; %m
Nz = 600; 
delta_z = (zend-zo)/(Nz-1); % step size
z=zo:delta_z:zend;
% z = linspace(zo,zend,Nz)
% z=z(1:100);

Lambda= 500e-9;
no=1;
ko=2*pi/Lambda*no;
n=1;
k=2*pi/Lambda*n;

xo=-0.001; %mvim
xend=-xo; %m
delta_x=(xend-xo)/(Nz-1);

x=xo:delta_x:xend;
wo=0.0001; %m
f=zeros(Nz);
f(:,1)=exp(-(x/wo).^2);

% 
% %% ANALYTIC
% zR = pi*wo^2/Lambda;
% Ia = zeros([length(x),length(z)]);
% 
% for i=1:length(z)
%     Ia(:,i) = (wo./width(z(i),zR,wo)).*exp(-x.^2/width(z(i),zR,wo).^2);
% end
%% NUMERIC
% Define L matrix
L = zeros(Nz);
for i=1:Nz
    L(i,i)=-2+(ko^2-k^2)*delta_x^2; % in free space^2 - in material^2
    end
for i=2:(Nz)
    L(i,i-1)=1;
    L(i-1,i)=1;
    end
L = 1/delta_x^2*L;
%
E = eye(Nz);

% Psi evolution
for i=1:Nz-1
% for i=1:2
    const_to_simplify = ((E-delta_z*L/(4*1j*ko))*pinv(E+delta_z*L/(4*1j*ko)));
    f(:,i+1)= const_to_simplify*f(:,i);

end

%% PLOT 3D
[Z, X] = meshgrid(delta_z*(1:Nz),delta_x*(1:Nz));
I = conj(f).*f;
figure,mesh(Z,X,I);
axis vis3d;
shading interp;
xlabel ('z (m)');
ylabel ('x (m)');
zlabel I;
rotate3d on
toc


%% FUNCTIONS 

function t = width(z,zR,wo)
    t = sqrt(1+(z/zR)^2)*wo;
end
function r = R(z,zR)
    r = z.*(1+(zR./z).^2);
end
function gouy = psi(z)
    gouy = atan(z/zR);
end