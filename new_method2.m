clear
%%              Beam Propagaton Method
    
%%              Initialization

zo = 0;                     % micrometer
zend = 250;                % micrometer
z_mesh = 1;                 % micrometer
ratio = 1;
xo = -400;                  % micrometer
xend = -xo;                 % micrometer
x_mesh = z_mesh/ratio;      % micrometer
Lambda = 10;              % micrometer
wo = 30;                    % micrometer
no = 1;



%%              Start calculation

[Ez,x,z,Nz,Nx] = solve_BPM(zo, zend, z_mesh, xo, xend, x_mesh, Lambda, wo, no);


[Z, X] = meshgrid(z,x);
I = conj(Ez).*Ez;
contourf(X,Z,(I))
view([90,90])

%%              Main BPM funciton 
function [f,x,z,Nz,Nx] = solve_BPM(zo, zend, z_mesh, xo, xend, x_mesh, Lambda, wo, no)

x = xo:x_mesh:xend-x_mesh;
z = zo:z_mesh:zend-z_mesh;
Nz = length(z);
Nx = length(x);

ko = 2*pi/Lambda*no;
n = 1;
k = 2*pi/Lambda*n*ones(1,Nx);
f = zeros(Nx,Nz);
f(:,1) = exp(-(x/wo).^2);

f(1,1) = 0;
f(end,1) = 0;

%%              Absorbtion rigion
% rigion = 1/4;
% k_grad = ((min(x)+max(x)*rigion-x(1:round(length(x)/8)))*rigion/2);
% k_abs = zeros(1,length(x));
% k_abs(1:length(k_grad))= k_grad;
% k_abs(end-length(k_grad)+1:end) = flip(k_grad);
% % figure, plot(x,k_abs)
% n = n + 1j*k_abs;
% k = 2*pi/Lambda*n;
%%              Definition of L matrix

L = zeros(Nx);
for i=1:Nx
    L(i,i)=-2+(ko^2-k(i)^2)*x_mesh^2; % in free space^2 - in material^2
end
for i=2:(Nx)
    L(i,i-1)=1;
    L(i-1,i)=1;
end
L = (1/x_mesh^2)*L;

%%              Psi evolution

% E = eye(Nx);
% const_to_simplify = (pinv(E-x_mesh*L/(4*1j*ko))*(E+x_mesh*L/(4*1j*ko)));
% % const_to_simplify(find(const_to_simplify<1e-14))=0
% for i=1:Nz-1
% 
% f(:,i+1)= const_to_simplify*f(:,i);
% % f(:,i+1) = f(:,i+1)*cos(ko*z) 
% 
% f(1,i+1) = 0;
% f(end,i+1) = 0;
alpha = 0.5;
aj = -alpha/x_mesh^2;
bj = 2*alpha/x_mesh^2-alpha*(k(1)^2-ko^2)+2j*ko/z_mesh;
cj = -alpha/x_mesh^2;

rj = zeros(1,Nx)
for j=2:Nx-1
   rj(j) = (1-alpha)/x_mesh^2*(f(j-1,1)+f(j+1,1))+((1-alpha)*(k(j)^2+ko^2)-2*(1-alpha)/x_mesh^2+2j*ko/z_mesh)*f(j,1);
end

beta = bj;

for j=2:Nx
   gama = cj/beta;
   beta = bj-aj*gama;
   f(:,j)
end




extra = repmat(exp(-1j*ko*z),size(f,1),1);
% f = f.*extra;
end
