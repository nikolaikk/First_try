%%          Beam Propagation Method

zo = 1;
zend = 1000;    % micrometer
z_mesh = 1;    % micrometer
z = zo:z_mesh:zend;
Nz = length(z);

Lambda = 0.8;   % micrometer
no = 1;
ko = 2*pi/Lambda*no;
n = 1;
k = 2*pi/Lambda*n;

xo = -60; %mvim
xend = -xo; %m
Nx = 101;
x_mesh = (xend - xo)/Nx;
x = linspace(xo,xend-x_mesh,Nx);
Nx = length(x);
wo = 8; % micrometer
f = zeros(Nx);
f(:,1) = exp(-(x/wo).^2);


%% Define L matrix
L = zeros(Nx);
for i=1:Nx
    L(i,i)=-2+(ko^2-k^2)*x_mesh^2; % in free space^2 - in material^2
end
for i=2:(Nx)
    L(i,i-1)=1;
    L(i-1,i)=1;
end
L = (1/x_mesh^2)*L;
%%
E = eye(Nx);


%% Psi evolution

const_to_simplify = (pinv(E-x_mesh*L/(4*1j*ko))*(E+x_mesh*L/(4*1j*ko)));
% const_to_simplify(find(const_to_simplify<1e-14))=0
for i=1:Nz-1
    f(:,i+1)= const_to_simplify*f(:,i);
%     disp(i)
end


[Z, X] = meshgrid(z,x);
I = conj(f).*f;

figure('Name','Numeric Gaussian Beam','NumberTitle','off');
mesh(Z,X,(I));
axis tight;
shading interp;
colormap jet;
xlabel ('z (\mum)');
ylabel ('x (\mum)');
zlabel Intensity;
rotate3d on

% I_analytic = function_analytic_BMP2(Lambda,wo,x,z,true);
% 
% figure('Name','Intensity Difference','NumberTitle','off');
% mesh(Z,X,((I-I_analytic).^2));
% axis tight;
% shading interp;
% xlabel ('z (\mum)');
% ylabel ('x (\mum)');
% zlabel Intensity;
% rotate3d on
% zlim([0,1]);
% xlim([0,zend]);
% ylim([xo,-xo]);
% 



