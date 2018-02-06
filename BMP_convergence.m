from = 1.0;
to = 0.5;
step = 0.001;
err = ones((to-from)/step,1)';
err1 = err;
err2 = err;
Number_Nz = ones(round((from-to)/step)+1,1)';
ii = 1;

zo = 0;
zend = 1000;    % micrometer


for z_mesh=from:-step:to
     
z = zo:z_mesh:zend-z_mesh;
Nz = length(z);

Lambda = 1.5;   % micrometer
no = 1;
ko = 2*pi/Lambda*no;
n = 1;
k = 2*pi/Lambda*n;

xo = -100; %mvim
xend = -xo; %m
x_mesh = z_mesh;
x = xo:x_mesh:xend-x_mesh;
Nx = length(x);
wo = 10; %m
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
I_analytic = function_analytic_BMP2(Lambda,wo,x,z);


err(ii) = sqrt(mean2((I-I_analytic).^2));
err1(ii) = norm(I-I_analytic);
err2(ii) = norm(I-I_analytic,1);


Number_Nz(ii) = Nz;

disp(strcat(num2str(ii),' / ', num2str(length(Number_Nz))));
ii=ii+1;

end

figure, plot(Number_Nz,err1)
figure, plot(Number_Nz,err2)

% figure,mesh(Z,X,(I));
% axis vis3d;
% shading interp;
% xlabel ('z (m)');
% ylabel ('x (m)');
% zlabel I;
% rotate3d on
% 
% figure,mesh(Z,X,((I-I_analytic).^2));
% axis normal;
% shading interp;
% xlabel ('z (m)');
% ylabel ('x (m)');
% zlabel I;
% zlim([0,1]);
% xlim([0,z_length]);
% ylim([-xo,xo]);
% rotate3d on
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
xlabel('Number of grid points')
ylabel('Error')

