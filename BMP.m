clear
%%              Beam Propagaton Method
    
%%              Initialization

zo = 0;                     % micrometer
zend = 1200;                % micrometer
z_mesh = 1;                 % micrometer
ratio = 1;
xo = -200;                  % micrometer
xend = -xo;                 % micrometer
x_mesh = z_mesh/ratio;      % micrometer
Lambda = 10;              % micrometer
wo = 40;                    % micrometer
no = 1;

[f,x,z] = solve_BPM(zo, zend, z_mesh, xo, xend, x_mesh, Lambda, wo, no);

%%              Convergence
 
    step = 0.002;
    err = [];
    err1 = err;
    err2 = err;
    % Number_Nz = ones(round((zend-zo)/step)+1,1)';
    Number_Nz = [];
    ii = 1;

    tic
    from = 1;
    to = 0.5;
    for i_z_mesh = z_mesh*(from:-step:to)
%         disp(i_z_mesh)
%         x_mesh = i_z_mesh/ratio;      % micrometer
        [f,x,z] = solve_BPM(zo, zend, i_z_mesh, xo, xend, x_mesh, Lambda, wo, no);

        I = conj(f).*f;
%         I_analytic = function_analytic_BMP(Lambda,wo,x,z);
        I_analytic = function_analytic_BMP(zo, zend, i_z_mesh, xo, xend, x_mesh, Lambda, wo, no, 'alll',false);

        Number_Nz(ii) = length(z);
        
%         err1(ii) = 1/Number_Nz(ii)*norm(I-I_analytic,2);
%         err2(ii) = 1/Number_Nz(ii)*norm(I-I_analytic,inf);
        
        err1(ii) = 1/Number_Nz(ii)*sqrt(sum(sum((I-I_analytic).^2)));
        err2(ii) = 1/Number_Nz(ii)*max(max(I-I_analytic));

        clc;
        disp(strcat(num2str(ii),' /  ', num2str(round((from-to)/step))));
        ii=ii+1;
    end
    toc

    subplot(2,1,1),plot(Number_Nz,err1),title('L2 norm error')
    xlabel('Number of grid points'),ylabel('Error')
    set(gca, 'XScale', 'log'),set(gca, 'YScale', 'log')

    subplot(2,1,2), plot(Number_Nz,err2),title('Lmax norm error')
    xlabel('Number of grid points'),ylabel('Error')
    set(gca, 'XScale', 'log'),set(gca, 'YScale', 'log')

%%              Plot


% [Z, X] = meshgrid(z,x);
% I = conj(f).*f;
% % I= real(f);
% 
% % Plot 
% fig1 = figure;
% set(fig1,'Position',[0,0,800,1200])
% subplot(3,1,1)
% mesh(Z,X,(I)),view([0,90])
% axis tight, shading interp, xlabel ('z (\mum)'), ylabel ('x (\mum)'), zlabel Intensity,...
%     rotate3d on, title('Crank–Nicolson BPM'), colormap jet;
% 
% subplot(3,1,2)
% I_analytic = function_analytic_BMP(zo, zend, z_mesh, xo, xend, x_mesh, Lambda, wo, no,'real',true);
% 
% subplot(3,1,3)
% % I_FFT = function_BMP_FFT(zo, zend, z_mesh, xo, xend, x_mesh, Lambda, wo, no,true);
% I_analytic = function_analytic_BMP(zo, zend, z_mesh, xo, xend, x_mesh, Lambda, wo, no,'alll',true);
% 
% % subplot(3,1,3)
% % mesh(Z,X,(I-I_analytic)),view([0,90])
% % axis tight, shading interp, xlabel ('z (\mum)'), ylabel ('x (\mum)');
% % zlabel Intensity, rotate3d on, zlim([0,1]), xlim([0,zend]), ylim([xo,-xo]),...
% %     title('Intensity Difference');
% 
% % disp(norm(I-I_analytic))
% 
% fig2 = copyobj(fig1,0);
% set(fig2, 'Position', [860,0,800,1200])
% subplot(3,1,1),view([0,-90,0]);
% subplot(3,1,2),view([0,-90,0]);
% subplot(3,1,3),view([0,-90,0]);
% 
% divergence_angle = theta_angle(I_analytic,Lambda,'plot');



%%              Main BPM funciton 
function [f,x,z] = solve_BPM(zo, zend, z_mesh, xo, xend, x_mesh, Lambda, wo, no)

x = xo:x_mesh:xend-x_mesh;
z = zo:z_mesh:zend-z_mesh;
Nz = length(z);
Nx = length(x);

ko = 2*pi/Lambda*no;
n = 1;
k = 2*pi/Lambda*n;
f = zeros(Nx);
f(:,1) = exp(-(x/wo).^2);

%%              Definition of L matrix

L = zeros(Nx);
for i=1:Nx
    L(i,i)=-2+(ko^2-k^2)*x_mesh^2; % in free space^2 - in material^2
end
for i=2:(Nx)
    L(i,i-1)=1;
    L(i-1,i)=1;
end
L = (1/x_mesh^2)*L;

%%              Psi evolution

E = eye(Nx);
const_to_simplify = (pinv(E-x_mesh*L/(4*1j*ko))*(E+x_mesh*L/(4*1j*ko)));
% const_to_simplify(find(const_to_simplify<1e-14))=0
for i=1:Nz-1
% for i=1:2

f(:,i+1)= const_to_simplify*f(:,i);
% disp(i)
end
end