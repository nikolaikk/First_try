clear
%%              Beam Propagaton Method
    
%%              Initialization

zo = 0;                     % micrometer
zend = 400;                % micrometer
z_mesh = 1;                 % micrometer
ratio = 1;
xo = -400;                  % micrometer
xend = -xo;                 % micrometer
x_mesh = z_mesh/ratio;      % micrometer
Lambda = 10;              % micrometer
wo = 30;                    % micrometer
no = 1;

[Ez,x,z,Nz,Nx,L] = solve_BPM(zo, zend, z_mesh, xo, xend, x_mesh, Lambda, wo, no);


%%              Start calculation
part_of_program = 0;
part_of_program = 1;

record_status = 0;
% record_status = 1;

movie_name = 'animate5.avi';

%%              Convergence
 
if part_of_program  == 0
    
    step = 0.01;
    err_l2 = [];
    err_max = [];
    err_l2_point = [];
    Number_Nz = [];
    ii = 1;

    tic
    from = 1;
    to = 0.7;
    range_z = from:-step:to;
    range_x = linspace(x_mesh,x_mesh-0.25,length(range_z));
    
    % Video
    if record_status == 1
        vid = VideoWriter(movie_name);
        vid.FrameRate = 6;
        vid.open()
    end
%     
    for i_z_mesh = z_mesh*(range_z)
        disp(i_z_mesh)

        x_mesh = range_x(ii);      % micrometer
        [f,x,z] = solve_BPM(zo, zend, i_z_mesh, xo, xend, x_mesh, Lambda, wo, no);

        I = conj(f).*f;
%         I_analytic = function_analytic_BMP(Lambda,wo,x,z);
        I_analytic = function_analytic_BMP(zo, zend, i_z_mesh, xo, xend, x_mesh, Lambda, wo, no,'real');

        Number_Nz(ii) = length(z);
        Number_Nz_Nx(ii) = length(z)*length(x);
        
%         err1(ii) = 1/Number_Nz(ii)*norm(I-I_analytic,2);
%         err2(ii) = 1/Number_Nz(ii)*norm(I-I_analytic,inf);
        
        I_diff = (I-I_analytic);
        err_l2(ii) = 1/Number_Nz_Nx(ii)*sqrt(sum(sum(I_diff).^2));
        err_max(ii) = max(max(I_diff));
        
        px = 0.2;
        py = 0.5;
        err_l2_point(ii) = abs(I_diff(round(py*size(I,1)),round(px*size(I,2)))); % change of error at point (px,py)

        clc;
        disp(strcat(num2str(ii),' /  ', num2str(round((from-to)/step)+1)));
        ii=ii+1;
        
        if record_status == 1
%             imshow(imresize(I,[800,1600]))
            contour((imresize(I,[800,1600])),25), set(gcf,'Position',[0,0,800,1600]);
            vid.writeVideo(getframe(gca));
        end
    end
    
    if record_status == 1
        vid.close()
    end
    toc

    fig = figure,set(fig,'Position',[0,0,900,1200]);
    
    subplot(2,1,1),plot(Number_Nz,err_l2),title('L2 norm error')
    xlabel('Number of grid points'),ylabel('Error')
    set(gca, 'XScale', 'log'),set(gca, 'YScale', 'log')
    hold on% f(1,i+1) = 0;
% f(end,i+1) = 0;

    
    i_max=(find(err_l2==max(err_l2)));
    Xx = [Number_Nz(i_max),Number_Nz(i_max)*2,Number_Nz(i_max)*4,Number_Nz(i_max)*8,Number_Nz(i_max)*16];
    Yy = [err_l2(i_max),err_l2(i_max)/2,err_l2(i_max)/4,err_l2(i_max)/8,err_l2(i_max)/16];
    
    
    
    subplot(2,1,2), plot(Number_Nz,err_max),title('Lmax norm error')
    xlabel('Number of grid points'),ylabel('Error')
    set(gca, 'XScale', 'log'),set(gca, 'YScale', 'log')

    % Plot error at particular point
    
%     subplot(3,1,3), plot(Number_Nz,err_l2_point),title('L2 point')
%     xlabel('Number of grid points'),ylabel('Error')
%     set(gca, 'XScale', 'log'),set(gca, 'YScale', 'log')
    
%     save('worlspace1II.mat')

%%              Plot

else
    
[Z, X] = meshgrid(z,x);
I = conj(Ez).*Ez;
% I= real(f);

% Plot 
fig1 = figure;
set(fig1,'Position',[0,0,1600,1200])

subplot(2,2,1)
% mesh(Z,X,(I));
mesh(Z,X,(real(Ez)));
view([0,90])
axis tight, shading interp, xlabel ('z (\mum)'), ylabel ('x (\mum)'), zlabel ('Real Part Electric Field'),...
    rotate3d on, title('Re{Ez} Crank–Nicolson BPM'), colormap jet,colorbar;

subplot(2,2,2)
mesh(Z,X,(I));
view([0,90])
axis tight, shading interp, xlabel ('z (\mum)'), ylabel ('x (\mum)'), zlabel Intensity,...
    rotate3d on, title('Intensity Crank–Nicolson BPM'), colormap jet,colorbar;

subplot(2,2,3)
I_analytic = function_analytic_BMP(zo, zend, z_mesh, xo, xend, x_mesh, Lambda, wo, no,'real');

subplot(2,2,4)
Ez_analytic = function_analytic_BMP(zo, zend, z_mesh, xo, xend, x_mesh, Lambda, wo, no,'I');


% subplot(3,1,3)
% I_FFT = function_BMP_FFT(zo, zend, z_mesh, xo, xend, x_mesh, Lambda, wo, no,true);

% subplot(3,1,3)
% mesh(Z,X,(I-I_analytic)),view([0,90])
% axis tight, shading interp, xlabel ('z (\mum)'), ylabel ('x (\mum)');
% zlabel Intensity, rotate3d on, zlim([0,1]), xlim([0,zend]), ylim([xo,-xo]),...
%     title('Intensity Difference');

% disp(norm(I-I_analytic))

% fig2 = copyobj(fig1,0);
% set(fig2, 'Position', [860,0,800,1200])
% subplot(3,1,1),view([0,-90,0]);
% subplot(3,1,2),view([0,-90,0]);
% subplot(3,1,3),view([0,-90,0]);

% divergence_angle = theta_angle(I_analytic,z_mesh,x_mesh,Lambda,'plot');
end



%%              Main BPM funciton 
function [f,x,z,Nz,Nx,L] = solve_BPM(zo, zend, z_mesh, xo, xend, x_mesh, Lambda, wo, no)

x = xo:x_mesh:xend-x_mesh;
z = zo:z_mesh:zend-z_mesh;
Nz = length(z);
Nx = length(x);

ko = 2*pi/Lambda*no;
n = 1;
k = 2*pi/Lambda*n;
f = zeros(Nx,Nz);
f(:,1) = exp(-(x/wo).^2);

f(1,1) = 0;
f(end,1) = 0;
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

f(:,i+1)= const_to_simplify*f(:,i);
% f(:,i+1) = f(:,i+1)*cos(ko*z) 

f(1,i+1) = 0;
f(end,i+1) = 0;

end

extra = repmat(exp(-1j*ko*z),800,1);
f = f.*extra;

end