% clear
%%              Beam Propagaton Method
    
%%              Initialization

epsilon = 8.854187817e-12;
c = 299792458;

zo = 0;                     % micrometer
zend = 300;                % micrometer
z_mesh = 1;                 % micrometer
ratio = 1;
xo = -100;                  % micrometer
xend = -xo;                 % micrometer
x_mesh = z_mesh/ratio;      % micrometer
Lambda = 10;              % micrometer
wo = 10;                    % micrometer
no = 1;



%%              Start calculation
part_of_program = 0;
part_of_program = 1;

record_status = 0;
% record_status = 1;

movie_name = 'animate4.avi';

%%              Convergence
 
if part_of_program  == 0
    
    step = 0.1;
    err_l2 = [];
    err_max = [];
    err_point = [];
    Number_Nz = [];
    ii = 1;

    tic
    from = 4;
    to = 0.5;
    range_z = from:-step:to;
%     range_z = [4,2];
    range_z = [2,1*sqrt(2),1,sqrt(2)*0.5,0.5];

    % Video
    if record_status == 1
        vid = VideoWriter(movie_name);
        vid.FrameRate = 6;
        vid.open()
    end
%     
    for i_z_mesh = range_z
        disp(i_z_mesh)

        x_mesh = range_z(ii);      % micrometer
        
        [Ez,x,z] = solve_BPM(zo, zend, i_z_mesh, xo, xend, x_mesh, Lambda, wo, no);

        I = epsilon*no*c*0.5*abs(Ez).^2;
%         I_analytic = function_analytic_BMP(Lambda,wo,x,z);
        [I_analytic, Ez_analytic] = function_analytic_BMP(zo, zend, i_z_mesh,...
            xo, xend, x_mesh, Lambda, wo, no);

        Number_Nz(ii) = length(z);
        Number_Nz_Nx(ii) = length(z)*length(x);
        
        Ez_diff = abs(abs(real(Ez))-abs(real(Ez_analytic)));
        err_l2(ii) = 1/Number_Nz_Nx(ii) * sqrt(sum(sum(Ez_diff.^2)));
        err_max(ii) = max(max(Ez_diff));
        
        % Error at particular point
        xp = 0.7;     % from 0 to 1
        zp = 0.5;    % from 0 to 1
        xp_coord = round(length(x)*xp)+1;
        zp_coord = round(length(z)*zp)+1;
        err_point(ii) = abs(real(Ez_analytic(xp_coord,zp_coord)-Ez(xp_coord,zp_coord)));
        % --------------------------
         
        clc;
        disp(strcat(num2str(ii),' /  ', num2str(length(range_z))));
        ii=ii+1;
        
        % Record
        if record_status == 1
%             imshow(imresize(I,[800,1600]))
            Ez_diff(xp_coord,zp_coord)=0.0001;
            contour((imresize(Ez_diff,[800,1600])),25),hold on
            scatter(zp*1600,xp*800),hold off;
            colorbar, set(gcf,'Position',[0,0,800,1600]);
            vid.writeVideo(getframe(gca));
        end
    end
    
    % Record
    if record_status == 1
        vid.close()
    end
    toc

    fig = figure,set(fig,'Position',[0,0,900,1200]);
    
    subplot(3,1,1),plot(Number_Nz_Nx,err_l2,'o-'),title('L2 norm error')
    xlabel('Number of grid points'),ylabel('Error')
    set(gca, 'XScale', 'log'),set(gca, 'YScale', 'log')
    hold on% f(1,i+1) = 0;
% f(end,i+1) = 0;

     
    
    subplot(3,1,2), plot(Number_Nz_Nx,err_max,'o-')
    title('Lmax norm error')
    xlabel('Number of grid points'),ylabel('Error')
    set(gca, 'XScale', 'log'),set(gca, 'YScale', 'log')
    
    subplot(3,1,3), plot(Number_Nz_Nx,err_point,'o-')
    title(strcat('Error at physical point (',num2str(xp),', ',num2str(zp),')'))
    xlabel('Number of grid points'),ylabel('Error')
    set(gca, 'XScale', 'log'),set(gca, 'YScale', 'log')
    
    
    % Show error and zrange    
    fprintf(strcat('delta_z values  :   ',num2str(range_z),'\n'))
%     fprintf(strcat('L2 normed error values  :   ',num2str(err_l2),'\n'))
    fprintf(strcat('Number of grid points devided by max  :   ',...
        num2str(Number_Nz_Nx/Number_Nz_Nx(1)),'   so increases x',...
        num2str(round(Number_Nz_Nx(2)/Number_Nz_Nx(1))),'\n'))
    fprintf(strcat('L2 normed error values devided by max  :   ',...
        num2str(err_l2/err_l2(1)),'   error decreases x',...
        num2str(round(err_l2(1)/err_l2(2))),'\n'))
    fprintf(strcat('Lmax normed error values devided by max  :   ',...
        num2str(err_max/err_max(1)),'   error decreases x',...
        num2str(round(err_max(1)/err_max(2))),'\n'))
    fprintf(strcat('Error at certain point devided by max  :   ',...
        num2str(err_point/err_point(1)),'   error decreases x',...
        num2str(round(err_point(1)/err_point(2))),'\n'))    
    
    % Plot 2
    
%     figure,set(gcf,'Position',[0,0,900,1200]), subplot(2,1,1), contourf(real(Ez))
%     subplot(2,1,2), contourf(real(Ez_analytic))

    % Plot error at particular point
    
%     subplot(3,1,3), plot(Number_Nz,err_l2_point),title('L2 point')
%     xlabel('Number of grid points'),ylabel('Error')
%     set(gca, 'XScale', 'log'),set(gca, 'YScale', 'log')
    
%     save('worlspace1II.mat')

%%              Plot

else

[Ez,x,z,Nz,Nx,L] = solve_BPM(zo, zend, z_mesh, xo, xend, x_mesh, Lambda, wo, no);
[I_analytic,Ez_analytic] = function_analytic_BMP(zo, zend, z_mesh, xo, xend, ...
    x_mesh, Lambda, wo, no);


[Z, X] = meshgrid(z,x);
I = epsilon*no*c*0.5*abs(Ez).^2;
% I= real(f);

% Plot 
fig1 = figure;
set(fig1,'Position',[0,0,2000,1200])

subplot(2,2,1)
plot_mesh_function(Z,X,(abs(real(Ez))),'z (\mum)','x (\mum)',...
    'Real Part Electric Field','\Ree(E_{z}) Crank–Nicolson BPM')

subplot(2,2,2)
plot_mesh_function(Z,X,I,'z (\mum)','x (\mum)',...
    'Intensity','Intensity Crank–Nicolson BPM')

subplot(2,2,3)
plot_mesh_function(Z,X,abs(real(Ez_analytic)),'z (\mum)','x (\mum)',...
    'Real Part Electric Field','\Ree(E_{z}) Analytic Gaussian Beam')

subplot(2,2,4)
plot_mesh_function(Z,X,I_analytic,'z (\mum)','x (\mum)','Intensity',...
    'Intensity Analytic Gaussian Beam')





fig2 = figure;
set(fig2,'Position',[0,0,2000,600])

subplot(1,2,1)
plot_mesh_function(Z,X,(abs(real(Ez_abs))),'z (\mum)','x (\mum)',...
    'Real Part Electric Field','b')

subplot(1,2,2)
plot_mesh_function(Z,X,(abs(real(Ez))),'z (\mum)','x (\mum)',...
    'Real Part Electric Field','b')


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
function [f,x,z,Nz,Nx,L] = solve_BPM(zo, zend, z_mesh, xo, xend, x_mesh, ...
    Lambda, wo, no)

x = xo:x_mesh:xend-x_mesh;
z = zo:z_mesh:zend-z_mesh;
Nz = length(z);
Nx = length(x);

ko = 2*pi/Lambda*no;
n = 1;
k = 2*pi/Lambda*n*ones(1,length(x));
f = zeros(Nx,Nz);
f(:,1) = exp(-(x/wo).^2);

f(1,1) = 0;
f(end,1) = 0;

%%              Absorbtion rigion
rigion = 1/2;
k_grad = ((min(x)+max(x)*rigion-x(1:round(length(x)*rigion/2)))*rigion/2).^2;
k_abs = zeros(1,length(x));
k_abs(1:length(k_grad))= k_grad;
k_abs(end-length(k_grad)+1:end) = flip(k_grad);
% figure, plot(x,k_abs)
n = n + 0.0013j*k_abs;
k = 2*pi/Lambda*n;
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

E = eye(Nx);
const_to_simplify = (pinv(E-x_mesh*L/(4*1j*ko))*(E+x_mesh*L/(4*1j*ko)));
% const_to_simplify(find(const_to_simplify<1e-14))=0
for i=1:Nz-1

f(:,i+1)= const_to_simplify*f(:,i);
% f(:,i+1) = f(:,i+1)*cos(ko*z) 

f(1,i+1) = 0;
f(end,i+1) = 0;

end
    
extra = repmat(exp(-1j*ko*z),size(f,1),1);
f = f.*extra;
end

function plot_mesh_function(X,Y,F,x_label,y_label,z_label,plot_title)

surf(X,Y,F);
view([0,90]);
axis tight, shading interp;
xlabel (x_label,'FontSize', 18);
ylabel (y_label,'FontSize', 18);
zlabel (z_label);
rotate3d on, 
% title(plot_title,'FontSize', 18);
colormap jet, colorbar;
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontName','Times','fontsize',18)


end