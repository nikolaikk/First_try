clear
%%              Beam Propagaton Method
    
%%              Initialization

epsilon = 8.854187817e-12;
c = 299792458;

zo = 0;                     % micrometer
zend = 600;                % micrometer
z_mesh = 20;                 % micrometer
ratio = 1;
xo = -400;                  % micrometer
xend = -xo;                 % micrometer
x_mesh = z_mesh/ratio;      % micrometer
Lambda = 10;              % micrometer
wo = 25;                    % micrometer
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
%     range_z = [2,1*sqrt(2),1,sqrt(2)*0.5,0.5,0.25*sqrt(2),0.25];
    range_z = 1.6:-0.05:0.4;
%     range_z = [1.5,1,0.75,0.5,];
%     range_z = [1,0.75,0.5,0.375,0.25];


    % Video
    if record_status == 1
%         vid = VideoWriter(movie_name);
%         vid.FrameRate = 6;
%         vid.open()

        h = figure;
        set(h,'color','w');
    end
%     
    for i_z_mesh = range_z

        x_mesh = 1;      % micrometer
        
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
        
        
        %%
%         fig1 = figure;
%         set(fig1,'Position',[0,0,1000,1200])
%         [Z, X] = meshgrid(z,x);
%         subplot(2,1,1)
%         plot_mesh_function(Z,X,abs(real(Ez)),'z (\mum)','x (\mum)',...
%             'Real Part Electric Field',...
%             {'\Ree(E_{z}) Crank–Nicolson BPM',strcat('x-meshsize =',num2str(x_mesh),...
%             '\mum and z-meshsize =',num2str(i_z_mesh),' \mum')})
%         subplot(2,1,2)
%         plot_mesh_function(Z,X,abs(real(Ez_analytic)),'z (\mum)','x (\mum)',...
%             'Real Part Electric Field','\Ree(E_{z}) Analytic Gaussian Beam')
        %%
        
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
%             Ez_diff(xp_coord,zp_coord)=0.0001;
            [Z, X] = meshgrid(z,x);
            ddx =1200,ddy =1200;
            subplot(1,2,1)
            plot_mesh_function(imresize(Z,[ddx,ddy]),imresize(X,[ddx,ddy]),...
                (abs(real(imresize(Ez,[ddx,ddy])))),'z (\mum)','x (\mum)',...
                'Real Part Electric Field',{'\Ree(E_{z}) Crank–Nicolson BPM',...
                strcat('z\_mesh = ', num2str(i_z_mesh),'(\mum), x\_mesh = ',...
                num2str(x_mesh),'(\mum), ratio_{z\_mesh/x\_mesh} =',num2str(i_z_mesh/x_mesh))})
            subplot(1,2,2)
            plot_mesh_function(imresize(Z,[ddx,ddy]),imresize(X,[ddx,ddy]),...
                (abs(real(imresize(Ez_analytic,[ddx,ddy])))),'z (\mum)','x (\mum)',...
                'Real Part Electric Field','\Ree(E_{z}) Analytic Gaussian Beam')
            hold on
%             scatter(zp*1600,xp*800),hold off;
            colorbar, set(gcf,'Position',[0,0,1600,700]);
            frame = getframe(h);
            im = frame2im(frame);
            [imind,cm] = rgb2ind(im,256);
            if i_z_mesh == range_z(1) 
              imwrite(imind,cm,'test.gif','gif', 'Loopcount',inf); 
            else 
              imwrite(imind,cm,'test.gif','gif','WriteMode','append')
            end
            


            
            
%             vid.writeVideo(getframe(gca));
        end
    end
    
    % Record
    if record_status == 1
%         vid.close()
    end
    toc

    fig = figure,set(fig,'Position',[0,0,1500,500]);
    
    Linewidth= 3;
    subplot(1,3,1),hold on, plot(Number_Nz_Nx,err_l2,'-b','Linewidth',Linewidth)
    plot(Number_Nz_Nx,err_l2,'*k','Linewidth',Linewidth),hold off;
    title('L_{2} norm error','FontSize', 18)
    xlabel('Number of grid points','FontSize', 18),
    ylabel('Error','FontSize', 18);
    set(gca, 'XScale', 'log'),set(gca, 'YScale', 'log')
    a = get(gca,'XTickLabel');
    set(gca,'XTickLabel',a,'FontName','Times','fontsize',18)
   
    subplot(1,3,2), hold on, plot(Number_Nz_Nx,err_max,'-b','Linewidth',Linewidth)
    plot(Number_Nz_Nx,err_max,'*k','Linewidth',Linewidth),hold off;
    title('L_{\infty} norm error','FontSize', 18);
    xlabel('Number of grid points','FontSize', 18),
    ylabel('Error','FontSize', 18);
    set(gca, 'XScale', 'log'),set(gca, 'YScale', 'log')
    a = get(gca,'XTickLabel');
    set(gca,'XTickLabel',a,'FontName','Times','fontsize',18)
    
    subplot(1,3,3), hold on, plot(Number_Nz_Nx,err_point,'-b','Linewidth',Linewidth)
    plot(Number_Nz_Nx,err_point,'*k','Linewidth',Linewidth),hold off;
    title({'Error at a particular physical point ',strcat('( x=',num2str(xp*(xend-xo)-xend),...
        '\mum, z=', num2str(zp*zend),'\mum)')},'FontSize', 18)
    xlabel('Number of grid points','FontSize', 18),
    ylabel('Error','FontSize', 18);
    set(gca, 'XScale', 'log'),set(gca, 'YScale', 'log')
    a = get(gca,'XTickLabel');
    set(gca,'XTickLabel',a,'FontName','Times','fontsize',18)
    
    
    
    
%     % Show error and zrange    
%     fprintf(strcat('delta_z values  :   ',num2str(range_z),'\n'))
% %     fprintf(strcat('L2 normed error values  :   ',num2str(err_l2),'\n'))
%     fprintf(strcat('Number of grid points devided by max  :   ',...
%         num2str(Number_Nz_Nx/Number_Nz_Nx(1)),'   so increases x',...
%         num2str(round(Number_Nz_Nx(2)/Number_Nz_Nx(1))),'\n'))
%     fprintf(strcat('L2 normed error values devided by max  :   ',...
%         num2str(err_l2/err_l2(1)),'   error decreases x',...
%         num2str(round(err_l2(1)/err_l2(2))),'\n'))
%     fprintf(strcat('Lmax normed error values devided by max  :   ',...
%         num2str(err_max/err_max(1)),'   error decreases x',...
%         num2str(round(err_max(1)/err_max(2))),'\n'))
%     fprintf(strcat('Error at certain point devided by max  :   ',...
%         num2str(err_point/err_point(1)),'   error decreases x',...
%         num2str(round(err_point(1)/err_point(2))),'\n'))    
    
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

%%              Absorbtion region
region = 1/2;
k_grad = ((min(x)+max(x)*region-x(1:round(length(x)*region/2)))*region/2).^2;
k_abs = zeros(1,length(x));
k_abs(1:length(k_grad))= k_grad;
k_abs(end-length(k_grad)+1:end) = flip(k_grad);
% figure, plot(x,k_abs)
n = n + 0.0013j*k_abs;
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

E = eye(Nx);
const_to_simplify = (pinv(E-x_mesh*L/(4*1j*ko))*(E+x_mesh*L/(4*1j*ko)));
% const_to_simplify(find(const_to_simplify<1e-14))=0
for i=1:Nz-1

f(:,i+1)= const_to_simplify*f(:,i);
% f(:,i+1) = f(:,i+1)*cos(ko*z) 

f(1,i+1) = 0;
f(end,i+1) = 0;

end
    
periodic = repmat(exp(-1j*ko*z),size(f,1),1);
f = f.*periodic;
end

function plot_mesh_function(X,Y,F,x_label,y_label,z_label,plot_title)

surf(X,Y,F);
view([0,90]);
axis tight, shading interp;
xlabel (x_label,'FontSize', 18);
ylabel (y_label,'FontSize', 18);
zlabel (z_label);
rotate3d on, 
title(plot_title,'FontSize', 18);
colormap jet, colorbar;
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontName','Times','fontsize',18)


end