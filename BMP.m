clear
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       BEAM PROPAGATION METHOD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%                      Initialization
epsilon = 8.854187817e-12;
c = 299792458;

zo = 0;                     % micrometer
zend = 600;                 % micrometer
z_mesh = 1;                 % micrometer
ratio = 1;
xo = -400;                  % micrometer
xend = -xo;                 % micrometer
x_mesh = z_mesh/ratio;      % micrometer
Lambda = 10;                % micrometer
wo = 25;                    % micrometer
no = 1;

% Creates x and x one dimensional grid
x = xo:x_mesh:xend-x_mesh;
z = zo:z_mesh:zend-z_mesh;


%%                      Start calculation
% part_of_program 0 runs convergence, 1 plots the Gaussian Beam
part_of_program = 0;
part_of_program = 1;

% record_status 1 saves covergence as gif, 0 does not
record_status = 0;
record_status = 1;


%%                      Convergence
 
if part_of_program  == 0
    
    step = 0.1;
    from = 4;
    to = 0.5;
%     range_z = from:-step:to;
%     range_z = [2,1*sqrt(2),1,sqrt(2)*0.5,0.5,0.25*sqrt(2),0.25];
%     range_z = 1.6:-0.05:0.4;
    range_z = [1.5,1,0.75,0.5,];
    err_l2 = zeros([1,length(range_z)]);
    err_max = err_l2;
    err_point = err_l2;
    Number_Nz_Nx = err_l2;
    ii = 1;


    % saves as gif (1/2)
    if record_status == 1
        h = figure;
        set(h,'color','w');
    end

    % z meshsize and x meshize are changed, BMP applied and error
    % ... shown
    
    for i_z_mesh = range_z
        
        % creates x and x one dimensional grid
        x_mesh = i_z_mesh;                  % micrometer
        x = xo:x_mesh:xend-x_mesh;
        z = zo:i_z_mesh:zend-i_z_mesh;
        
        % derives numerical nad analytical Ez
        Ez = solve_BPM(z, x, Lambda, wo, no);
        Ez_analytic = function_analytic_BMP(z, x, Lambda, wo, no);
        
        % finds difference between numerical and analytical Ez
        Ez_diff = abs(abs(real(Ez))-abs(real(Ez_analytic)));

        % calculates L2 and Lmax errors
        Number_Nz_Nx(ii) = length(z)*length(x);
        err_l2(ii) = 1/Number_Nz_Nx(ii) * sqrt(sum(sum(Ez_diff.^2)));
        err_max(ii) = max(max(Ez_diff));
        
        %%
        
        % error at particular point at physical coordinate (zp,xp)
        zp = 0.5;    % from 0 to 1
        xp = 0.7;     % from 0 to 1
        point = [zp,xp];
       
        % transforms physical point (zp,xp) into numerical point for 
        % ... any grid size
        xp_coord = round(length(x)*xp)+1;
        zp_coord = round(length(z)*zp)+1;
        err_point(ii) = abs(real(Ez_analytic(xp_coord,zp_coord)-Ez(xp_coord,zp_coord)));
        
        % saves as gif (2/2)
        if record_status == 1
            [Z, X] = meshgrid(z,x);
            
            % resize Z, X and Ez, this is neccessary for gif saving
            ddx =1200,ddy =1200;
            Z = imresize(Z,[ddx,ddy]); Ez = imresize(Ez,[ddx,ddy]);
            X =imresize(X,[ddx,ddy]); Ez_analytic = imresize(Ez_analytic,[ddx,ddy]);
            
            % plot Ez and Ez_analytic
            subplot(1,2,1)
            plot_mesh_function(Z,X,abs(real(Ez)),'z (\mum)','x (\mum)',...
                'Real Part Electric Field',{'\Ree(E_{z}) Crank–Nicolson BPM',...
                strcat('z\_mesh = ', num2str(i_z_mesh),'(\mum), x\_mesh = ',...
                num2str(x_mesh),'(\mum), ratio_{z\_mesh/x\_mesh} =',num2str(i_z_mesh/x_mesh))})
            
            subplot(1,2,2)
            plot_mesh_function(Z,X,abs(real(Ez_analytic)),'z (\mum)','x (\mum)',...
                'Real Part Electric Field','\Ree(E_{z}) Analytic Gaussian Beam')
            colorbar, set(gcf,'Position',[0,0,1600,700]);
            
            % save each frame
            frame = getframe(h);
            im = frame2im(frame);
            [imind,cm] = rgb2ind(im,256);
            if i_z_mesh == range_z(1) 
              imwrite(imind,cm,'test2.gif','gif', 'Loopcount',inf); 
            else 
              imwrite(imind,cm,'test2.gif','gif','WriteMode','append')
            end
        end
        
        % shows the current proccesing i_z_mesh out of range_z
        clc;
        disp(strcat(num2str(ii),' /  ', num2str(length(range_z))));
        ii=ii+1;
    end

    % shows error convergence as figure and also in command window
    data = show_and_plot_error(Number_Nz_Nx,err_l2,err_max,err_point,x,z,range_z,point);    
   
%%                        Comparison of Ez and Ez_analytic

else
    % uses solve_BPM and function_analytic_BMP functions
    [Ez, L, const_to_simplify] = solve_BPM(z, x, Lambda, wo, no);
    Ez_analytic = function_analytic_BMP(z, x, Lambda, wo, no);
    
    % solves corresponding intensity
    I = epsilon*no*c*0.5*abs(Ez).^2;
    I_analytic = epsilon*no*c*0.5*abs(Ez_analytic).^2;

    [Z, X] = meshgrid(z,x);

    % creates structure for more convenient and fast plot
    plot_data.zx = {(abs(real(Ez))),I,abs(real(Ez_analytic)),I_analytic};
    plot_data.axis3_titles = {'Real Part Electric Field','Intensity',...
        'Real Part Electric Field','Intensity'};
    plot_data.titles = {'\Ree(E_{z}) Crank–Nicolson BPM','Intensity Crank–Nicolson BPM',...
        '\Ree(E_{z}) Analytic Gaussian Beam','Intensity Analytic Gaussian Beam'};
    
    % figure for plot
    figG = figure;
    set(figG,'Position',[0,0,2000,1200])
    
    % fills figure with subplots
    for i=1:length(plot_data.titles)
        subplot(2,2,i)
        plot_mesh_function(Z,X,plot_data.zx{i},'z (\mum)','x (\mum)',...
        plot_data.axis3_titles{i},plot_data.titles{i})
    end
end

%% ==================================================================== %%
%                              FUNCTIONS
%  ====================================================================  %

%%                       Main Numeric BPM funciton 
function [Ez,L,const_to_simplify] = solve_BPM(z, x, Lambda, wo, no)

    x_mesh = x(2)-x(1);

    Nz = length(z);
    Nx = length(x);

    ko = 2*pi/Lambda*no;
    n = 1;
    k = 2*pi/Lambda*n*ones(1,length(x));
    f = zeros(Nx,Nz);
    f(:,1) = exp(-(x/wo).^2);

    f(1,1) = 0;
    f(end,1) = 0;

    % creates absorbtion region, imaginary refractive index in specified
    % ... region is growing quadraticaly, outside the regio is equal to 0
    region = 1/2
    k_grad = ((min(x)+max(x)*region-x(1:round(length(x)*region/2)))*region/2).^2;
    k_abs = zeros(1,length(x));
    k_abs(1:length(k_grad))= k_grad;
    k_abs(end-length(k_grad)+1:end) = flip(k_grad);
    % figure, plot(x,k_abs)
    n = n + 0.0013j*k_abs;
    % k = 2*pi/Lambda*n;
    
    % Thomas algorithm, or tridiagonal matrix
    L = zeros(Nx);
    for i=1:Nx
        L(i,i)=(1/x_mesh^2)*(-2+(ko^2-k(i)^2)*x_mesh^2); % in free space^2 - in material^2
    end
    for i=2:(Nx)
        L(i,i-1)=(1/x_mesh^2)*1;
        L(i-1,i)=(1/x_mesh^2)*1;
    end

    Identity = eye(Nx);
    const_to_simplify = (Identity+x_mesh* L/(4*1j*ko))/(Identity-x_mesh*L/(4*1j*ko));
    
    % calculates f column at each z point
    for i=1:Nz-1

        f(:,i+1)= const_to_simplify*f(:,i);

    % no boundary conditions (reflection is expected)
    f(1,i+1) = 0;
    f(end,i+1) = 0;

    end

    % addes periodicity
    periodic = repmat(exp(-1j*ko*z),size(f,1),1);
    Ez = f.*periodic;
end

%%                       Main Analytic BPM funciton
function Ez = function_analytic_BMP(z,x, Lambda, wo, no)

    E = 1;
    k = 2* pi/Lambda*no;
    zR = pi*wo^2/Lambda;
    Nx = length(x);
    Nz = length(z);
    Ez = zeros([Nx,Nz]);
    
    % solves analytically beam propagation at each z point
    for i=1:Nz
        w = wo*sqrt(1+(z(i)/zR).^2);
        r = z(i)+(zR^2/z(i));
        gouy = atan(z(i)/zR);

        % Amplitude Normalized
        Ez(:,i) = E * (wo./w).^0.5.*exp(-(x.^2)./w^2).*exp(1j*(k*z(i)+k*(x.^2)./(2*r)-gouy));

    end
end

%%                       Plot Function for Ez
function plot_mesh_function(X,Y,E,x_label,y_label,z_label,plot_title)

    surf(X,Y,E);
    view([0,90]);
    axis tight, shading interp;
    xlabel (x_label,'FontSize', 18);
    ylabel (y_label,'FontSize', 18);
    zlabel (z_label);
    rotate3d on, 
    title(plot_title,'FontSize', 18);
    colormap jet; colorbar;
    a = get(gca,'XTickLabel');
    set(gca,'XTickLabel',a,'FontName','Times','fontsize',18)

end

%%                        Plot Function for Error 
function error_data = show_and_plot_error(Number_Nz_Nx,err_l2,err_max,...
    err_point,x,z,range_z,point)

    zp = point(1);
    xp = point(2);
    zend = z(end);
    
    xo= x(1);
    xend = -xo;

    % creates structure for more convenient and fast plot
    error_data.error = {err_l2,err_max,err_point};
    error_data.titles = {'L_{2} norm error','L_{\infty} norm error',...
        {'Error at a particular physical point ',...
        strcat('( x=',num2str(xp*(xend-xo)-xend),'\mum, z=', num2str(zp*zend),'\mum)')}};
    
    % figure for plot
    figG = figure;
    set(figG,'Position',[0,0,1500,500])
    
    % fills figure with subplots
    for i=1:length(error_data.titles)
        Linewidth= 3;
        subplot(1,3,i)
        hold on, plot(Number_Nz_Nx,error_data.error{i},'-b','Linewidth',Linewidth)
        plot(Number_Nz_Nx,error_data.error{i},'*k','Linewidth',Linewidth),hold off;
        title(error_data.titles{i},'FontSize', 18)
        xlabel('Number of grid points','FontSize', 18),
        ylabel('Error','FontSize', 18);
        set(gca, 'XScale', 'log'),set(gca, 'YScale', 'log')
        a = get(gca,'XTickLabel');
        set(gca,'XTickLabel',a,'FontName','Times','fontsize',18)
    end
    
    % Show error and zrange    
    fprintf(strcat('delta_z values  :   ',num2str(range_z),'\n'))
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
end