function Output = function_BMP_FFT(zo, zend, z_mesh, xo, xend, x_mesh, Lambda, wo, no, bla)
% % % zo = 0;
% % % zend = 800;    % micrometer
% % % z_mesh = 1;    % micrometer
% % % z = zo:z_mesh:zend-z_mesh;
% % % Nz = length(z);
% % % ratio = 1;
% % % 
% % % Lambda = 0.6;   % micrometer
% % % no = 1;
% % % 
% % % xo = -200; %mvim
% % % xend = -xo; %m
% % % x_mesh = z_mesh/ratio;
% % % x = xo:x_mesh:xend-x_mesh;
% % % Nx = length(x);
% % % wo = 10; % micrometer

    x = xo:x_mesh:xend-x_mesh;
    z = zo:z_mesh:zend-z_mesh;
    Nx = length(x);
    Nz = length(z);

    n = 1;
    k = 2*pi/Lambda*n;
    ko = 2*pi/Lambda*no;
    f = zeros(Nx);
    
    dkx = pi/(xo)*(length(x)-1)/length(x);
    kx = ( (-length(x)/2):(length(x)/2-1) )*dkx;
    N = 0;

    %************Input Field  Parameters****************%

    fz=exp(-0.5*(x/wo).^2);        % initial field distribution
    Output = zeros(Nx,Nz);
    Output(:,1) = abs(fz).^2;

    %%%%%%%%%%%%%%%%%%%%%%%%
    %************PROPAGATION *****************%
    %%%%%%%%%%%%%%%%%%%%%%%%

    Fz=ifftshift(ifft(ifftshift(fz)));
    Fz_d=Fz.*exp(1i*(kx.^2)./(2*ko*(no))*z_mesh/2);
    fz=fftshift(fft(fftshift(Fz_d)));

    for itr=1:Nz                                % was itr=1:Nz-1  !!!!!!!!!!

        fz_N=fz.*exp(-1i*N*ko*z_mesh);
        Fz=ifftshift(ifft(ifftshift(fz_N)));
        D=exp(1i*(kx.^2)./(2*ko*(no))*z_mesh);
        Fz_D=Fz.*D;
        fz=fftshift(fft(fftshift(Fz_D)));

        Output(:,itr)=abs(fz).^2;

    end

    if bla == true  
        [Z, X] = meshgrid(z,x);
        mesh(Z,X,(Output));
        view([0,90]);
        title('FFT BPM'),
        axis tight;
        shading interp;
        xlabel ('z (\mum)');
        ylabel ('x (\mum)');
        zlabel Intensity;
        colormap jet;
        rotate3d on
    end
    
% % %     [Z, X] = meshgrid(z,x);
% % %     figure
% % %     surf(Z,X,Output)
% % %     shading flat
% % %     xlabel ('z (\mum)');
% % %     ylabel ('x (\mum)');
% % %     zlabel Intensity;
% % %     rotate3d on
% % %     view([0,90])
% % % 
% % %     % imshow(imbinarize(Output,edge))
end