zo = 0;
zend = 3000;    % micrometer
z_mesh = 0.3;    % micrometer
z = zo:z_mesh:zend-z_mesh;
Nz = length(z);

Lambda = 0.6;   % micrometer
no = 1;
ko = 2*pi/Lambda*no;
n = 1;
k = 2*pi/Lambda*n;

xo = -100; %mvim
xend = -xo; %m
x_mesh = z_mesh;
x = xo:x_mesh:xend-x_mesh;
Nx = length(x);
wo = 10; % micrometer
% f = zeros(Nx);

dkx = pi/(xo)*(length(x)-1)/length(x);
kx = ( (-length(x)/2):(length(x)/2-1) )*dkx;


%************Input Field  Parameters****************%

fz=exp(-0.5*(x/wo).^2);        % initial field distribution
figure, plot(real(fz))
Output = zeros(Nx,Nz);
Output(:,1) = abs(fz).^2;

%%%%%%%%%%%%%%%%%%%%%%%%
%************PROPAGATION *****************%
%%%%%%%%%%%%%%%%%%%%%%%%

Fz=ifftshift(ifft(ifftshift(fz)));
figure, plot(real(Fz))
Fz_d=Fz.*exp(1i*(kx.^2)./(2*ko*(no))*z_mesh/2);
figure, plot(real(Fz))
fz=fftshift(fft(fftshift(Fz_d)));
figure, plot(fz)

for itr=1:Nz-1
    
    fz_N=fz.*exp(-1i*N*ko*z_mesh);
    Fz=ifftshift(ifft(ifftshift(fz_N)));
    D=exp(1i*(kx.^2)./(2*ko*(no))*z_mesh);
    Fz_D=Fz.*D;
    fz=fftshift(fft(fftshift(Fz_D)));
    
    Output(:,itr)=abs(fz).^2;
    
end

[Z, X] = meshgrid(z,x);
surf(Z,X,Output)
shading flat
xlabel ('z (\mum)');
ylabel ('x (\mum)');
zlabel Intensity;
rotate3d on