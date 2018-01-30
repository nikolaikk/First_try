from = 10
to = 800
err = ones(to-from,1)';
ii = 1;

for Nz=from:to-1
     
    zo = 0;
    z_length = 0.3; %m
%     Nz = 200; 
    delta_z = (z_length-zo)/Nz; % step size
    z=zo:delta_z:z_length;
    % z=z(1:100);

    Lambda = 500e-9;
    no=1;
    ko=2*pi/Lambda*no;
    n=1;
    k=2*pi/Lambda*n;

    xo=-0.001; %mvim
    xend=-xo; %m
    delta_x=(xend-xo)/(Nz-1);

    x=xo:delta_x:xend;
    wo=0.0001; %m
    f=zeros(Nz);
    f(:,1)=exp(-(x/wo).^2);


    %% Define L matrix
    L = zeros(Nz);
    for i=1:Nz
        L(i,i)=-2+(ko^2-k^2)*delta_x^2; % in free space^2 - in material^2
        end
    for i=2:(Nz)
        L(i,i-1)=1;
        L(i-1,i)=1;
        end
    L = (1/delta_x^2)*L;
    %%
    E = eye(Nz);


    %% Psi evolution

    const_to_simplify = ((E+delta_z*L/(4*1j*ko))*pinv(E-delta_z*L/(4*1j*ko)));
    % const_to_simplify(find(const_to_simplify<1e-14))=0
    for i=1:Nz-1
    % for i=1:2

        f(:,i+1)= const_to_simplify*f(:,i);

    end


    [Z, X] = meshgrid(delta_z*(1:Nz),delta_x*(1:Nz));
    I = conj(f).*f;

%     I(find(I<0.01)) = -0.2;
    I_analytic = function_analytic_BMP(Lambda,Nz,wo,z_length,1);
%     figure,mesh(Z,X,(I));
%     figure,mesh(Z,X,((I-I_analytic).^2));
%     axis vis3d;
%     shading interp;
%     xlabel ('z (m)');
%     ylabel ('x (m)');
%     zlabel I;
%     rotate3d on

    err(ii) = mean2((I-I_analytic).^2);

    ii=ii+1
end


figure, plot(from:to-1,err)

figure,mesh(Z,X,(I));
axis vis3d;
shading interp;
xlabel ('z (m)');
ylabel ('x (m)');
zlabel I;
rotate3d on

figure,mesh(Z,X,((I-I_analytic).^2));
axis normal;
shading interp;
xlabel ('z (m)');
ylabel ('x (m)');
zlabel I;
zlim([0,1]);
xlim([0,z_length]);
ylim([-xo,xo]);
rotate3d on
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
xlabel('Number of grid points')
ylabel('Error')




