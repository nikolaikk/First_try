
n=1;
x=linspace(-400,399,800);
rigion = 1/2;
k_grad = ((min(x)+max(x)*rigion-x(1:round(length(x)*rigion/2)))*rigion/2).^2;

k_abs = zeros(1,length(x));
k_abs(1:length(k_grad))= k_grad;

k_max = 0.0013;
k_abs(end-length(k_grad)+1:end) = flip(k_grad);
% figure, plot(x,k_abs)
n = n + k_max*j*k_abs;
% k = 2*pi/Lambda*n;

fig1 = figure;
set(fig1,'Position',[0,0,1000,600])
plot(x,k_abs, 'linewidth',3)

xlabel ('x (\mum)','FontSize', 18);
ylabel ('\kappa(x)','FontSize', 18);

% title('Absorption Index','FontSize', 18);

a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontName','Times','fontsize',28)



%% For absorption 


epsilon = 8.854187817e-12;
c = 299792458;

zo = 0;                     % micrometer
zend = 800;                % micrometer
z_mesh = 1;                 % micrometer
ratio = 1;
xo = -400;                  % micrometer
xend = -xo;                 % micrometer
x_mesh = z_mesh/ratio;      % micrometer
Lambda = 20;              % micrometer
wo = 50;                    % micrometer
no = 1;
