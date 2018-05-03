
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

xlabel ('x (\mu m)','FontSize', 18);
ylabel ('\kappa','FontSize', 18);

title('Absorption Index','FontSize', 18);

a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontName','Times','fontsize',18)



%% For absorption 

zo = 0;                     % micrometer
zend = 500;                % micrometer
z_mesh = 1;                 % micrometer
ratio = 1;
xo = -100;                  % micrometer
xend = -xo;                 % micrometer
x_mesh = z_mesh/ratio;      % micrometer
Lambda = 10;              % micrometer
wo = 15;                    % micrometer
no = 1;