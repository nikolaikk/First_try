%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Written for Course :- Computational Electromagnetics, Fall 2011
%                       Department of Electrical Engineering
%                       Indian Institute of Technology Madras
%                       Chennai - 600036, India
%
% Authors            :- Sathya Swaroop Ganta, B.Tech., M.Tech. Electrical Engg.
%                       Kayatri, M.S. Engineering Design
%                       Pankaj, M.S. Electrical Engg.
%                       Sumanthra Chaudhuri, M.S. Electrical Engg.
%                       Projesh Basu, M.S. Electrical Engg.
%                       Nikhil Kumar CS, M.S. Electrical Engg.
%
% Instructor :- Ananth Krishnan
%               Assistant Professor
%               Department of Electrical Engineering
%               Indian Institute of Technology Madras
%
% Any correspondance regarding this program may be addressed to
% Prof. Ananth Krishnan at 'computational.em.at.iit.madras@gmail.com'
%
% Copyright/Licensing :- For educational and research purposes only. No
% part of this program may be used for any financial benefit of any kind
% without the consent of the instructor of the course.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% "1D FDTD solution for Mur's Absorbing Boundary Condition (ABC)"
% 
% Objective of the program is to solve for the Maxwell's equation for an
% x-directed z-polarized TEM wave containg the y-directed magnetic
% field Hy and z-directed electric field Ez. The fields are updated at every 
% timestep, in a space, where all physical parameters of free space are
% not normalized to 1 (given real and known values), using standard update 
% equations obtained from the difference form of Maxwell's curl equations 
% with very low electric and magnetic conductivities of 4 x10^(-4) units 
% incorporated in them. The field points are defined in a grid described 
% by Yee's algorithm. The H fields are defined at every half coordinate of 
% a spacestep and  E field is defined at every coordinate point. Also, 
% the time update is done using Leapfrog time-stepping. Here, the H fields 
% are updated every half time-step and E fields are updated every full 
% time-step. This is shown by two alternating vector updates spanning only 
% a part of spatial grid where the wave, starting from source, has reached  
% at that particular time instant avoiding field updates at all points in 
% grid which is unnecessary at that time instant. These spatial updates are 
% inside the main for-loop for time update, spanning the entire time grid. 
% Also here, the space-step length is taken as 1 micron instead of 1 unit 
% in unitless domain assumed in previous programs. Also, here, the vectors
% used as multiplication factors for update equations are initialized
% before the loop starts to avoid repeated calculation of the same in every
% loop iteration, a minor attempt at optimization. The boundary condition
% here is Mur's Absorbing Boundary Condition (ABC) where the fields at the 
% grid points have electric field values formulated using Engquist Majda one
% way wave equations [1] where the boundaries give a sense of absorbing the 
% total field incident on them and reflecting none back to the domain.
%
% A source of electric field is defined at the center of the spatial domain 
% which is a hard source, in that it does not change its value due to 
% interference from external fields i.e in other words, the source is a 
% perfect electric conductor. The form of the source can be %     E(:,i) = (1/R(z(i),zR)-1j*Lambda/(pi*width(z(i),zR,wo)))...
%         *exp(-1j*k.*(z(i)+x.^2./(2*R(z(i),zR)*pi*width(z(i),zR,wo)/(pi*width(z(i),zR,wo)-1j*Lambda*R(z(i),zR)))));
varied using 
% the variables sine, gaussian and impulse. The source is available in four 
% standard forms- Unit-time step, Impulse, Gausian and Sinusoidal forms. 
% For the sinusoidal source, the user is allowed to give a frequency of
% his/her choice in Hz. The plot of Ez field v/s spacesteps is shown at 
% every time step. The simulation can be ended by closing this plot 
% window or by waiting till all the time step updates are completed.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Clearing variables in memory and Matlab command screen
clear all;
clc;

% Grid Dimension in x (xdim) direction
xdim=200;

%Total no of time steps
time_tot=350;

%Position of the source (center of the domain)
xsource=100;

%Courant stability factor
S=1;

% Parameters of free space (permittivity and permeability and speed of
% light) are all not 1 and are given real values
epsilon0=(1/(36*pi))*1e-9;
mu0=4*pi*1e-7;
c=3e+8;

% Spatial grid step length (spatial grid step = 1 micron)
delta=1e-6;
% Temporal gris step obtained using Courant condition
deltat=S*delta/c;

% Initialization of field vectors
Ez=zeros(1,xdim);
Hy=zeros(1,xdim);

% Initialization of permittivity and permeability vectors
epsilon=epsilon0*ones(1,xdim);
mu=mu0*ones(1,xdim);

% Initializing electric and magnetic field vectors
sigma=4e-4*ones(1,xdim);
sigma_star=4e-4*ones(1,xdim);

% Choice of nature of source
gaussian=0;
sine=1;
% The user can give a frequency of his choice for sinusoidal (if sine=1 above) waves in Hz 
frequency=0.75e+13;
impulse=0;
%Choose any one as 1 and rest as 0. Default (when all are 0): Unit time step

%Multiplication factor vectors for H vector update to avoid being calculated many times 
%in the time update loop so as to increase computation speed
A=((mu-0.5*deltat*sigma_star)./(mu+0.5*deltat*sigma_star)); 
B=(deltat/delta)./(mu+0.5*deltat*sigma_star);
                          
%Multiplication factor vectors for E vector update to avoid being calculated many times 
%in the time update loop so as to increase computation speed                          
C=((epsilon-0.5*deltat*sigma)./(epsilon+0.5*deltat*sigma)); 
D=(deltat/delta)./(epsilon+0.5*deltat*sigma);                     
                   

% Update loop begins
for n=1:1:time_tot
    
    % if source is impulse or unit-time step 
    if gaussian==0 && sine==0 && n==1
        Ez(xsource)=1;
    end
    
    % Setting time dependent boundaries to update only relevant parts of the 
    % vector where the wave has reached to avoid unnecessary updates.
    if n<xsource-2
        n1=xsource-n-1;
    else
        n1=1;
    end
    if n<xdim-1-xsource
        n2=xsource+n;
    else
        n2=xdim-1;
    end
    
    % Vector Update instead of for loop for Hy field incorporating magnetic
    % conductivity
    Hy(n1:n2)=A(n1:n2).*Hy(n1:n2)+B(n1:n2).*(Ez(n1+1:n2+1)-Ez(n1:n2));
    
    % Vector Update instead of for loop for Ez field incorporating magnetic
    % conductivity
    Ez(n1+1:n2+1)=C(n1+1:n2+1).*Ez(n1+1:n2+1)+D(n1+1:n2+1).*(Hy(n1+1:n2+1)-Hy(n1:n2));
    
    % Source conditions
    if impulse==0
        % If unit-time step
        if gaussian==0 && sine==0
            Ez(xsource)=1;
        end
        %if sine
        if sine==1
            tstart=1;
            N_lambda=c/(frequency*delta);
            Ez(xsource)=sin(((2*pi*(c/(delta*N_lambda))*(n-tstart)*deltat)));
        end
        %if gaussian
        if gaussian==1
            if n<=42
                Ez(xsource)=(10-15*cos(n*pi/20)+6*cos(2*n*pi/20)-cos(3*n*pi/20))/32;
            end
        end
    else
        %if impulse
        Ez(xsource)=0;
    end
    
    % Mur's Absorbing boundary condition
    if n>1
        Ez(xdim)=dum2+((S-1)/(S+1))*(Ez(xdim-1)-Ez(xdim));
        Ez(2)=dum1+((S-1)/(S+1))*(Ez(2)-Ez(1));
        Ez(1)=Ez(2);
    end
    dum1=Ez(3);
    dum2=Ez(xdim-1);
        
    
    plot((1:1:xdim)*delta,Ez,'color','k');
    titlestring=['FDTD Mur Absorbing Boundary Condition at time = ',num2str(round(n*deltat/10e-15)),' fs'];
    title(titlestring,'color','k');
    xlabel('x in m');
    ylabel('Ez in V/m');
    axis([0 xdim*delta -3 3]);
    getframe;
end
% [1] "Computational Electrodynamics - The Finite Difference Time Domain
%      Method" - Allen Taflove, Susan.c.Hagness, Third Edition, Artech House
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% END OF PROGRAM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%