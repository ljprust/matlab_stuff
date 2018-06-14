% Aer E 322 - Lab 2 - Section 4 - Group 5

clear,clc;

% initialize vectors
n=1000;
radius = linspace(1,5,n);
sigmar = zeros(1,n);
sigmat = zeros(1,n);
tau = zeros(1,n);

% calculate stresses
for i = 1:n
    sigmar(i) = 0.5*(1 - 1/radius(i)^2 - (1 - 4/radius(i)^2 + 3/radius(i)^4)*cos(2*0));
    sigmat(i) = 0.5*(1 + 1/radius(i)^2 + (1 + 3/radius(i)^4)*cos(2*0));
    tau(i) = 0.5*(1 + 2/radius(i)^2 - 3/radius(i)^4)*sin(2*0);
end

% generate plots
figure(1);
plot(radius,sigmar);
title('Radial Stress');
xlabel('Radius in units of a');
ylabel('Stress in units of sigma0');

figure(2);
plot(radius,sigmat);
title('Azimuthal Stress');
xlabel('Radius in units of a');
ylabel('Stress in units of sigma0');

figure(3);
plot(radius,tau);
title('Shear Stress');
xlabel('Radius in units of a');
ylabel('Stress in units of sigma0');

clear;

% generate grid
x = -10:0.1:10;
y = -12:0.1:12;

% set hold size and applied stress to arbitrary values
a = 2;
sigma0 = 1;

% do calculations
for i = 1:size(x,2);
    for j = 1:size(y,2);
        % find polar coordinates
        r = sqrt(x(i)^2 + y(j)^2);
        theta = atan2(y(j),x(i));
        
        % set stresses to zero inside hole
        if (r < a)
            sigmar(i,j) = 0;
            sigmat(i,j) = 0;
            tau(i,j) = 0;
        else
            % calculate stresses
            sigmar(i,j) = sigma0/2*(1-a^2/r^2-(1-4*a^2/r^2+3*a^4/r^4)*cos(2*theta));
            sigmat(i,j) = sigma0/2*(1+a^2/r^2+(1+3*a^4/r^4)*cos(2*theta));
            tau(i,j) = sigma0/2*(1+2*a^2/r^2-3*a^4/r^4)*sin(2*theta);
        end
    end
end

% calculate diameter of Mohr's circle
deltasigma = 2*sqrt(((sigmar-sigmat)/2).^2 + tau.^2);

% generate plots
figure(4);
imagesc([-10 10],[-12 12],sigmar);
title('Radial Stress');
xlabel('x (arbitrary units)');
ylabel('y (arbitrary units)');

figure(5);
imagesc([-10 10],[-12 12],sigmat);
title('Azimuthal Stress');
xlabel('x (arbitrary units)');
ylabel('y (arbitrary units)');

figure(6);
imagesc([-10 10],[-12 12],tau);
title('Shear Stress');
xlabel('x (arbitrary units)');
ylabel('y (arbitrary units)');

figure(7);
imagesc([-10 10],[-12 12],deltasigma);
colormap([jet;jet;jet]);
title('Diameter of Mohrs Circle');
xlabel('x (arbitrary units)');
ylabel('y (arbitrary units)');