% Logan Prust - Aer E 322 - Section 4 - Prelab 7

clear;

% get type of section (C-channel or open pipe)
type = input('Type of Section (c or p): ');

if (type == 'c')
    
    % get specs of section
    h = input('Height: ');
    b = input('Base: ');
    t = input('Thickness: ');
    
    % calculate stuff
    I = 1/12*t*h^3 + 2*(h/2)^2*b*t;
    e = h^2*b^2*t/4/I;
    
else
    
    % get specs of section
    D = input('Diameter: ');
    r = D/2; % radius
    theta0 = input('Theta0: ');
    theta0 = theta0*pi/180; % radians
    t = input('Thickness: ');
    
    % calculate stuff
    I = t*r^3*(pi - theta0 - 0.5*sin(2*theta0));
    e = 2*r*(cos(theta0)*(2*pi-2*theta0)+2*sin(theta0))/(2*pi-2*theta0+sin(2*theta0));
    
end

% output results
fprintf('%20s %.15f \n','I =',I);
fprintf('%20s %.15f \n','e =',e);