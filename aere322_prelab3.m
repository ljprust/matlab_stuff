% Aer E 322 - Lab 3 - Section 4 - Group 5

clear;

% get info about column
type = input('Type of column "c" or "r": ');
E = input('E: ');
L = input('L: ');
ends = input('1 pivot-pivot, 2 fixed-fixed, 3 fixed-free, 4 pivot-fixed: ');

% get relevant dimensions
if (type == 'c')
    R = input('R: ');
else
    B = input('B (larger): ');
    H = input('H (smaller): ');
end

% set effective length
if (ends == 1)
    Leff = L;
elseif (ends == 2)
    Leff = L/2;
elseif (ends == 3)
    Leff = 2*L;
else
    Leff = 0.7*L;
end

% calculate Pcr and slenderness ratio
if (type == 'c')
    Pcr = pi^3*E*R^4/4/Leff^2;
    slend = 2*Leff/R;
else
    Pcr = pi^2*E*B*H^3/12/Leff^2;
    slend = 2*sqrt(3)*Leff/H;
end

% output results
disp('Critical Load Pcr:');
disp(Pcr);
disp('Slenderness Ratio:');
disp(slend);