% Logan Prust - Aer E 322 - Section 4 - Prelab 6

clear;

% define constants
L = 0.5; % m
E = 68.9E9; % Pa
g = 9.81; % m/s^2

% get input parameters
m = input('Mass of Weight (g): ');
m = m/1000; % kg
B = input('Base of Beam (mm): ');
B = B/1000; % m
H = input('Height of Beam (mm): ');
H = H/1000; % m

% calculate stuff
I = 1/12*B*H^3;
w = m*g;
y = -w/E/I*14/81*L^3;

% output results
fprintf('%20s %.15f %s \n','y =',y*1000,'mm');