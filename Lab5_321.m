clear,clc;

% twothetadeg=input('2 theta (degrees) = ');
n=input('n = ');

h=6.626E-34;
c=2.998E8;

d=2.015E-10;

% theta=twothetadeg/2*pi/180;
% 
% lambda=2/n*d*sin(theta);
% sigmalambda=1/n*d*cos(theta)*0.5*pi/180;
% 
% E=h*c/lambda/(1.6E-19)
% sigmaE=h*c/lambda^2*sigmalambda/(1.6E-19)

E=input('E (keV) = ');
E=E*1000*(1.6E-19);
lambda=h*c/E;
twotheta=asin(n*lambda/2/d)*180/pi*2;