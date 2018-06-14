%Astro 342 Project Part 4
%12/4/2015
%Alicia Ekman
%***********************************************************************

clear, clc
close all

data=dlmread('summary1.txt');
% 
% 1 step #
% 2 Age years
% 3 mass msolar
% 4 luminosity ls
% 5 radius rs
% 6 surface temp K
% 7 Central temp K
% 8 Central density
% 9 central pressure
% 10 central electron degeneracy parameter

G=6.673*10.^-11; %Nm^2kg^-2
sig=5.67*10.^(-8); %W/m^2K^4
P=3.70*24*60*60; %seconds
Ab=.2;
em=.96;

step=data(:,1);
age=data(:,2); %years

mass=data(:,3); %mass_solar
mass_initial=data(1,3); %mass_solar
%mass_initial=0.912 %mass_solar

lumin=10.^(data(:,4)); %lum_solar
lumin_W=lumin*3.839*10.^26; %W
rad=10.^(data(:,5)); %rad_solar

rad_orbital_m=((G*mass*1.989*10.^30)./(4*pi.^2)*P.^2).^(1./3); %m
rad_orbital=rad_orbital_m./(6.955*10.^8); %rad_solar
rad_orbital_initial=(((G*mass_initial*1.989*10.^30)./(4*pi.^2)*P.^2).^(1./3))./(6.955*10.^8); %rad_solar
%rad_orbital_initial=9.7659 %rad_solar
rad_orbital_change=rad_orbital_initial*(mass./mass_initial)./(2*mass./mass_initial-1); %rad_solar

surf_temp=10.^(data(:,6)); %K
surf_temp_p=((1-Ab)./(em*sig)*(lumin_W)./(16*pi*(rad_orbital_m).^2)).^(1./4); %K
cen_temp=10.^(data(:,7));

% plot(age,mass)
% title('Mass over time')
% xlabel('Age (years)')
% ylabel('Mass (Solar Mass)')
% axis([0 16*10.^9 0.91 0.9125])

% plot(age,mass)
% title('Mass over time')
% xlabel('Age (years)')
% ylabel('Mass (Solar Mass)')
% axis([16*10.^9 17.5*10.^9 0.59 0.915])

% figure
% plot(age,lumin)
% title('Luminosity over time')
% xlabel('Age (years)')
% ylabel('Luminosity (Solar Luminosity)')
% axis([0 16*10.^9 0.45 2])

% figure
% plot(age,lumin)
% title('Luminosity over time')
% xlabel('Age (years)')
% ylabel('Luminosity (Solar Luminosity)')
% axis([16*10.^9 17.5*10.^9 0.9 6000])

% figure
% plot(age,rad)
% title('Radius over time')
% xlabel('Age (years)')
% ylabel('Radius (Solar Radius)')
% axis([0 16*10.^9 0.8 2])

% figure
% plot(age,rad)
% title('Radius over time')
% xlabel('Age (years)')
% ylabel('Radius (Solar Radius)')
% axis([16*10.^9 17.5*10.^9 2 250])

% figure
% plot(age,surf_temp)
% title('Surface Temperature over time')
% xlabel('Age (years)')
% ylabel('Surface Temp (K)')
% axis([0 16*10.^9 5250 5650])

% figure
% plot(age,surf_temp)
% title('Surface Temperature over time')
% xlabel('Age (years)')
% ylabel('Surface Temp (K)')
% axis([16*10.^9 17.5*10.^9 3000 5650])



figure
plot(age,rad_orbital,age,rad)
title('Orbital vs Stellar Radius over time')
xlabel('Stellar Age (years)')
ylabel('Orbital and Stellar Radius (Solar Radius)')





% figure
% plot(age,surf_temp_p)
% title('Planet Surface Temperature over time')
% xlabel('Stellar Age (years)')
% ylabel('Surface Temperature (K)')
% axis([0 16*10.^9 1020 1400])
% 
% figure
% plot(age,surf_temp_p)
% title('Planet Surface Temperature over time')
% xlabel('Stellar Age (years)')
% ylabel('Surface Temperature (K)')
% axis([16*10.^9 17.5*10.^9 1400 12000])