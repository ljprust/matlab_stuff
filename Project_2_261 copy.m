% Logan Prust Aer E 261 Project 2 Main Program
clear,clc; % clear all variables as well as the screen
R=dlmread('reentryconstants.txt',',',[3,0,3,0]); % input specific gas constant
gamma=dlmread('reentryconstants.txt',',',[6,0,6,0]); % input specific heat ratio
TSL=dlmread('reentryconstants.txt',',',[9,0,9,0]); % input temperature at sea level
rohSL=dlmread('reentryconstants.txt',',',[12,0,12,0]); % input density at sea level
g0=dlmread('reentryconstants.txt',',',[18,0,18,0]); % input acceleration of gravity at sea level
L=dlmread('reentryconstants.txt',',',[21,0,21,0]); % input specific length of space shuttle
C1=dlmread('reentryconstants.txt',',',[24,0,24,0]); % input constant 1 for Sutherland's law
C2=dlmread('reentryconstants.txt',',',[27,0,27,0]); % input constant 2 for Sutherland's law
rE=dlmread('reentryconstants.txt',',',[30,0,30,0]); % input radius of Earth
lapserate=dlmread('reentryconstants.txt',',',[33,0,33,0]); % input lapse rate
hI=dlmread('reentryconstants.txt',',',[36,0,36,0]); % input altitude at bottom of isothermal layer
reentryinnput=dlmread('reentryinput.txt',' '); % input geometric altitude, velocity, and time matrix
hG=reentryinnput(:,1); % create geometric altitude matrix
v=reentryinnput(:,2); % create velocity matrix
for i=1:length(hG) % for-loop with length(hG) iterations
    h(i)=geopotentialaltitude(hG(i),rE); % create geopotential altitude matrix
end % end for-loop
TI0=Ttransition(TSL,lapserate,hI); % calculate temperature at bottom of isothermal layer
rohI0=rohtransition(rohSL,TI0,TSL,g0,lapserate,R); % calculate density at bottom of isothermal layer
fid=fopen('Project_2_Table.txt','w+'); % open output file
fprintf(fid,'%24s %20s %17s \n \n','Geometric Altitude (m)','Reynolds Number','Mach Number'); % write column headings to output file
for i=1:length(h) % for-loop for creating output table
    if h(i)>hI % if the shuttle is in the isothermal layer...
        T(i)=TI0; % the temperature is constant (and equal to the temperature at the bottom of the isothermal layer)
        roh(i)=rohisothermal(rohI0,g0,R,T(i),h(i),hI); % calculate density in isothermal layer
        a(i)=speedofsound(gamma,R,T(i)); % calculate speed of sound in isothermal layer
        mu(i)=dynamicviscosity(C1,T(i),C2); % calculate dynamic viscosity in isothermal layer
        Re(i)=reynolds(roh(i),v(i),L,mu(i)); % calculate Reynolds number in isothermal layer
        M(i)=mach(v(i),a(i)); % calculate Mach number in isothermal layer
    else % if the shuttle is in the gradient layer...
        T(i)=Tgradient(TSL,lapserate,h(i)); % calculate temperature in gradient layer
        roh(i)=rohgradient(rohSL,T(i),TSL,g0,lapserate,R); % calculate density in gradient layer
        a(i)=speedofsound(gamma,R,T(i)); % calculate speed of sound in gradient layer
        mu(i)=dynamicviscosity(C1,T(i),C2); % calculate dynamic viscosity in gradient layer
        Re(i)=reynolds(roh(i),v(i),L,mu(i)); % calculate Reynolds number in gradient layer
        M(i)=mach(v(i),a(i)); % calculate Mach number in gradient layer
    end % end if-else
    fprintf(fid,'%17.4f %24.0f %17.4f \n',hG(i),Re(i),M(i)); % write output table to output file
end % end for-loop for creating output file
fclose(fid); % close output file
figure(1) % create figure 1
plot(Re,hG) % plot Reynolds number vs. geometric altitude
title('Reynolds Number vs. Geometric Altitude') % put title on figure 1
xlabel('Reynolds Number') % label x-axis
ylabel('Geometric Altitude (m)') % label y-axis
figure(2) % create figure 2
plot(M,hG) % plot Mach number vs. geometric altitude
title('Mach Number vs. Geometric Altitude') % put title on figure 2
xlabel('Mach Number') % label x-axis
ylabel('Geometric Altitude (m)') % label y-axis