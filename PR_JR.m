clear, clc

%AerE 351 Rocket Optimization Program
%Josh Roney, Created 10/26/16

re= 6378; %radius of Earth in km
rp= 72+re; %radius of perigee in km
ra= 35786+re; %radius of apogee in km
mu= 3.986e5; %gravitational parameter of Earth in km^2/s^3
eT= .734644; %eccentricity of transfer orbit
latitude= 28.5; %Geodetic latitude of KSC in degrees
agto= (rp+ra)/2; %Semi-major axis of transfer orbit
% V0= .34*cosd(latitude); %Initial speed when launching eastward
V0 = 2*pi/24/3600*re*cosd(latitude);
ge= .00981; %km/s^2
m03= 2600; %Total mass of the 3rd stage (motor + payload)

%Structure ratios of all stages of the rocket
e(1)= .12;
e(2)= .14;
e(3)= .15;

%Isp values for all stages
Isp(1)= 300;
Isp(2)= 455;
Isp(3)= 240;

for i= 1:3
    Ve(i)= Isp(i)*ge;
end

%Radius where vehicle enters GTO trajectory
r1= (24307*(1-eT^2))/(1+eT*cosd(15));
V(1)= sqrt((2*mu/r1)-(mu/agto)); %Velocity at GTO trajectory entrance
V(2)= sqrt((2*mu/ra)-(mu/agto)); %Velocity at apogee of GTO
V(3)= sqrt(mu/ra); %Velocity of final geostationary circular orbit

%DeltaV required by apogee motor to circularize and change plane
deltaV3= sqrt((V(3))^2+(V(2))^2-(2*V(3)*V(2)*cosd(latitude)));

%Derived from structure ratio, payload ratio, and Vbo equation
payload3= (exp(-deltaV3/Ve(3))-e(3))/(1-e(3));
m04= payload3*m03;
ms3= m04*e(3)*((1/payload3)-1);
mp3= m04*((1-e(3))*((1/payload3)-1));
m3= m04+ms3+mp3;

%Creating arrays to run payload ratio through to find best combination
payload1= [0.0001:.0001:.3];
payload2= [0.0001:.0001:.3];
%payload1=linspace(0,0.122582,3000);
%payload1=0.123*ones(1,length(3000));
%payload2=0.074*ones(1,length(3000));

%Running through every element of both arrays to find best match that has
%sufficient delta V
for i= 1:length(payload1)
    deltaV1(i)= -Ve(1)*log(e(1)+((1-e(1))*payload1(i)));
    for j= 1:length(payload2)
        deltaV2(j)= -Ve(2)*log(e(2)+((1-e(2))*payload2(j)));
        if (deltaV1(i)+deltaV2(j)) >= (V(1)+1.6-V0)
            m02(j)= m03/payload2(j);
            m01(j)= m02(j)/payload1(i);
            payloadplot1(j)=m02(j)/m01(j);
            payloadplot2(j)=m03/m02(j);
        end
    end
end

%Algorithm to find smallest Stage 1 mass
m01min= m01(1);
for i= 1:length(m01)
    if m01min > m01(i)
        m01min= m01(i);
        payload1min= payloadplot1(i);
    end
end

%Calculations after finding minimal Stage 1 mass
m02min= payload1min*m01min;
payload2min= m03/m02min;
ms2= e(2)*m03*((1/payload2min)-1);
mp2= m03*((1-e(2))*((1/payload2min)-1));
ms1= e(1)*m02min*((1/payload1min)-1);
mp1= m02min*((1-e(1))*((1/payload1min)-1));

%Delta V calculations from previously, with minimal mass payload ratio
%found from optimization algorithm.
deltaV1actual= -Ve(1)*log(e(1)+((1-e(1))*payload1min));
deltaV2actual= -Ve(2)*log(e(2)+((1-e(2))*payload2min));

%Plotting Stage 1 mass vs. Payload ratio
plot(payloadplot1,m01,payload1min,m01min,'o')
ylim([3e5 8e5])
ylabel('Mass of 1st stage (kg)')
xlabel('\pi')

disp('Stage 1 total mass (kg)')
disp(m01min)
disp('Stage 1 structure mass (kg)')
disp(ms1)
disp('Stage 1 propellant mass (kg)')
disp(mp1)
disp('Stage 2 total mass (kg)')
disp(m02min)
disp('Stage 2 structure mass (kg)')
disp(ms2)
disp('Stage 2 propellant mass (kg)')
disp(mp2)
disp('Stage 3 total mass (kg)')
disp(m03)
disp('Stage 3 structure mass (kg)')
disp(ms3)
disp('Stage 3 propellant mass (kg)')
disp(mp3)
disp('Payload total mass (kg)')
disp(m04)
disp('Payload Ratio 1')
disp(payload1min)
disp('Payload Ratio 2')
disp(payload2min)
disp('Payload Ratio 3')
disp(payload3)
disp('Delta V 1')
disp(deltaV1actual)
disp('Delta V 2')
disp(deltaV2actual)
disp('Delta V 3')
disp(deltaV3)
