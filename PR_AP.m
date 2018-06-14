clear, clc
% Alexis Palliser %
% October 28, 2016 %
% AerE 351 %
% Section 1 %
% Project Rocket %

m_f = 2600; %total mass of stage 3 in kg
%m_0 = m_star + m_s + m_p; %Initial mass = payload mass + structure mass + propellant mass
%m_f = m_star + m_s; %Final mass = payload mass + structure mass
%m_p = m_0 - m_f; %Mass of propellant
GTO_palt = 72; %Geostationary transfer orbit altitude at perigee in km GTO_aalt = 35786; %Geostationary transfer orbbit altitude at apogee in km r_earth = 6378; %Radius of Earth in km
GTO_aalt = 35786; %Geostationary transfer orbbit altitude at apogee in km r_earth = 6378; %Radius of Earth in km
r_earth = 6378;
r_p = r_earth + GTO_palt; %Radius at perigee of the geostationary transfer orbit
r_a = r_earth + GTO_aalt; %Radius at apogee of the geostationary transfer orbit
a = (r_p + r_a) / 2; %Semi-major axis
ecc = (r_a - r_p) / (r_a + r_p); %Solving for eccentricity
e_1 = 0.12; %First stage structure ratio
Isp_1 = 300; %First stage engine's specific impulse in seconds
e_2 = 0.14; %Second stage structure ratio
Isp_2 = 455; %Second stage engine's specific impulse in seconds
Isp_3 = 240; %Apogee motor's specific impulse in seconds
e_3 = 0.15; %Structure ratio for apogee motor, tank, fuel, and associated structures
nu_GTO = 15; %True anomaly in degrees when joining GTO trajectory
g = 0.00981; %Earth's gravitational constant in km/s^2
%V_e = Isp*g; %Exhaust velocity
V_e_1 = Isp_1 * g; %Exhaust velocity at stage 1
V_e_2 = Isp_2 * g; %Exhaust velocity at stage 2
V_e_3 = Isp_3 * g; %Exhaust velocity at stage 3
%pie = m_star / m_0; %Payload ratio
%e = m_s / (m_s + m_p); %Structure ratio
%V_bo = (-V_e)*ln(e+(1-e)*pie); %Rocket equation
%when payload vanishes, pi goes to zero
%0.15 = m_s /(m_s + m_p)
mu = 3.986*10^5; %Defining the gravitational parameter
a_1 = r_p; %Semi-major axis 1
a_2 = r_a; %Semi-major axis 2
V_tp = (((2 * mu) / a_1) - ((2 * mu) / (a_1 + a_2)))^(1/2); %Tangential velocity at perigee
V_ta = (((2 * mu) / a_2) - ((2 * mu) / (a_1 + a_2)))^(1/2); %Tangential velocity at apogee
V_c1 = (mu / a_1)^(1/2);
V_c2 = (mu / a_2)^(1/2);
i = 28.5; %Inclination angle
delta_V3_1 = V_tp - V_c1;
delta_V3_2 = ((V_c2)^2 + (V_ta)^2 - (2 * V_c2 * V_ta * cosd(i)))^(1/2); %hange in velocity for stage 3 pie_3 = ((exp(delta_V3_total / -V_e_3)) - 0.15) / 0.85; %Payload ratio at stage 3
delta_V3_total = delta_V3_1 + delta_V3_2; %Change in velocity for stage 3
pie_3 = ((exp(delta_V3_total / -V_e_3)) - 0.15) / 0.85; %Payload ratio at stage 3
e_t = (a_2 - a_1)/(a_1 + a_2); %Tangential eccentricity
a_t = (a_1 + a_2)/2; %Tangential semi-major axis
p_t = (2 * a_1 * a_2)/(a_1 + a_2); %Tangential parameter (semi-lactus rectum) r_t = p_t / (1 + e_t * cosd(nu_GTO)); %Tangential radius
r_t = p_t / (1 + e_t * cosd(nu_GTO)); %Tangential radius
v_t = ((2 * mu / r_t) - (mu / a_t))^(1/2); %Tangential velocity
c_earth = 2 * pi * r_earth; %Circumference of Earth
sidereal_day = 86164; %The amount of seconds in 1 sidereal day
v_earth = (c_earth / sidereal_day) * cosd(i); %Velocity of Earth's rotation in km/s
v_adjust = 1.6; %Adjustment for velocity to accomodate drag and other variables
v_t = v_t - v_earth + v_adjust; %Tangential velocity corrected
%v_t is the velocity that stage 1 and 2 need to reach
%v_t = delta_V1_total + delta_V2_total
m_0_3 = 2600; %Initial mass or total mass at stage 3
m_star_3 = pie_3 * m_0_3; %Mass of payload at stage 3
m_s_3 = m_star_3 * 14.8849; %Mass of structure at stage 3
m_p_3 = m_star_3 * 84.3483; %Mass of propellant at stage 3
%m_01 = m_02 / (-0.170819 - (0.904013 * ((2600 / m_02)^1.51665))) %plot(m_02,m_02_m_01)
m_02 = (0:0.0001:1);
pie_2 = (0:0.0001:1);
for i = 1:10001
pie_1(i,:) = (0.0231535- 0.136364*(0.86*pie_2(i)+0.14).^1.51665)/((0.86*pie_2(i)+0.14).^1.51665);
m_0_3 = 2600; %Initial mass of stage 3
m_star_2 = m_0_3; %Payload mass at stage 2
m_0_2 = m_star_2 / pie_2(i); %Initial mass of stage 2
m_star_1 = m_0_2; %Payload mass at stage 1
m_0_1(i) = m_star_1 / pie_1(i); %Initial mass of stage 1
end
plot(pie_2,m_0_1);
axis([0 1 -10^6 10^6]);
%plot(pie_2, pie_1); %Graphing Payload ratios to determine the optimum values for stage 1 and 2
xlabel('pi_2'); %Labeling the x axis
ylabel('pi_1'); %Labeling the y axis
title('pi_2 vs. pi_1'); %Giving the graph a title
M = max(pie_1);
N = min(pie_1);
pie_1 = 0.8847; %Payload ratio at stage 1
pie_2 = 0.0960; %Payload ratio at stage 2

delta_V_1 = -V_e_1*log(e_1+(1-e_1)*pie_1); %Change in velocity for stage 1 delta_V_2 = -V_e_2*log(e_2+(1-e_2)*pie_2); %Change in velocity for stage 2 m_s_2 = ((1-pie_2)*m_star_2)/(pie_2 + (pie_2*((1-e_2)/e_2))); %Structure mass at stage 2
delta_V_2 = -V_e_2*log(e_2+(1-e_2)*pie_2); %Change in velocity for stage 2
m_s_2 = ((1-pie_2)*m_star_2)/(pie_2 + (pie_2*((1-e_2)/e_2))); %Structure mass at stage 2
m_p_2 = ((1 - e_2) * m_s_2) / e_2; %Propellant mass at stage 2
m_s_1 = ((1-pie_1)*m_star_1)/(pie_1 + (pie_1*((1-e_1)/e_1))); %Structure mass at stage 1
m_p_1 = ((1 - e_1) * m_s_1) / e_1; %Propellant mass at stage 1
%V_star = -V_e_1*ln(e_1+(1-e_1)*pie_1)-V_e_2*ln(e_2+(1-e_2)*pie_2)- V_e_3*ln(e_3+(1-e_3)*pie_3)