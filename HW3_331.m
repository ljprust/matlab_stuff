% Logan Prust - Aer E 331 - Homework 3

clear,clc;

% 1a

s=tf('s');
G=(s+4)/s/(s+2);
rlocus(G);

% 2b

% clear;
% s=tf('s');
% G=1/s/(s+5)/(s+10);
% rlocus(G);

% 3b

% s=tf('s');
% L=(s+2)/s^2/(s+4);
% rlocus(L);

% K=7.19;
% s=tf('s');
% G1=1/(s^2+K*(s+2)/(s+4));
% G2=K/(s^2*(s+4)/(s+2)+K);
% figure(1);
% step(G1*14.38/4);
% figure(2);
% step(G2);

% 4c

% s=tf('s');
% G=4.6/(s^2+0.76*s+4.55);
% rlocus(G);
% 
% G2=4.6*(s+4.919)/(s^2+0.76*s+4.55);
% rlocus(G2);
% step(G2/4.6);