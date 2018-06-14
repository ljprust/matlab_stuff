clear,clc;
t=linspace(0,10,1000);
for i=1:1000
    vy(i)=2*exp(-t(i))+1-exp(-t(i));
end
plot(t,vy);
title('Velocity vs. Time');
xlabel('Time');
ylabel('Velocity');