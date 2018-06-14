% Logan Prust - Aer E 322 - Group 4 - Lab 9 Plot

clear,clc;

% get data
loads = [3.14 6.14];
strain = [-11 -25 15 -1 -18 -18 7 -40.5 4 29; -20 -48 28 0 -34.5 -33 12 -63 5 57];

% fit lines to data
for i=1:10
    m(i) = (strain(2,i)-strain(1,i))/3;
    x0(i) = strain(2,i)-loads(2)*m(i);
end

% calculate lines
n=1000;
for j=1:n
    x(j) = j/100;
    for i=1:10
        line(i,j) = x(j)*m(i)+x0(i);
    end
end

% plot everything
figure(1);
hold on;
plot(loads,strain(:,1),'b*');
plot(loads,strain(:,2),'g*');
plot(loads,strain(:,3),'r*');
plot(loads,strain(:,4),'c*');
plot(loads,strain(:,5),'m*');
plot(loads,strain(:,6),'k*');
plot(loads,strain(:,7),'bo');
plot(loads,strain(:,8),'go');
plot(loads,strain(:,9),'ro');
plot(loads,strain(:,10),'co');
plot(x,line(1,:),'b');
plot(x,line(2,:),'g');
plot(x,line(3,:),'r');
plot(x,line(4,:),'c');
plot(x,line(5,:),'m');
plot(x,line(6,:),'k');
plot(x,line(7,:),'b');
plot(x,line(8,:),'g');
plot(x,line(9,:),'r');
plot(x,line(10,:),'c');
legend('Gauge 1','Gauge 2','Gauge 3','Gauge 4','Gauge 5','Gauge 6','Gauge 7','Gauge 8','Gauge 9','Gauge 10');
xlabel('Load (lbf)');
ylabel('Strain (microstrain)');
title('Linearity of Strain Readings with Load');