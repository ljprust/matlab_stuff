% Logan Prust - Math 407 - HW 2 - Problem 2.1.2

A=[1,3;1,0]; % define matrix

vals=eig(A); % get eigenvalues
xvals=real(vals); % put real values on the x-axis
yvals=imag(vals); % put imaginary values on the y-axis

rowRadii = sum(abs(A),2) - abs(diag(A)); % sum across rows, neglecting the diagonal elements

figure;
hold on;
plot(xvals,yvals,'k*'); % plot eigenvalues
plot_circle(gca,real(A(1,1)),imag(A(1,1)),rowRadii(1)); % plot disks
plot_circle(gca,real(A(2,2)),imag(A(2,2)),rowRadii(2));
hold off;
axis equal; % set aspect ratio
draw_axes(gca); % draw axes on the plot
title('Row Disks for Problem 2.1.1');
xlabel('Real');
ylabel('Imaginary');