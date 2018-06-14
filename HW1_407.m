% Logan Prust - Math 407 - HW 2 - Problem 2.1.1

A=[7+3j,-1,2;3,-5,3;1j,2,5+8j];
vals=eig(A); % get eigenvalues
xvals=real(vals); % put real values on the x-axis
yvals=imag(vals); % put imaginary values on the y-axis

rowRadii = sum(abs(A),2) - abs(diag(A)); % sum across rows, neglecting the diagonal elements

figure;
hold on;
plot(xvals,yvals,'k*'); % plot eigenvalues
plot_circle(gca,real(A(1,1)),imag(A(1,1)),rowRadii(1));
plot_circle(gca,real(A(2,2)),imag(A(2,2)),rowRadii(2)); % plot disks
plot_circle(gca,real(A(3,3)),imag(A(3,3)),rowRadii(3));
hold off;
axis equal; % set aspect ratio
draw_axes(gca); % draw axes on plot
title('Row Disks for Problem 2.1.1');
xlabel('Real');
ylabel('Imaginary');

At=transpose(A); % flip matrix
colRadii = sum(abs(At),2) - abs(diag(At)); % sum across rows (which were the columns of the original matrix)

figure;
hold on;
plot(xvals,yvals,'k*'); % plot eignvalues
for ctr=1:3
plot_circle(gca,real(At(ctr,ctr)),imag(At(ctr,ctr)),colRadii(ctr)); % plot disks
end
hold off;
axis equal; % set aspect ratio
draw_axes(gca); % draw axes on plot
title('Column Disks for Problem 2.1.1');
xlabel('Real');
ylabel('Imaginary');