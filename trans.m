function trusstrans = trans(t)

T = zeros(6);

T(1,1) = cos(t);
T(1,2) = sin(t);
T(2,1) = -T(1,2);
T(2,2) = T(1,1);
T(3,3) = 1;
T(4,4) = T(1,1);
T(4,5) = T(1,2);
T(5,4) = T(2,1);
T(5,5) = T(4,4);
T(6,6) = 1;

trusstrans = T;