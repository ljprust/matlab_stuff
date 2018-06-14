function dthetadt = func(t,theta)
dthetadt = [theta(2); -theta(1)*theta(2)-1];