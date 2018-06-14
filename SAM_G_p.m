function[p]=SAM_G_p(p1,T,T1,g0,a,R)
p=p1*(T/T1).^(-g0/a/R);
end