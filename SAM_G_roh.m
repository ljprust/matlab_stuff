function[rohG]=SAM_G_roh(roh1,T,T1,g0,a,R)
rohG=roh1*(T/T1).^(-g0/a/R-1);
end