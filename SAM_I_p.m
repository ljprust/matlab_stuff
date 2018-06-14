function[pI]=SAM_I_p(pG,g0,R,TG,hI,I)
pI=pG(length(pG))*exp(-g0/R/TG(length(TG))*(hI-I)*1000);
end