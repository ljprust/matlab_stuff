function[rohI]=SAM_I_roh(rohG,g0,R,TG,hI,I)
rohI=rohG(length(rohG))*exp(-g0/R/TG(length(TG))*(hI-I)*1000);
end