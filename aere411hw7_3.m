% Logan Prust - Aer E 411 - Homework 7 - Problem 3

clear,clc;

% define constants
T0 = 220;
gammac = 1.4;
Cpc = 1000;
gammat = 1.3;
Cpt = 1240;
pib = 0.98;
etab = 0.99;
etaAB = 0.96;
etam = 0.99;
ec = 0.92;
et = 0.90;
h = 4.5E7;
piAB = 1;

% Part (b), afterburning case

n = 1000;
M0 = linspace(1,3.5,n);

pin = 0.95;
taulambda = 7;
taulambdaAB = 8.5;
pic = 15;

for i=1:n
    taur = 1+(gammac-1)/2*M0(i)^2;
    pir = taur^(gammac/(gammac-1));
    pid = 1-0.015*M0(i)^2;
    tauc = pic^((gammac-1)/gammac/ec);
    f = (tauc*taur-taulambda)/(taulambda-etab*h/Cpc/T0);
    taut = 1-1/etam/(1+f)*taur/taulambda*(tauc-1);
    pit = taut^(gammat*et/(gammat-1));
    fAB = (1+f)*(taulambdaAB-taut*taulambda)/(etaAB*h/Cpc/T0-taulambdaAB);
    Rc = (gammac-1)/gammac*Cpc;
    Rt = (gammat-1)/gammat*Cpt;
    a0 = sqrt(gammac*Rc*T0);
    Pratio = pin*piAB*pit*pib*pic*pid*pir;
    gammaAB = gammat;
    RAB = Rt;
    this = sqrt(2*taulambdaAB/(gammaAB-1)*gammaAB*RAB/gammac/Rc*(1-Pratio^((1-gammaAB)/gammaAB)));
    spcthrust_b_AB(i) = a0*((1+f+fAB)*this - M0(i));
    S_b_AB(i) = (f+fAB)/spcthrust_b_AB(i);
end

% Part (b), non-afterburning case

pin = 0.98;
taulambdaAB = taulambda;

for i=1:n
    taur = 1+(gammac-1)/2*M0(i)^2;
    pir = taur^(gammac/(gammac-1));
    pid = 1-0.015*M0(i)^2;
    tauc = pic^((gammac-1)/gammac/ec);
    f = (tauc*taur-taulambda)/(taulambda-etab*h/Cpc/T0);
    taut = 1-1/etam/(1+f)*taur/taulambda*(tauc-1);
    pit = taut^(gammat*et/(gammat-1));
    fAB = 0;
    Rc = (gammac-1)/gammac*Cpc;
    Rt = (gammat-1)/gammat*Cpt;
    a0 = sqrt(gammac*Rc*T0);
    Pratio = pin*piAB*pit*pib*pic*pid*pir;
    gammaAB = gammat;
    RAB = Rt;
    this = sqrt(2*taulambdaAB/(gammaAB-1)*gammaAB*RAB/gammac/Rc*(1-Pratio^((1-gammaAB)/gammaAB)));
    spcthrust_b(i) = a0*((1+f+fAB)*this - M0(i));
    S_b(i) = (f+fAB)/spcthrust_b(i);
end

% Part (c), afterburning case

pic = linspace(4,73,n);

pin = 0.95;
taulambda = 7;
taulambdaAB = 8.5;
M0 = 2.2;

for i=1:n
    taur = 1+(gammac-1)/2*M0^2;
    pir = taur^(gammac/(gammac-1));
    pid = 1-0.015*M0^2;
    tauc = pic(i)^((gammac-1)/gammac/ec);
    f = (tauc*taur-taulambda)/(taulambda-etab*h/Cpc/T0);
    taut = 1-1/etam/(1+f)*taur/taulambda*(tauc-1);
    pit = taut^(gammat*et/(gammat-1));
    fAB = (1+f)*(taulambdaAB-taut*taulambda)/(etaAB*h/Cpc/T0-taulambdaAB);
    Rc = (gammac-1)/gammac*Cpc;
    Rt = (gammat-1)/gammat*Cpt;
    a0 = sqrt(gammac*Rc*T0);
    Pratio = pin*piAB*pit*pib*pic(i)*pid*pir;
    gammaAB = gammat;
    RAB = Rt;
    this = sqrt(2*taulambdaAB/(gammaAB-1)*gammaAB*RAB/gammac/Rc*(1-Pratio^((1-gammaAB)/gammaAB)));
    spcthrust_c_AB(i) = a0*((1+f+fAB)*this - M0);
    S_c_AB(i) = (f+fAB)/spcthrust_c_AB(i);
end

% Part (c), non-afterburning case

pin = 0.98;
taulambdaAB = taulambda;

for i=1:n
    taur = 1+(gammac-1)/2*M0^2;
    pir = taur^(gammac/(gammac-1));
    pid = 1-0.015*M0^2;
    tauc = pic(i)^((gammac-1)/gammac/ec);
    f = (tauc*taur-taulambda)/(taulambda-etab*h/Cpc/T0);
    taut = 1-1/etam/(1+f)*taur/taulambda*(tauc-1);
    pit = taut^(gammat*et/(gammat-1));
    fAB = 0;
    Rc = (gammac-1)/gammac*Cpc;
    Rt = (gammat-1)/gammat*Cpt;
    a0 = sqrt(gammac*Rc*T0);
    Pratio = pin*piAB*pit*pib*pic(i)*pid*pir;
    gammaAB = gammat;
    RAB = Rt;
    this = sqrt(2*taulambdaAB/(gammaAB-1)*gammaAB*RAB/gammac/Rc*(1-Pratio^((1-gammaAB)/gammaAB)));
    spcthrust_c(i) = a0*((1+f+fAB)*this - M0);
    S_c(i) = (f+fAB)/spcthrust_c(i);
end

% plot results

figure(1);
hold on;
M0 = linspace(1,3.5,n);
plot(M0,spcthrust_b_AB,'r');
plot(M0,spcthrust_b);
legend('Afterburning','Non-Afterburning');
title('Specific Thrust vs Mach Number');
xlabel('Mach Number');
ylabel('Specific Thrust (m/s)');

figure(2);
hold on;
plot(M0,S_b_AB,'r');
plot(M0,S_b);
legend('Afterburning','Non-Afterburning');
title('Specific Fuel Consumption vs Mach Number');
xlabel('Mach Number');
ylabel('Specific Fuel Consumption (s/m)');

pic1 = 72.45;
pic2 = 71.16;

figure(3);
plot(pic,spcthrust_c_AB);
axis([4 pic1 0 1100]);
title('Specific Thrust vs Compressor Pressure Ratio (Afterburning Engine)');
xlabel('Compressor Total Pressure Ratio');
ylabel('Specific Thrust (m/s)');

figure(4);
plot(pic,S_c_AB);
axis([4 pic1 0 0.0002]);
title('Specific Fuel Consumption vs Compressor Pressure Ratio (Afterburning Engine)');
xlabel('Compressor Total Pressure Ratio');
ylabel('Specific Fuel Consumption (s/m)');

figure(5);
plot(pic,spcthrust_c);
axis([4 pic2 0 800]);
title('Specific Thrust vs Compressor Pressure Ratio (Non-Afterburning Engine)');
xlabel('Compressor Total Pressure Ratio');
ylabel('Specific Thrust (m/s)');

figure(6);
plot(pic,S_c);
axis([4 pic2 0 0.00003]);
title('Specific Fuel Consumption vs Compressor Pressure Ratio (Non-Afterburning Engine)');
xlabel('Compressor Total Pressure Ratio');
ylabel('Specific Fuel Consumption (s/m)');