% Logan Prust - Astro 346 Project
clear,clc;
planetinfo=dlmread('planetinfo.txt');
pmass=planetinfo(:,1);
pradius=planetinfo(:,2);
sradius=planetinfo(:,3);
Teff=planetinfo(:,4);
a=planetinfo(:,5);
e=planetinfo(:,6);
smass=planetinfo(:,8);
Eradius=6.378e6;
pmass=pmass*1.8986e27;
a=a*1.5e11;
pradius=pradius*69911000;
sradius=sradius*6.955e8;
smass=smass*1.989e30;
albedo=zeros(1,length(planetinfo));
moon=zeros(1,length(planetinfo));
L=zeros(1,length(planetinfo));
peri=zeros(1,length(planetinfo));
apo=zeros(1,length(planetinfo));
tmoonhab=zeros(1,length(planetinfo));
imoonhab=zeros(1,length(planetinfo));
HZinner=zeros(1,length(planetinfo));
HZouter=zeros(1,length(planetinfo));
innerpdist=zeros(1,length(planetinfo));
outerpdist=zeros(1,length(planetinfo));
phab=zeros(1,length(planetinfo));
tHZinner=zeros(1,length(planetinfo));
tHZouter=zeros(1,length(planetinfo));
innertdist=zeros(1,length(planetinfo));
outertdist=zeros(1,length(planetinfo));
iHZinner=zeros(1,length(planetinfo));
iHZouter=zeros(1,length(planetinfo));
inneridist=zeros(1,length(planetinfo));
outeridist=zeros(1,length(planetinfo));
pdensity=zeros(1,length(planetinfo));
for i=1:length(planetinfo)
    pdensity(i)=pmass(i)./(4/3*pi*pradius(i)^3);
    if pdensity(i)>=4000
        albedo(i)=.25;
    elseif pradius(i)<2.5*Eradius
        albedo(i)=.99;
    else
        moon(i)=1;
    end
    L(i)=4*pi*(5.67e-8)*sradius(i)^2*Teff(i)^4;
    peri(i)=a(i)*(1-e(i)^2)/(1+e(i));
    apo(i)=a(i)*(1-e(i)^2)/(1-e(i));
    if moon(i)==0
        tmoonhab(i)=0/0;
        imoonhab(i)=0/0;
        HZinner(i)=sqrt((1-albedo(i))*L(i)/16/pi/(5.67e-8)/.61/373.15^4);
        HZouter(i)=sqrt((1-albedo(i))*L(i)/16/pi/(5.67e-8)/.61/273.15^4);
        innerpdist(i)=HZinner(i)-peri(i);
        if innerpdist(i)<0
            innerpdist(i)=0;
        end
        outerpdist(i)=apo(i)-HZouter(i);
        if outerpdist(i)<0
            outerpdist(i)=0;
        end
        phab(i)=max(innerpdist(i),outerpdist(i))*(-1);
    else
        phab(i)=0/0;
        tHZinner(i)=sqrt((1-.25)*L(i)/16/pi/(5.67e-8)/.61/373.15^4);
        tHZouter(i)=sqrt((1-.25)*L(i)/16/pi/(5.67e-8)/.61/273.15^4);
        innertdist(i)=tHZinner(i)-peri(i);
        if innertdist(i)<0
            innertdist(i)=0;
        end
        outertdist(i)=apo(i)-tHZouter(i);
        if outertdist(i)<0
            outertdist(i)=0;
        end
        tmoonhab(i)=max(innertdist(i),outertdist(i))*(-1);
        iHZinner(i)=sqrt((1-.99)*L(i)/16/pi/(5.67e-8)/.61/373.15^4);
        iHZouter(i)=sqrt((1-.99)*L(i)/16/pi/(5.67e-8)/.61/273.15^4);
        inneridist(i)=iHZinner(i)-peri(i);
        if inneridist(i)<0
            inneridist(i)=0;
        end
        outeridist(i)=apo(i)-iHZouter(i);
        if outeridist(i)<0
            outeridist(i)=0;
        end
        imoonhab(i)=max(inneridist(i),outeridist(i))*(-1);
    end
end
phab=phab/(1.5e11);
tmoonhab=tmoonhab/(1.5e11);
imoonhab=imoonhab/(1.5e11);
smass=smass/(1.989e30);
output=[smass,phab',tmoonhab',imoonhab'];
dlmwrite('planetresults.txt',output,',');