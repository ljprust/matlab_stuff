clear;
input('Max HP ');
hpmax=ans;
input('Current HP ');
hpcur=ans;
input('Catch Rate ');
pokerate=ans;
input('Ball Bonus ');
ball=ans;
input('Status Multiplier ');
stat=ans;
the=(3*hpmax-2*hpcur)*pokerate*ball/(3*hpmax)*stat;
the2=1048560*(16711680/the)^(1/4);
disp(the);
disp(the2);