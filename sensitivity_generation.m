clear all
size=1000000; %number of scenario (10,000 scenario take around 1hour to run)

%to generate random number in a given range

%here you must chose which external source of low-δ15N nitrate: N2 fixation, atmospheric deposition, or partial DON breakdown 
a = 0.001./14.*10^15;% range between 0.001 and 10 Tg N yr-1, external source input
b = 10./14.*10^15;
r = (b-a).*rand(size,1) + a;
sensitivity(:,1) = r; % external N source input in mmol yr-1
a = 16; %non-relevant for N2 fixation, (a = 100 and b = 160) for atmospheric deposition, (a = 16 and b = 25) for partial DON breakdown
b = 25;
p = (b-a).*rand(size,1) + a;
sensitivity(:,2) = 0; % to run with (sensitivity(:,2) = r or without (sensitivity(:,2) = 0) benthic denitrification (mmol yr-1), assuming that benthic denitrifcation = external N input 
sensitivity(:,4) = 0.0000000000000001;% external P source input in mmol yr-1; (sensitivity(:,4) = r./p) for atmospheric deposition or partial DON breakdown; (sensitivity(:,4) = 0.0000000000000001) for N2 fixation

%d15N of the external source
a = -2;%(a = -2 and b = 0) for N2 fixation, (a = -4 and b = -2) for atmospheric deposition, (a = 2 and b = 4) for partial DON breakdown
b = 0;
r = (b-a).*rand(size,1) + a;
sensitivity(:,3) = r; % d15N of the external source

%nitrate concentration in the Atlantic inflow
a = 0.5;%
b = 5;
r = (b-a).*rand(size,1) + a;
sensitivity(:,5) = r; % nitrate concentration in the Atlantic inflow in mmol m-3

%phosphate concentration in the Atlantic inflow
a = 0.03;%
b = 0.35;
r = (b-a).*rand(size,1) + a;
sensitivity(:,6) = r; % phosphate concentration in the Atlantic inflow in mmol m-3

%nitrate d15N in the Atlantic inflow
a = 4;%
b = 5;
r = (b-a).*rand(size,1) + a;
sensitivity(:,7) = r; % nitrate d15N in the Atlantic inflow in per mil

%advective water flow at the Gibraltar Strait (ω – anti-estuarine overturning circulation)
a = 5e5*3600*24*365;% from 0.5 to 1 Sv 
b = 10e5*3600*24*365;
r = (b-a).*rand(size,1) + a;
sensitivity(:,8) = r; % Advective water flow at the Gibraltar Strait in m3 s-1

%mixing between the two deep Mediterranean basins (mD) 
a = 0e5*3600*24*365;% from 0.0 to 5 Sv 
b = 50e5*3600*24*365;
r = (b-a).*rand(size,1) + a;
sensitivity(:,9) = r; % mixing between the two deep Mediterranean basins in m3 s-1 

%mixing between surface and deep waters in each basin (mS) 
a = 0e5*3600*24*365;% from 0.0 to 5 Sv 
b = 50e5*3600*24*365;%
r = (b-a).*rand(size,1) + a; 
sensitivity(:,10) = r;% mixing between surface and deep waters in each basin in m3 s-1 

%External N supply by river
a = 0.0000000000000001% (a = 0.0000000000000001 and b = 0.0000000000000001) for no river supply; (a = 0.1/14.*10^15 and b = 0.1/14.*10^15) for natural river supply; and (a = 0.6/14.*10^15 and b = 0.6/14.*10^15) for anthropogenic river supply 
b = 0.0000000000000001;
r = (b-a).*rand(size,1) + a;
sensitivity(:,11) = r; 
sensitivity(:,12) = r./6; % river P input in mmol yr-1 ; (r./6) for natural river and (r./13.2) for anthropogenic river 

%River nitrate d15N
a = 5;%
b = 15;
r = (b-a).*rand(size,1) + a;
sensitivity(:,13) = r; % River nitrate d15N in per mil