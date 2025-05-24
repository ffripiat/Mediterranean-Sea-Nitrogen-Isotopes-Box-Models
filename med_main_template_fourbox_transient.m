%%%%%%%     MODEL MAIN TEMPLATE    %%%%%%%

%A
%1 = NO3_surf_W;2 = NO3_deep_W;3 = NO3_surf_E;4 = NO3_deep_E,; 
%5 = PO4_surf_W;6 = PO4_deep_W;7 = PO4_surf_E;8 = PO4_deep_E; 
%9 = 15NO3_surf_W;10 = 15NO3_deep_W,11 = 15NO3_surf_E,12 = 15NO3_deep_E; 
%13 = d15N_surf_W;14 = d15N_deep_W,15 = d15N_surf_E;16 = d15N_deep_E; 
%17 = NO3_total; 18 = PO4_total; 19 = d15N_total
%20 = external N input (mmol yr-1); 21 = benthic denitrification (mmol yr-1);
%22 = external N d15N input; 23 = external P input (mmol yr-1);
%24 = NO3 concentration in the Atlantic inflow; 25 = PO4 concentration in the Atlantic inflow;
%26 = nitrate d15N in the Atlantic inflow; 27 = w; 28 = mD; 29 = mS 
%30 = External N supply by river (mmol yr-1); 31 = External P supply by river (mmol yr-1)
%32 = external N input (Tg yr-1); 33 = benthic denitrification (Tg yr-1);
%34 = external P input (Tg yr-1); 35 = external N/P input; 36 = Mediterrenean Sea N/P;
%37 = w/mD, 38 = mS (same as 29); 39 = export production from the Mediterreanean sea (Tg yr-1)
%40 = Nitrate supply from the Atlantic (Tg N yr-1); 41 = External N supply by river (Tg N yr-1)
%42 = External P supply by river (Tg P yr-1)

%==============================================
%INITIAL CONDITIONS and time span of simulation
%===============================================

RrefN        =  3.6782e-3;

tspan = (0:1:70); %70 years 

%=======================
%Parameters declariation
%=======================
%fixed model parameters (i.e., not changing)
V1           = 6.1600e+13; %volume surface box in the western basin in m3 
V2           = 1.7884e+15; %volume deep box in the western basin in m3
V3           = 1.3353e+14; %volume surface box in the eastern basin in m3 
V4           = 3.8765e+15; %volume deep box in the easter basin in m3
eps_deni     = 0; %isotope effect in p.mil for benthic denitrification
K_upt        = 0.9999; %degree of nitrate consumption in surface (1 = full consumption)
K_rem        = 1; %degree of remineralization in deep box (1 = full remineralization)
eps          = 5; %isotope effect in p.mil for nitrate assimilation

%create a matrix of zero value with a number of rows being equal to the
%sensitivity matrix and a number of column being equal to the model
%variables. To be filled with model output
x0 = zeros(1,12);
xend = zeros(length(sensitivity(:,1)),length(x0));

%the following loop generated one scenario (with a given set of model parameters) for each row of the sensitivity matrix
%in addition, it generate the initial box conditions in the range given for the Atlantic Ocean
for i=1:length(sensitivity(:,1))
x0 = [0.001 sensitivity(i,5) 0.001 sensitivity(i,5) 0.001 sensitivity(i,6) 0.001 sensitivity(i,6) ((((sensitivity(i,7))./1000)+1).*RrefN ).*0.001 ((((sensitivity(i,7))./1000)+1).*RrefN ).*sensitivity(i,5) ((((sensitivity(i,7))./1000)+1).*RrefN ).*0.001 ((((sensitivity(i,7))./1000)+1).*RrefN ).*sensitivity(i,5)];    
input_N     = sensitivity(i,1);
deni        = sensitivity(i,2);
d15N_input  = sensitivity(i,3);
input_P     = sensitivity(i,4);
N15_input   = ((((d15N_input)/1000)+1)*RrefN ).*input_N;   
NO3_atl     = sensitivity(i,5); 
PO4_atl     = sensitivity(i,6); 
d15N_atl    = sensitivity(i,7); 
N15_atl     =((((d15N_atl)/1000)+1)*RrefN).*NO3_atl; 
M           = sensitivity(i,8); 
M_h         = sensitivity(i,9); 
M_SD        = sensitivity(i,10);
river_input_N = sensitivity(i,11);
river_input_N_W = river_input_N.*0.55;
river_input_N_E = river_input_N.*0.45;
river_input_P = sensitivity(i,12);
river_input_P_W = river_input_P.*0.55;
river_input_P_E = river_input_P.*0.45;
d15N_river    = sensitivity(i,13); 
N15_river_W   = ((((d15N_river)/1000)+1)*RrefN ).*river_input_N_W;  
N15_river_E   = ((((d15N_river)/1000)+1)*RrefN ).*river_input_N_E; 

%=================
%SOLVING THE ODE
%=================
[t,x]       = ode15s(@(t,x)med_fourbox_mixingWE_ode(t,x,V1,V2,V3,V4,NO3_atl,PO4_atl,N15_atl,M,M_SD,M_h,eps_deni,K_upt,K_rem,eps, input_N,deni,N15_input,input_P,river_input_N_W,river_input_N_E,river_input_P_W,river_input_P_E,N15_river_W,N15_river_E),tspan,x0,[]);

% to convert 15N concentration in delta values
d15N        = zeros(length(tspan),2);
d15N(:,1)   = (((x(:,9)./x(:,1))./RrefN )-1)*1000;
d15N(:,2)   = (((x(:,10)./x(:,2))./RrefN )-1)*1000;
d15N(:,3)   = (((x(:,11)./x(:,3))./RrefN )-1)*1000;
d15N(:,4)   = (((x(:,12)./x(:,4))./RrefN )-1)*1000;
% to compute weighted values for the Mediterreanean sea
total       = zeros(length(tspan),3);
total(:,1)  = (x(:,1).*V1+x(:,2).*V2+x(:,3).*V3+x(:,4).*V4)./(V1+V2+V3+V4);% nitrate concentration
total(:,2)  = (x(:,5).*V1+x(:,6).*V2+x(:,7).*V3+x(:,8).*V4)./(V1+V2+V3+V4);% phosphate concentration
total(:,3)  = (d15N(:,1).*x(:,1).*V1+d15N(:,2).*x(:,2).*V2+d15N(:,3).*x(:,3).*V3+d15N(:,4).*x(:,4).*V4)./((V1+V2+V3+V4).*total(:,1));% nitrate d15N
% xend is the matrix with only final solution (i.e., last row in x, when there is no more varation with time)
xend(i,1:12)=x(length(tspan),1:12);
xend(i,13:16)  =d15N(length(tspan),1:4);
xend(i,17:19)  =total(length(tspan),1:3);
end
% to create the final matrix A
A = [xend,sensitivity]; % merging the model outputs (only final solution) with sensitivity matrix (i.e., model input parameters)
%offline calculation
A(:,32) = sensitivity(:,1).*14./10^15; % external N input (Tg yr-1)
A(:,33) = sensitivity(:,2).*14./10^15; % benthic denitrification (Tg yr-1)
A(:,34) = sensitivity(:,4).*31./10^15; % external P input (Tg yr-1)
A(:,35) = sensitivity(:,1)./sensitivity(:,4); %external N/P input
A(:,36) = A(:,17)./A(:,18); %Mediterrenean Sea N/P
A(:,37) = A(:,28)./A(:,27); %w/mD
A(:,38) = sensitivity(:,10); % mS
A(:,39) = (A(:,38).*K_upt.*A(:,2)+A(:,27).*K_upt.*A(:,24)+A(:,38).*K_upt.*A(:,4)).*14./10^15; %export production from the Mediterreanean sea (Tg yr-1)
A(:,40) = A(:,24).*A(:,27).*14./10^15; %Nitrate supply from the Atlantic (Tg N yr-1)
A(:,41) = A(:,29).*14./10^15; % External N supply by river (Tg N yr-1) 
A(:,42) = A(:,30).*31./10^15; % External P supply by river (Tg P yr-1) 
