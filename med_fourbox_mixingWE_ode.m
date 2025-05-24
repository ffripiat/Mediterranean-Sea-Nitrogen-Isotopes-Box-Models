function dx = med_fourbox_mixing_ode(t,x,V1,V2,V3,V4,NO3_atl,PO4_atl,N15_atl,M,M_SD, M_h,eps_deni,K_upt,K_rem,eps, input_N,deni,N15_input,input_P, river_input_N_W,river_input_N_E,river_input_P_W,river_input_P_E,N15_river_W,N15_river_E)

%================
% DEFINE YOUR SYSTEM
%================
dx=zeros(9,1);

%NO3 conc
dx(1) = M/V1*NO3_atl-M/V1*x(1)+M_SD/V1*(x(2)-x(1))-K_upt*(M/V1*NO3_atl+M_SD/V1*x(2));
dx(2) = M/V2*x(4)-M/V2*x(2)+M_h/V2*(x(4)-x(2))+M_SD/V2*(x(1)-x(2))+K_upt*K_rem*(M/V2*NO3_atl+M_SD/V2*x(2))+(input_N*0.34)/V2-(deni*0.34)/V2+river_input_N_W/V2;
dx(3) = M/V3*x(1)-M/V3*x(3)+M_SD/V3*(x(4)-x(3))-K_upt*(M/V3*x(1)+M_SD/V3*x(4));
dx(4) = M/V4*x(3)-M/V4*x(4)+M_h/V4*(x(2)-x(4))+M_SD/V4*(x(3)-x(4))+K_upt*K_rem*(M/V4*x(1)+M_SD/V4*x(4))+(input_N*0.66)/V4-(deni*0.66)/V4+river_input_N_E/V4;

%PO4 conc
dx(5) = M/V1*PO4_atl-M/V1*x(5)+M_SD/V1*(x(6)-x(5))-K_upt*(M/V1*PO4_atl+M_SD/V1*x(6));
dx(6) = M/V2*x(8)-M/V2*x(6)+M_h/V2*(x(8)-x(6))+M_SD/V2*(x(5)-x(6))+K_upt*K_rem*(M/V2*PO4_atl+M_SD/V2*x(6))+(input_P*0.34)/V2+river_input_P_W/V2;
dx(7) = M/V3*x(5)-M/V3*x(7)+M_SD/V3*(x(8)-x(7))-K_upt*(M/V3*x(5)+M_SD/V3*x(7));
dx(8) = M/V4*x(7)-M/V4*x(8)+M_h/V4*(x(6)-x(8))+M_SD/V4*(x(7)-x(8))+K_upt*K_rem*(M/V4*x(5)+M_SD/V4*x(7))+(input_P*0.66)/V4+river_input_P_E/V4;

%15NO3 conc
dx(9) = M/V1*NO3_atl*N15_atl/NO3_atl-M/V1*x(1)*x(9)/x(1)+M_SD/V1*(x(2)*x(10)/x(2)-x(1)*x(9)/x(1))-K_upt*M/V1*NO3_atl*N15_atl/NO3_atl*(1-(1-K_upt)^(1-eps/1000))/K_upt-K_upt*M_SD/V1*x(2)*x(10)/x(2)*(1-(1-K_upt)^(1-eps/1000))/K_upt;
dx(10) = M/V2*x(4)*x(12)/x(4)-M/V2*x(2)*x(10)/x(2)+M_h/V2*(x(4)*x(12)/x(4)-x(2)*x(10)/x(2))+M_SD/V2*(x(1)*x(9)/x(1)-x(2)*x(10)/x(2))+K_upt*K_rem*M/V2*NO3_atl*N15_atl/NO3_atl*(1-(1-K_upt)^(1-eps/1000))/K_upt+K_upt*K_rem*M_SD/V2*x(2)*x(10)/x(2)*(1-(1-K_upt)^(1-eps/1000))/K_upt+(input_N*0.34)/V2*(N15_input*0.34)/(input_N*0.34)-(deni*0.34)/V2*x(10)/x(2)*(1-eps_deni/1000)+river_input_N_W/V2*N15_river_W/river_input_N_W;
dx(11) = M/V3*x(1)*x(9)/x(1)-M/V3*x(3)*x(11)/x(3)+M_SD/V3*(x(4)*x(12)/x(4)-x(3)*x(11)/x(3))-K_upt*M/V3*x(1)*x(9)/x(1)*(1-(1-K_upt)^(1-eps/1000))/K_upt-K_upt*M_SD/V3*x(4)*x(12)/x(4)*(1-(1-K_upt)^(1-eps/1000))/K_upt;
dx(12) = M/V4*x(3)*x(11)/x(3)-M/V4*x(4)*x(12)/x(4)+M_h/V4*(x(2)*x(10)/x(2)-x(4)*x(12)/x(4))+M_SD/V4*(x(3)*x(11)/x(3)-x(4)*x(12)/x(4))+K_upt*K_rem*M/V4*x(1)*x(9)/x(1)*(1-(1-K_upt)^(1-eps/1000))/K_upt+K_upt*K_rem*M_SD/V4*x(4)*x(12)/x(4)*(1-(1-K_upt)^(1-eps/1000))/K_upt+(input_N*0.66)/V4*(N15_input*0.66)/(input_N*0.66)-(deni*0.66)/V4*x(12)/x(4)*(1-eps_deni/1000)+river_input_N_E/V4*N15_river_E/river_input_N_E;


