tic
clear;
clc;

% IMPORT CHARACTERISTIC TEMPERATURES for P_FET=1, P_Shunt=0, P_PCB=0
opts1 = spreadsheetImportOptions("NumVariables", 4);
opts1.Sheet = "Char_FET";
opts1.DataRange = "A2:D27";
opts1.VariableNames = ["t1", "T_FET", "T_Shunt", "T_PCB"];
opts1.SelectedVariableNames = ["t1", "T_FET", "T_Shunt", "T_PCB"];
opts1.VariableTypes = ["double", "double", "double", "double"];
tbl1 = readtable("ECU.xlsx", opts1, "UseExcel", false);
t1 = tbl1.t1;
T_FET_100 = tbl1.T_FET;
T_Shunt_100 = tbl1.T_Shunt;
T_PCB_100 = tbl1.T_PCB;
clear opts tbl1

% IMPORT CHARACTERISTIC TEMPERATURES for P_FET=0, P_Shunt=1, P_PCB=0
opts2 = spreadsheetImportOptions("NumVariables", 4);
opts2.Sheet = "Char_Shunt";
opts2.DataRange = "A2:D27";
opts2.VariableNames = ["t1", "T_FET", "T_Shunt", "T_PCB"];
opts2.SelectedVariableNames = ["t1", "T_FET", "T_Shunt", "T_PCB"];
opts2.VariableTypes = ["double", "double", "double", "double"];
tbl2 = readtable("ECU.xlsx", opts2, "UseExcel", false);
t2 = tbl2.t1;
T_FET_010 = tbl2.T_FET;
T_Shunt_010 = tbl2.T_Shunt;
T_PCB_010 = tbl2.T_PCB;
clear opts tbl2

% IMPORT THERMAL PROPERTIES
opts3 = spreadsheetImportOptions("NumVariables", 1);
opts3.Sheet = "Properties";
opts3.DataRange = "J19:J19";
opts3.VariableNames = "Cth";
opts3.SelectedVariableNames = "Cth";
opts3.VariableTypes = "double";
tbl3 = readtable("ECU.xlsx", opts3, "UseExcel", false);
C_T = tbl3.Cth;
clear opts tbl3

% IMPORT TRANSIENT POWER INPUTS 
opts4 = spreadsheetImportOptions("NumVariables", 3);
opts4.Sheet = "Input";
opts4.DataRange = "A2:C27";
opts4.VariableNames = ["t1", "P1", "P2"];
opts4.SelectedVariableNames = ["t1", "P1", "P2"];
opts4.VariableTypes = ["double", "double", "double"];
tbl4 = readtable("ECU.xlsx", opts4, "UseExcel", false);
t_tr = tbl4.t1;
P1_tr = tbl4.P1;
P2_tr = tbl4.P2;
clear opts tbl4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CALCULATE THERMAL CAPACITANCE
PoC_T=1/C_T;

T0=20;
% Linear Temperature curve
TL=T0+PoC_T*t1;
% Deviations
TD1_100=abs(T_FET_100-TL);
TD2_100=abs(T_Shunt_100-TL);
TD3_100=abs(T_PCB_100-TL);

TD1_010=abs(T_FET_010-TL);
TD2_010=abs(T_Shunt_010-TL);
TD3_010=abs(T_PCB_010-TL);

% ESTIMATE Rth, k via Non-linear Optimization
T_model=@(MC) (MC(1)*(1-exp(-MC(2)*t1)));
initial_guess=[0.1 0.01];
Loss=@(T_model,MC,y_true) sum((T_model(MC)-y_true).^2,'all');

%TM1_100
MC_fit=fminsearch(@(MC) Loss(T_model,MC,TD1_100),initial_guess);
R1B1=MC_fit(1) * sign(T_FET_100(end)-TL(end)); % For +/- values of Rth: sign(TB - TL)
k1_100=MC_fit(2);

%TM2_100
MC_fit=fminsearch(@(MC) Loss(T_model,MC,TD2_100),initial_guess);
R1B2=MC_fit(1) * sign(T_Shunt_100(end)-TL(end));
k2_100=MC_fit(2);

%TM3_100
MC_fit=fminsearch(@(MC) Loss(T_model,MC,TD3_100),initial_guess);
R1B3=MC_fit(1) * sign(T_PCB_100(end)-TL(end));
k3_100=MC_fit(2);
%------------------------
%TM1_010
MC_fit=fminsearch(@(MC) Loss(T_model,MC,TD1_010),initial_guess);
R2B1=MC_fit(1) * sign(T_FET_010(end)-TL(end)); % For +/- values of Rth: sign(TB - TL)
k1_010=MC_fit(2);

%TM2_010
MC_fit=fminsearch(@(MC) Loss(T_model,MC,TD2_010),initial_guess);
R2B2=MC_fit(1) * sign(T_Shunt_010(end)-TL(end));
k2_010=MC_fit(2);

%TM3_010
MC_fit=fminsearch(@(MC) Loss(T_model,MC,TD3_010),initial_guess);
R2B3=MC_fit(1) * sign(T_PCB_010(end)-TL(end));
k3_010=MC_fit(2);


K_FET = k1_100;
K_Shunt = k2_010;
K_PCB = 0.5*(k3_100 + k3_010);

[TLtr,TF1,TF2,TF3] = f_ECU(P1_tr,P2_tr,t_tr,R1B1,R1B2,R1B3,R2B1,R2B2,R2B3,K_FET,K_Shunt,K_PCB,C_T);

figure
plot(t_tr,TF1,t_tr,TF2,'--',t_tr,TF3,'-+')
legend('Body 1','Body 2','Body 3')
title('Model Temperature')

toc