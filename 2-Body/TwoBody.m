tic
clear;
clc;

% IMPORT CHARACTERISTIC TEMPERATURES for P1=1, P2=0
opts1 = spreadsheetImportOptions("NumVariables", 3);
opts1.Sheet = "T1";
opts1.DataRange = "A2:C102";
opts1.VariableNames = ["t", "TB1_10", "TB2_10"];
opts1.SelectedVariableNames = ["t", "TB1_10", "TB2_10"];
opts1.VariableTypes = ["double", "double", "double"];
tbl1 = readtable("2_Body.xlsx", opts1, "UseExcel", false);

t = tbl1.t;
TB1_10 = tbl1.TB1_10;
TB2_10 = tbl1.TB2_10;
% IMPORT CHARACTERISTIC TEMPERATURES for P1=0, P2=1
opts2 = spreadsheetImportOptions("NumVariables", 3);
opts2.Sheet = "T2";
opts2.DataRange = "A2:C102";
opts2.VariableNames = ["t", "TB1_01", "TB2_01"];
opts2.SelectedVariableNames = ["t", "TB1_01", "TB2_01"];
opts2.VariableTypes = ["double", "double", "double"];
tbl2 = readtable("2_Body.xlsx", opts2, "UseExcel", false);

t2 = tbl2.t;
TB1_01 = tbl2.TB1_01;
TB2_01 = tbl2.TB2_01;
% IMPORT THERMAL PROPERTIES
opts3 = spreadsheetImportOptions("NumVariables", 8);

opts3.Sheet = "properties";
opts3.DataRange = "A2:H2";

opts3.VariableNames = ["V1", "rho1", "cp1", "k1", "V2", "rho2", "cp2", "k2"];
opts3.SelectedVariableNames = ["V1", "rho1", "cp1", "k1", "V2", "rho2", "cp2", "k2"];
opts3.VariableTypes = ["double", "double", "double", "double", "double", "double", "double", "double"];

tbl3 = readtable("2_Body.xlsx", opts3, "UseExcel", false);

V1 = tbl3.V1;
rho1 = tbl3.rho1;
cp1 = tbl3.cp1;
k1_10 = tbl3.k1;
V2 = tbl3.V2;
rho2 = tbl3.rho2;
cp2 = tbl3.cp2;
k2_10 = tbl3.k2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear opts1 tbl1
clear opts2 tbl2
clear opts3 tbl3

% IMPORT TRANSIENT POWER INPUTS 
opts4 = spreadsheetImportOptions("NumVariables", 3);
opts4.Sheet = "Inputs";
opts4.DataRange = "A2:C27";
opts4.VariableNames = ["t1", "P1", "P2"];
opts4.SelectedVariableNames = ["t1", "P1", "P2"];
opts4.VariableTypes = ["double", "double", "double"];

tbl4 = readtable("2_Body.xlsx", opts4, "UseExcel", false);

t_tr = tbl4.t1;
P1_tr = tbl4.P1;
P2_tr = tbl4.P2;

clear opts tbl4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CALCULATE THERMAL CAPACITANCE
Cth1=rho1*cp1*V1;
Cth2=rho2*cp2*V2;
C_T=Cth1+Cth2;
PoC_T=1/C_T;

T0=20;
% Linear Temperature curve
TL=T0+PoC_T*t;
% Deviations
TD1_10=abs(TB1_10-TL);
TD2_10=abs(TB2_10-TL);
TD1_01=abs(TB1_01-TL);
TD2_01=abs(TB2_01-TL);

% Estimate Rth and k via GRG
T_model=@(MC) (MC(1)*(1-exp(-MC(2)*t)));
initial_guess=[0.1 0.01];
Loss=@(T_model,MC,y_true) sum((T_model(MC)-y_true).^2,'all');
%TM1_10
MC_fit=fminsearch(@(MC) Loss(T_model,MC,TD1_10),initial_guess);
R1B1=MC_fit(1) * sign(TB1_10(end)-TL(end)); % For +/- values of Rth: sign(TB - TL)
k1_10=MC_fit(2);

%TM2_10
MC_fit=fminsearch(@(MC) Loss(T_model,MC,TD2_10),initial_guess);
R1B2=MC_fit(1) * sign(TB2_10(end)-TL(end));
k2_10=MC_fit(2);

%TM1_01
MC_fit=fminsearch(@(MC) Loss(T_model,MC,TD1_01),initial_guess);
R2B1=MC_fit(1) * sign(TB1_01(end)-TL(end));
k1_01=MC_fit(2);

%TM2_01
MC_fit=fminsearch(@(MC) Loss(T_model,MC,TD2_01),initial_guess);
R2B2=MC_fit(1) * sign(TB2_01(end)-TL(end));
k2_01=MC_fit(2);


K1 = 0.5*(k1_10+k1_01);
K2 = 0.5*(k2_10+k2_01);

[TLtr,TF1,TF2] = f_TB(P1_tr,P2_tr,t_tr,R1B1,R1B2,R2B1,R2B2,K1,K2,C_T);

figure
plot(t_tr,TF1,'-+',t_tr,TF2,'--')
legend('Body 1','Body 2')
title('Model Temperature')

toc