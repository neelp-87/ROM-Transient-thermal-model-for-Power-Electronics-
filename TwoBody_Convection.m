tic
clear;
clc;

% IMPORT CHARACTERISTIC TEMPERATURES for P1=1
opts1 = spreadsheetImportOptions("NumVariables", 2);
opts1.Sheet = "T1";
opts1.DataRange = "A2:B26";
opts1.VariableNames = ["t", "TB1_10"];
opts1.SelectedVariableNames = ["t", "TB1_10"];
opts1.VariableTypes = ["double", "double"];
tbl1 = readtable("2_Body_Convection.xlsx", opts1, "UseExcel", false);

t = tbl1.t;
TB1 = tbl1.TB1_10;

% IMPORT THERMAL PROPERTIES
opts3 = spreadsheetImportOptions("NumVariables", 8);

opts3.Sheet = "properties";
opts3.DataRange = "A2:H2";

opts3.VariableNames = ["V1", "rho1", "cp1", "k1", "V2", "rho2", "cp2", "k2"];
opts3.SelectedVariableNames = ["V1", "rho1", "cp1", "k1", "V2", "rho2", "cp2", "k2"];
opts3.VariableTypes = ["double", "double", "double", "double", "double", "double", "double", "double"];

tbl3 = readtable("2_Body_Convection.xlsx", opts3, "UseExcel", false);

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
clear opts3 tbl3

% IMPORT TRANSIENT POWER INPUTS 
opts4 = spreadsheetImportOptions("NumVariables", 2);
opts4.Sheet = "Inputs";
opts4.DataRange = "A2:B27";
opts4.VariableNames = ["t1", "P1"];
opts4.SelectedVariableNames = ["t1", "P1"];
opts4.VariableTypes = ["double", "double"];

tbl4 = readtable("2_Body_Convection.xlsx", opts4, "UseExcel", false);

t_tr = tbl4.t1;
P1_tr = tbl4.P1;

clear opts tbl4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CALCULATE THERMAL CAPACITANCE
Cth1=rho1*cp1*V1;
Cth2=rho2*cp2*V2;
C_T=Cth1+Cth2;
PoC_T=1/C_T;

T0=20;

% Estimate Rth and k via GRG
T_model=@(MC) T0+(MC(1)*(1-exp(-MC(2)*t)));
initial_guess=[100 1];

Loss=@(T_model,MC,y_true) sum((T_model(MC)-y_true).^2,'all');
MC_fit=fminsearch(@(MC) Loss(T_model,MC,TB1),initial_guess);
R1B1=MC_fit(1);
K1=MC_fit(2);
TM = T0 + R1B1*(1-exp(-K1*t));


TF1 = f_TB_conv(P1_tr,t_tr,R1B1,K1,C_T);

figure
plot(t_tr,TF1)
title('Model Temperature')

toc