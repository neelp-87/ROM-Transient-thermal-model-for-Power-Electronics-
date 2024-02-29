tic
clear;
clc;

% Material - Ag 99.9%
rho=10524;
cp=236;
k=411;
Lx=5e-3; 
Ly=5e-3;
Lz=5e-3;

V=Lx*Ly*Lz;

C_th=rho*cp*V;

%%%%%%%%%%%%
opts = spreadsheetImportOptions("NumVariables", 2);
opts.Sheet = "IP";
opts.DataRange = "A2:B27";
opts.VariableNames = ["t1", "P"];
opts.SelectedVariableNames = ["t1", "P"];
opts.VariableTypes = ["double", "double"];

tbl = readtable("1_Body.xlsx", opts, "UseExcel", false);
t1 = tbl.t1;
P = tbl.P;

clear opts tbl
%%%%%%%%%%%%

Tf=f_OB(P,t1,C_th);

plot(t1,Tf)
title('Model temperature')

toc