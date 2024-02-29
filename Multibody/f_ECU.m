function [TLf,Tfinal1,Tfinal2,Tfinal3]=f_ECU(P1,P2,t,R1B1,R1B2,R1B3,R2B1,R2B2,R2B3,K1,K2,K3,C_T)
P1(1)=0;
P2(1)=0;
%% 
% i - indices where P changes
% N - total number of changes

% Assumption: P can change at every time step
% If P changes, PDelta = P-P_prev; t0 = t_prev;
% else PDelta = PDelta_prev; to = to_prev;
%% Constants
PoC=(P1+P2)/C_T;

T1C=P1*R1B1+P2*R2B1;
T2C=P1*R1B2+P2*R2B2;
T3C=P1*R1B3+P2*R2B3;

T1 = zeros(size(P1));
T2 = zeros(size(P2));
T3 = zeros(size(P2));

TL = zeros(size(PoC));
T0 = 20;

idxChange0 = find(diff(PoC)~=0)+1;
%% Compute T=f(P)
for ct=1:numel(idxChange0)
    PoCDelta = PoC(idxChange0(ct)) - PoC(idxChange0(ct)-1);
    T1CDelta = T1C(idxChange0(ct)) - T1C(idxChange0(ct)-1);
    T2CDelta = T2C(idxChange0(ct)) - T2C(idxChange0(ct)-1);
    T3CDelta = T3C(idxChange0(ct)) - T3C(idxChange0(ct)-1);

    t0 = t(idxChange0(ct)-1);
    
    for i=idxChange0(ct):numel(PoC)       
        TL(i) = TL(i) + PoCDelta*(t(i)-t0);
        T1(i) = T1(i) + T1CDelta*(1-exp(-K1*(t(i)-t0)));
        T2(i) = T2(i) + T2CDelta*(1-exp(-K2*(t(i)-t0)));
        T3(i) = T3(i) + T3CDelta*(1-exp(-K3*(t(i)-t0)));
    end
end


TLf = T0 + TL;
Tfinal1 = T0 + TL + T1;
Tfinal2 = T0 + TL + T2;
Tfinal3 = T0 + TL + T3;
end
