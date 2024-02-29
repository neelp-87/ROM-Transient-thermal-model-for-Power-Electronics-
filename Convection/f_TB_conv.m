function [TF1]=f_TB_conv(P1,t,R1B1,K1,C_T)
P1(1)=0;
%% 
% i - indices where P changes
% N - total number of changes

% Assumption: P can change at every time step
% If P changes, PDelta = P-P_prev; t0 = t_prev;
% else PDelta = PDelta_prev; to = to_prev;
T1C=P1*R1B1;


T1 = zeros(size(P1));
T0 = 20;

idxChange1 = find(diff(T1C)~=0)+1;

% Formula
for ct=1:numel(idxChange1)
    T1CDelta = T1C(idxChange1(ct)) - T1C(idxChange1(ct)-1);
    t0 = t(idxChange1(ct)-1);
    for i=idxChange1(ct):numel(P1)       
        T1(i) = T1(i) + T1CDelta*(1-exp(-K1*(t(i)-t0)));
    end
end


TF1 = T0 + T1;

end
