function [TLf]=f_OB(P1,t,C_T)
P1(1)=0;
%% 
% i - indices where P changes
% N - total number of changes

% Assumption: P can change at every time step
% If P changes, PDelta = P-P_prev; t0 = t_prev;
% else PDelta = PDelta_prev; to = to_prev;

% Constant
PoC=(P1)/C_T;

TL = zeros(size(PoC));
T0 = 20;


idxChange0 = find(diff(PoC)~=0)+1;

% Formula
for ct=1:numel(idxChange0)
    PoCDelta = PoC(idxChange0(ct)) - PoC(idxChange0(ct)-1);

    t0 = t(idxChange0(ct)-1);
    for i=idxChange0(ct):numel(PoC)       
        TL(i) = TL(i) + PoCDelta*(t(i)-t0);
    end
end


TLf = T0 + TL;

end
