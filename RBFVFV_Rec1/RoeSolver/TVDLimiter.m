function phi = TVDLimiter(dright,dleft)
% Van Leer
negative = dright.*dleft <= 0;
phi = 2*dright./(dleft + dright);
phi(negative) = 0;

% SB
% r = dright./dleft;
% phi = max(max(zeros(size(dright)),min(2*r,1)),min(r,2));

% MM
% r = dright./dleft;
% phi = max(zeros(size(dright)),min(r,1));

end