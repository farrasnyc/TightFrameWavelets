function [h0,h1,h2,h3,g0,g1,g2,g3] = DesignBiframeSym(h0,g0)

% 30-Apr-2006
% function to generate a set of four symmetric 
% frame filters
% M=2, N=4 (with redundant highpass filters)
%
% Input: h0 and g0 lowpass filters satisfying
% inequality...
%
% functions called:
% 
% makeh1
% biframe_info
%
% Last modified: 01/13/2009
% Note to self: K is used twice, that's ok. 

h0 = h0(:);
g0 = g0(:);
Lh = length(h0);
Lg = length(g0);
dl = abs(Lh - Lg);
%LL = max(length(h0),length(g0));

if abs(h0(1)) > 1E-5 && abs(g0(1)) > 1E-5
K = abs(h0(1)*g0(1));
elseif abs(h0(1)) == 0 && abs(g0(1)) ~= 0
    K = abs(h0(2)*g0(1));
elseif abs(h0(1)) ~= 0 && abs(g0(1)) == 0
    K = abs(h0(1)*g0(2));
 elseif abs(h0(1)) == 0 && abs(g0(1)) == 0
   K = abs(h0(2)*g0(2));
end
K = sqrt(K);

if Lh > Lg
    g0 = [zeros(1,dl/2); g0; zeros(1,dl/2)];
elseif Lh < Lg
    h0 = [zeros(1,dl/2); h0; zeros(1,dl/2)];
end

%g0 = h0;

r = conv(makeh1(h0),g0);
% if r(1) < 1E-5
%     r = r(2:end-1);
% elseif r(1) < 1E-5 && r(2) < 1E-5
%     r = r(3:end-2);
% end

rts = roots(r);
ll = length(rts);

l = 1:ll;
l = l(:);
disp('Roots of H0(-z)G0(z)');
disp([l rts]);

h1 = input('Choose zeros of h1:');
h1 = poly(h1);
h1 = real(h1);
h1 = makeh1(h1);
h1 = h1(:);
h1 = h1*K;
g1 = input('Choose zeros of g1:');
g1 = poly(g1);
g1 = real(g1);
g1 = g1(:);
%g1 = h1;

g1 = g1*K;

% if h1(1)*g1(1) < 0
%     g1 = -g1;
% end

% Normalize h1 & g1:
dl = abs(length(h1) - length(g1));
if length(h1) > length(g1)
    g1 = [zeros(1,dl/2); g1; zeros(1,dl/2)];
elseif length(h1) < length(g1)
    h1 = [zeros(1,dl/2); h1; zeros(1,dl/2)];
end

% DOUBLE check, is it necessary?
ChkErr = conv(altsign(h0),g0)+conv(altsign(h1),g1);
if sum(abs(ChkErr)/length(ChkErr) > 1E-5)
    g1 = -g1;
end

%h1 = [h1;zeros(LL - length(h1),1)];
%g1 = [0;0;g1;zeros(1,LL - 2 - length(g1))];

r = conv(h0,g0) + conv(h1,g1);
%K0 = max(abs(r));
r = -r;
%r((length(r)+1)/2) = r((length(r)+1)/2) + 2;
rr = find(abs(r) == max(abs(r)));
r(rr) = r(rr) + 2;
K0 = max(abs(r));

if (r(1) < 1E-5) && (r(2) < 1E-5)
    r = r(3:end-2);
elseif (r(1) < 1E-5)
    r = r(2:end-1);
end

rts = roots(r);
l = 1:length(rts);
l = l';
disp('Roots of 2 - H0(z)G0(z) - z^2H1(z)G1(z)');
disp([l rts]);

h2 = input('Choose zeros of h2:');
h2 = poly(h2);
g2 = input('Choose zeros of g2:');
g2 = poly(g2);

h2 = real(h2);
g2 = real(g2);
h2 = h2(:);
g2 = g2(:);
r2 = conv(h2,g2);
kh = sum(altsign(h2));
kg = sum(altsign(g2));

h2 = h2/abs(kh);
g2 = g2/abs(kg);

    
%K = norm(h2)/norm(g2);
%K = sqrt(K);
%h2 = h2*K;
%g2 = g2*K;

% rr = norm(h2)/norm(g2);
% rr = sqrt(rr);
% 
% h2 = h2/rr;
% g2 = g2*rr;
% 
% h2 = h2(:);
% g2 = g2(:);
%h2 = [h2;zeros(LL-length(h2),1)];
%g2 = [0;0;g2;zeros(LL-length(g2)-2,1)];

%h3 = [0;h2(1:end)];
%g3 = [g2(2:end);0];

%Lmax = max([length(h0) length(h2) length(h3) length(g2) length(g3)]);

% h0 = [h0;zeros(Lmax - length(h0),1)];
% h1 = [h1;zeros(Lmax - length(h1),1)];
% h2 = [h2;zeros(Lmax - length(h2),1)];
%h3 = [h3;zeros(Lmax - length(h3),1)];

% g0 = [g0;zeros(Lmax - length(g0),1)];
% g1 = [g1;zeros(Lmax - length(g1),1)];
% g2 = [g2;zeros(Lmax - length(g2),1)];
%g3 = [g3;zeros(Lmax - length(g3),1)];

h3 = h2;
g3 = g2;
%A = [h0 h1 h2 h3];
%B = [g0 g1 g2 g3];

%biframe_info(A,B,2);

