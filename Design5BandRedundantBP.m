function [h0,h1,h2,h3,h4,g0,g1,g2,g3,g4] = Design5BandRedundantBP(h0)

% 27-Sep-2011
% Design 5-band sibling frame with a given lowpass filter h0.
% The synthesis and analysis highpass filters are redundant. 
% The bandpass filters are modulated versions of one another.
%

% Generate LP filter h0 using symflat.m

h0 = h0(:);
g0 = flip(h0);

%PR cond 2:
rh0 = conv(altsign(h0),flip(h0));
rts = roots(rh0);

% Find A(z)
ll = length(rts);

l = 1:ll;
l = l(:);
disp('Roots of H0(-z)G0(z)');
disp([l rts]);

rtsh1 = input('Choose zeros of A(z):');
h1 = poly(rtsh1);
h1 = real(h1);
h1 = conv(h1,flip(h1));
h1 = h1(:);

h2 = altsign(h1);

g1 = -h1;
g2 = -h2;

%
% Place an if statement for obtaining h1 and h2 when h0 is of even length
%

rh1g1 = conv(altsign(h1),g1) + conv(altsign(h2),g2);

K = sqrt(abs(rh0(1)/(rh1g1(1))));

h1 = h1*K;
h2 = h2*K;
g1 = g1*K;
g2 = g2*K;

% Find the highpass filters...
r1 = conv(h0,g0)+conv(h1,g1)+conv(h2,g2);

if (r1(1) < 1e-6 && r1(2) < 1e-6)
    r1 = r1(3:end-2);
elseif (r1(1) < 1e-6)
    r1 = r1(2:end-1);
end

r1 = -r1;
ll = 1:length(r1);
ll = median(ll);
r1(ll) = 2 + r1(ll);
rts1 = roots(r1)

rtsh3 = input('Choose zeros of H3(z):');
h3 = poly(rtsh3);
rtsg3 = input('Choose zeros of G3(z):');
g3 = poly(rtsg3);

h3 = h3(:);
g3 = g3(:);

rhg = conv(h3,g3);
K1 = sqrt(abs(r1(ll)/(2*rhg(ll))));

h3 = K1*h3;
g3 = K1*g3;



%%

% for the case of redundant bandpass filters, try the following:
% from pr condition 2, assume H4(z) = H0(-z),
% find H1(z) (bandpass) vis spectral factorization
% find remaining H1(z) and H2(z) (redundant) via spectral factorization
% of pr condition 1, given h0, h4, and h1.


%%
% Example of sibling frames with h0 of odd length:
H = [
    -0.0625   -0.0442   -0.0442    0.1250
         0   -0.1768    0.1768    0.3750
    0.5625    0.0442    0.0442    0.1250
    1.0000    0.3536   -0.3536   -1.2500
    0.5625    0.0442    0.0442    0.1250
         0   -0.1768    0.1768    0.3750
   -0.0625   -0.0442   -0.0442    0.1250
   ];

G = [
   -0.0625    0.0442    0.0442         0
         0    0.1768   -0.1768         0
    0.5625   -0.0442   -0.0442    0.5000
    1.0000   -0.3536    0.3536   -1.0000
    0.5625   -0.0442   -0.0442    0.5000
         0    0.1768   -0.1768         0
   -0.0625    0.0442    0.0442         0
   ];


