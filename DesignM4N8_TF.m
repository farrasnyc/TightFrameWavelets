function [h0,h1,h2,h3,h4,h5,h6,h7,A,B] = DesignM4N8_TF(h0) 

% 26-Dec-2005
% function to generate a set of eight symmetric tight 
% frame filters starting with a lowpass filter h0.
% M=4, N=8
%
% Last modified: 09/14/2010

%L = length(h0);
if mod(length(h0),4) == 2
    h0 = [0;h0;0];
end
L = length(h0);

% Find H0 polyphase components for M=4 then
% generate h1, h2 & h3

h00 = h0(1:4:end);
h01 = h0(2:4:end);
h02 = h0(3:4:end);
h03 = h0(4:4:end);

h00 = up(h00,4);
h01 = up(h01,4);
h01 = [0; h01(1:end-1)];
h02 = up(h02,4);
h02 = [0; 0; h02(1:end-2)];
h03 = up(h03,4);
h03 = [0; 0; 0; h03(1:end-3)];

h1 = h00 - h01 + h02 - h03;
h2 = h00 + h01 - h02 - h03;
h3 = -h00 + h01 + h02 - h03;

% Find A(z)A(1/z) & B(z)B(1/z)

kk = j.^(1:L);
ll = (-j).^(1:L);

r00 = conv(h0 + h0.*kk' + h0.*ll',flip(h0));
r01 = conv(h1 + h1.*kk' + h1.*ll',flip(h1));
r02 = conv(h2 + h2.*kk' + h2.*ll',flip(h2));
r03 = conv(h3 + h3.*kk' + h3.*ll',flip(h3));

r10 = conv(h0 - h0.*kk' - h0.*ll',flip(h0));
r11 = conv(h1 - h1.*kk' - h1.*ll',flip(h1));
r12 = conv(h2 - h2.*kk' - h2.*ll',flip(h2));
r13 = conv(h3 - h3.*kk' - h3.*ll',flip(h3));

r0 = r00 + r01 + r02 + r03;
r1 = r10 + r11 + r12 + r13;

r0 = -r0;
r1 = -r1;

r0((length(r0) + 1)/2) = 4 + r0((length(r0) + 1)/2);
r1((length(r1) + 1)/2) = 4 + r1((length(r1) + 1)/2);

r0 = r0/16;
r1 = r1/16;

% "downsample" r0 and r1 to simplify the resulting zeros.

r0 = r0(4:4:end);
r1 = r1(4:4:end);

rts0 = roots(r0);
rts1 = roots(r1);

% generate filter h4
disp('r0 zeros (rts0)');
disp(roots(r0))
z0 = input('Choose zeros of A(z): ')
rootsA= [z0];

disp('r1 zeros (rts1)');
disp(roots(r1))
z1 = input('Choose zeros of B(z): ')
rootsB = [z1];

A = poly(rootsA);
B = poly(rootsB);

% A & B normalization

AA = conv(A,flip(A));
BB = conv(B,flip(B));

Ka = sqrt(r0((length(r0) + 1)/2)/AA((length(AA) + 1)/2));
Kb = sqrt(r1((length(r1) + 1)/2)/BB((length(BB) + 1)/2));

A = A*Ka;
B = B*Kb;

A = up(A,4);
B = up(B,4);
B = [0 B(1:end-1)];

h4 = A + B + flip(B) + flip(A);

h40 = up(h4(1:4:end),4);
h41 = up(h4(2:4:end),4);
h41 = [0 h41(1:end-1)];
h42 = up(h4(3:4:end),4);
h42 = [0 0 h42(1:end-2)];
h43 = up(h4(4:4:end),4);
h43 = [0 0 0 h43(1:end-3)];

% generate remaining filters

h5 = -h40 + h41 + h42 - h43;
h6 = -h40 + h41 - h42 + h43;
h7 = -h40 - h41 + h42 + h43;

h4 = h4(:);
h5 = h5(:);
h6 = h6(:);
h7 = h7(:);

AA = [h0 h1 h2 h3 h4 h5 h6 h7];
plot8bandM4(AA)

function y = up(x,M)
% y = up(x,M); 
% Upsample x by M
%
% Input:  x    - input signal (vector!)
%         M    - upsampling factor
%
% Output: y - output signal

s = size(x);
if s(1) > s(2)
   y = zeros(M*s(1),s(2));
else
   y = zeros(s(1),M*s(2));
end
ylen = length(y);
y(1:M:ylen) = x;
%y = y(1:length(y)-1);

y = y(1:length(y));
