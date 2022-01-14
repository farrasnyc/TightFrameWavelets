function W = WavMatM4N8(H, N, k0, shift)
% WavMat -- Transformation Matrix of FWT_PO
%  Usage
%    W = WavMat(h, N, k0, shift)
%  Inputs
%    H      Matrix containing tight frame filters with M=4
%    N      size of matrix/length of data. Should be multiple of power of 4.
%      
%    k0     depth of transformation. Ranges from 1 to J=log4(N).
%           Default is J. 
%    shift  the matrix is not unique an any integer shift gives
%           a valid transformation. Default is 2.
%  Outputs
%    W      N x L transformation matrix  
%
%  Description
%    For a quadrature mirror filter h (low pass) the wavelet
%    matrix is formed. The algorithm is described in 
%    [BV99] Vidakovic, B. (1999). Statistical Modeling By Wavelets, Wiley,
%    on pages 115-116.
%    Any shift is valid. 
%
%  Usage
%    We will mimic the example 4.3.1 from [BV99] page 112.
%   > dat=[1 0 -3 2 1 0 1 2];
%   > W = WavMat(MakeONFilter('Haar',99),2^3,3,2);
%   > wt = W * dat' %should be [sqrt(2)  |  -sqrt(2) |   1 -1  | ...         
%              %  1/sqrt(2) -5/sqrt(2) 1/sqrt(2) - 1/sqrt(2) ]
%   > data = W' * wt % should return you to the 'dat'
%
%  See Also
%    FWT_PO, IWT_PO, MakeONFilter
%
if nargin==3
    shift = 4;
end
%J = log2(N);
M = 4;
J = k0;
L = N/M^J; % Data must of length multiple of M^J

%2^J*ceil(N/2^J)-N;
if (L ~= floor(L) )
    error('Data must be of length L*4^J, where L is an integer and J is the number of filterbank stages used.')
end

if nargin==2
    shift = 4; % shift = 4?
    k0 = J;
end

h0 = H(1,:); 
h1 = H(2,:); 
h2 = H(3,:); 
h3 = H(4,:); 
h4 = H(5,:); 
h5 = H(6,:); 
h6 = H(7,:); 
h7 = H(8,:); 


h0=[h0,zeros(1,N)]; %extend filters by 0's to sample by modulus
h1=[h1,zeros(1,N)]; 
h2=[h2,zeros(1,N)];
h3=[h3,zeros(1,N)];
h4=[h4,zeros(1,N)];
h5=[h5,zeros(1,N)];
h6=[h6,zeros(1,N)];
h7=[h7,zeros(1,N)];

oldmat = speye(M^(J-k0)*L); 
for k= k0:-1:1
    %clear gmat; clear hmat;
    clear h0mat h1mat h2mat h3mat h4mat h5mat h6mat h7mat
         ubJk = M^(J-k)*L; ubJk1 = M^(J-k+1)*L;
   for  jj= 1:ubJk
       for ii=1:ubJk1
           modulus = mod(N+ii-M*jj+shift,ubJk1);
           modulus = modulus + (modulus == 0)*ubJk1;
           h0mat(ii,jj) = h0(modulus);
           h1mat(ii,jj) = h1(modulus);
           h2mat(ii,jj) = h2(modulus);
           h3mat(ii,jj) = h3(modulus);
           h4mat(ii,jj) = h4(modulus);
           h5mat(ii,jj) = h5(modulus);
           h6mat(ii,jj) = h6(modulus);
           h7mat(ii,jj) = h7(modulus);
       end
   end
  W = [oldmat * h0mat'; h1mat' ; h2mat' ; h3mat' ; h4mat' ; h5mat' ; h6mat' ; h7mat'];
   oldmat = sparse(W);
end

W = sparse(W);
%
% 
% Copyright (c) 2004. Brani Vidakovic
%        
%  
% ver 1.0 Built 8/24/04; ver 1.2 Built 12/1/2004
% This is Copyrighted Material
% Comments? e-mail brani@isye.gatech.edu
%   
