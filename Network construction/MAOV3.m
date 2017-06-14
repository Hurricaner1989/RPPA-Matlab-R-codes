function [MAOV3, E, dfE] = MAOV3 (X, alpha)
% MAOV3 Three-way multivariate analysis of variances test.
% Computes a three-way multivariate analysis of variance for equal sample
% sizes.
% Considering every possible interaction
% y = m + A + B + C + AB + AC + BC + ABC + E

% Data example:
%
% ------------------------------------------------------
%    Factor                           Variable
% ------------------------------------------------------
%  A     B     C       X1      X2
% ------------------------------------------------------
%  1     1     1       
%  1     1     1
%  1     1     1
%  1     1     2
%  1     1     2
%  1     1     2
%  1     1     3
%  1     1     3
%  1     1     3
%  1     2     1
%  1     2     1
%  1     2     1
%  1     2     2
%  1     2     2
%  1     2     2
%  1     2     3
%  1     2     3
%  1     2     3
%  1     3     1
%  1     3     1
%  1     3     1
%  1     3     2
%  1     3     2
%  1     3     2
%  1     3     3
%  1     3     3
%  1     3     3
%  1     4     1
%  1     4     1
%  1     4     1
%  1     4     2
%  1     4     2
%  1     4     2
%  1     4     3
%  1     4     3
%  1     4     3
%  1     5     1
%  1     5     1
%  1     5     1
%  1     5     2
%  1     5     2
%  1     5     2
%  1     5     3
%  1     5     3
%  1     5     3
%  2     1     1       
%  2     1     1
%  2     1     1
%  2     1     2
%  2     1     2
%  2     1     2
%  2     1     3
%  2     1     3
%  2     1     3
%  2     2     1
%  2     2     1
%  2     2     1
%  2     2     2
%  2     2     2
%  2     2     2
%  2     2     3
%  2     2     3
%  2     2     3
%  2     3     1
%  2     3     1
%  2     3     1
%  2     3     2
%  2     3     2
%  2     3     2
%  2     3     3
%  2     3     3
%  2     3     3
%  2     4     1
%  2     4     1
%  2     4     1
%  2     4     2
%  2     4     2
%  2     4     2
%  2     4     3
%  2     4     3
%  2     4     3
%  2     5     1
%  2     5     1
%  2     5     1
%  2     5     2
%  2     5     2
%  2     5     2
%  2     5     3
%  2     5     3
%  2     5     3
% -----------------------------------------------------------

% Result
% -----------------------------------------------------------------------
% Test        df            Statistic       F-stat   p-value   Conclusion
% -----------------------------------------------------------------------
% A           a-1       det(E)/det(E+HA)
% B           b-1       det(E)/det(E+HB)
% C           c-1       det(E)/det(E+HC)
% A*B      (a-1)(b-1)   det(E)/det(E+HAB)
% A*C      (a-1)(c-1)   det(E)/det(E+HAC)
% B*C      (b-1)(c-1)   det(E)/det(E+HBC)
% A*B*C (a-1)(b-1)(c-1) det(E)/det(E+HABC)
% Error  abc(d-1) 
% Total    n-1
% -----------------------------------------------------------------------
% n=abcd

if nargin < 2
    alpha = 0.05;
end

Y = X;

a = max(X(:,1));          % # of levels in factor A
b = max(X(:,2));          % # of levels in factor B
c = max(X(:,3));          % # of levels in factor C

% r = size(X,1);

X = X(:,4:size(X,2));
[n, p] = size(X);     % n = total # of samples
                      % p = # of variables

dfA = a-1;
dfB = b-1;
dfC = c-1;
dfAB = (a-1)*(b-1);
dfAC = (a-1)*(c-1);
dfBC = (b-1)*(c-1);
dfABC = (a-1)*(b-1)*(c-1);
dfT = n-1;
dfE = dfT-dfA-dfB-dfC-dfAB-dfAC-dfBC-dfABC;

if p > dfE;
   fprintf('Warning: likelihood ratio procedures require the number of independent variables: %2i\n', p);
   fprintf('must be <= error degrees of freedom: %2i\n', dfE);
   disp('So that error sum of squares matrix will be positive definite.');
   return;
end

ym = mean(X, 1);
dT = (X-repmat(ym,n,1))'*(X-repmat(ym,n,1));

ymA = [];
ymB = [];
ymC = [];
ymAB = [];
ymAC = [];
ymBC = [];
ymABC = [];
for i = 1:a
    ymA(i,:) = mean(X(find(Y(:,1)==i),:),1);
    for j = 1:b
        ymB(j,:) = mean(X(find(Y(:,2)==j),:),1);
        ymAB(i,j,:) = mean(X(find(Y(:,1)==i&Y(:,2)==j),:),1);
        for k = 1:c
            ymC(k,:) = mean(X(find(Y(:,3)==k),:),1);
            ymAC(i,k,:) = mean(X(find(Y(:,1)==i&Y(:,3)==k),:),1);
            ymBC(j,k,:) = mean(X(find(Y(:,2)==j&Y(:,3)==k),:),1);
            ymABC(i,j,k,:) = mean(X(find(Y(:,1)==i&Y(:,2)==j&Y(:,3)==k),:),1);
        end
    end
end

HA = zeros(p,p);
HB = zeros(p,p);
HC = zeros(p,p);
HAB = zeros(p,p);
HAC = zeros(p,p);
HBC = zeros(p,p);
HABC = zeros(p,p);
E = zeros(p,p);
for i = 1:a
    for j = 1:b
        for k = 1:c
            f = find(Y(:,1)==i&Y(:,2)==j&Y(:,3)==k);
            z = X(f,:);
            for l = 1:numel(f)
                HA = HA+(ymA(i,:)-ym)'*(ymA(i,:)-ym);
                HB = HB+(ymB(j,:)-ym)'*(ymB(j,:)-ym);
                HC = HC+(ymC(k,:)-ym)'*(ymC(k,:)-ym);
                HAB = HAB+(reshape(ymAB(i,j,:),1,p)-ymA(i,:)-ymB(j,:)+ym)'*(reshape(ymAB(i,j,:),1,p)-ymA(i,:)-ymB(j,:)+ym);
                HAC = HAC+(reshape(ymAC(i,k,:),1,p)-ymA(i,:)-ymC(k,:)+ym)'*(reshape(ymAC(i,k,:),1,p)-ymA(i,:)-ymC(k,:)+ym);
                HBC = HBC+(reshape(ymBC(j,k,:),1,p)-ymB(j,:)-ymC(k,:)+ym)'*(reshape(ymBC(j,k,:),1,p)-ymB(j,:)-ymC(k,:)+ym);
                HABC = HABC+(reshape(ymABC(i,j,k,:),1,p)-reshape(ymAB(i,j,:),1,p)-reshape(ymAC(i,k,:),1,p)-reshape(ymBC(j,k,:),1,p)...
                    +ymA(i,:)+ymB(j,:)+ymC(k,:)-ym)'*(reshape(ymABC(i,j,k,:),1,p)...
                    -reshape(ymAB(i,j,:),1,p)-reshape(ymAC(i,k,:),1,p)-reshape(ymBC(j,k,:),1,p)...
                    +ymA(i,:)+ymB(j,:)+ymC(k,:)-ym);
                E = E+(z(l,:)-reshape(ymABC(i,j,k,:),1,p))'*(z(l,:)-reshape(ymABC(i,j,k,:),1,p));
            end
        end
    end
end

LA = det(E)/det(E+HA);
LB = det(E)/det(E+HB);
LC = det(E)/det(E+HC);
LAB = det(E)/det(E+HAB);
LAC = det(E)/det(E+HAC);
LBC = det(E)/det(E+HBC);
LABC = det(E)/det(E+HABC);

if p == 2
    FA = ((1-sqrt(LA))/sqrt(LA))*((dfE-1)/dfA);
    FB = ((1-sqrt(LB))/sqrt(LB))*((dfE-1)/dfB);
    FC = ((1-sqrt(LC))/sqrt(LC))*((dfE-1)/dfC);
    FAB = ((1-sqrt(LAB))/sqrt(LAB))*((dfE-1)/dfAB);
    FAC = ((1-sqrt(LAC))/sqrt(LAC))*((dfE-1)/dfAC);
    FBC = ((1-sqrt(LBC))/sqrt(LBC))*((dfE-1)/dfBC);
    FABC = ((1-sqrt(LABC))/sqrt(LABC))*((dfE-1)/dfABC);
    FdfAE = 2*(dfE-1);
    FdfBE = 2*(dfE-1);
    FdfCE = 2*(dfE-1);
    FdfABE = 2*(dfE-1);
    FdfACE = 2*(dfE-1);
    FdfBCE = 2*(dfE-1);
    FdfABCE = 2*(dfE-1);
    FdfA = 2*dfA;
    FdfB = 2*dfB;
    FdfC = 2*dfC;
    FdfAB = 2*dfAB;
    FdfAC = 2*dfAC;
    FdfBC = 2*dfBC;
    FdfABC = 2*dfABC;
elseif p == 1
    FA = ((1-LA)/LA)*(dfE/dfA);
    FB = ((1-LB)/LB)*(dfE/dfB);
    FC = ((1-LC)/LC)*(dfE/dfC);
    FAB = ((1-LAB)/LAB)*(dfE/dfAB);
    FAC = ((1-LAC)/LAC)*(dfE/dfAC);
    FBC = ((1-LBC)/LBC)*(dfE/dfBC);
    FABC = ((1-LABC)/LABC)*(dfE/dfABC);
    FdfAE = dfE;
    FdfBE = dfE;
    FdfCE = dfE;
    FdfABE = dfE;
    FdfACE = dfE;
    FdfBCE = dfE;
    FdfABCE = dfE;
    FdfA = dfA;
    FdfB = dfB;
    FdfC = dfC;
    FdfAB = dfAB;
    FdfAC = dfAC;
    FdfBC = dfBC;
    FdfABC = dfABC;
else
    if dfA == 1
        FA = ((1-LA)/LA)*((dfE-p+1)/p);
        FdfA = p;
        FdfAE = dfE-p+1;
    elseif dfA == 2
        FA = ((1-sqrt(LA))/sqrt(LA))*((dfE-p+1)/p);
        FdfA = 2*p;
        FdfAE = 2*(dfE-p+1);
    end
    if dfB == 1
        FB = ((1-LB)/LB)*((dfE-p+1)/p);
        FdfB = p;
        FdfBE = dfE-p+1;
    elseif dfB == 2
        FB = ((1-sqrt(LB))/sqrt(LB))*((dfE-p+1)/p);
        FdfB = 2*p;
        FdfBE = 2*(dfE-p+1);
    end
    if dfC == 1
        FC = ((1-LC)/LC)*((dfE-p+1)/p);
        FdfC = p;
        FdfCE = dfE-p+1;
    elseif dfC == 2
        FC = ((1-sqrt(LC))/sqrt(LC))*((dfE-p+1)/p);
        FdfC = 2*p;
        FdfCE = 2*(dfE-p+1);
    end
    if dfAB == 1
        FAB = ((1-LAB)/LAB)*((dfE-p+1)/p);
        FdfAB = p;
        FdfABE = dfE-p+1;
    elseif dfAB == 2
        FAB = ((1-sqrt(LAB))/sqrt(LAB))*((dfE-p+1)/p);
        FdfAB = 2*p;
        FdfABE = 2*(dfE-p+1);
    end
    if dfAC == 1
        FAC = ((1-LAC)/LAC)*((dfE-p+1)/p);
        FdfAC = p;
        FdfACE = dfE-p+1;
    elseif dfAC == 2
        FAC = ((1-sqrt(LAC))/sqrt(LAC))*((dfE-p+1)/p);
        FdfAC = 2*p;
        FdfACE = 2*(dfE-p+1);
    end
    if dfBC == 1
        FBC = ((1-LBC)/LBC)*((dfE-p+1)/p);
        FdfBC = p;
        FdfBCE = dfE-p+1;
    elseif dfBC == 2
        FBC = ((1-sqrt(LBC))/sqrt(LBC))*((dfE-p+1)/p);
        FdfBC = 2*p;
        FdfBCE = 2*(dfE-p+1);
    end
    if dfABC == 1
        FABC = ((1-LABC)/LABC)*((dfE-p+1)/p);
        FdfABC = p;
        FdfABCE = dfE-p+1;
    elseif dfABC == 2
        FABC = ((1-sqrt(LABC))/sqrt(LABC))*((dfE-p+1)/p);
        FdfABC = 2*p;
        FdfABCE = 2*(dfE-p+1);
    end
end

pvA = 1-fcdf(FA, FdfA, FdfAE);
pvB = 1-fcdf(FB, FdfB, FdfBE);
pvC = 1-fcdf(FC, FdfC, FdfCE);
pvAB = 1-fcdf(FAB, FdfAB, FdfABE);
pvAC = 1-fcdf(FAC, FdfAC, FdfACE);
pvBC = 1-fcdf(FBC, FdfBC, FdfBCE);
pvABC = 1-fcdf(FABC, FdfABC, FdfABCE);
    
if pvA < alpha
    CA = 'S';
else
    CA = 'NS';
end
if pvB < alpha
    CB = 'S';
else
    CB = 'NS';
end
if pvC < alpha
    CC = 'S';
else
    CC = 'NS';
end
if pvAB < alpha
    CAB = 'S';
else
    CAB = 'NS';
end
if pvAC < alpha
    CAC = 'S';
else
    CAC = 'NS';
end
if pvBC < alpha
    CBC = 'S';
else
    CBC = 'NS';
end
if pvABC < alpha
    CABC = 'S';
else
    CABC = 'NS';
end

MAOV3 = {'Test'; 'A'; 'B'; 'C'; 'A*B'; 'A*C'; 'B*C'; 'A*B*C'};
MAOV3 = [MAOV3 {'Statistic'; FA; FB; FC; FAB; FAC; FBC; FABC}];
MAOV3 = [MAOV3 {'df1'; FdfA; FdfB; FdfC; FdfAB; FdfAC; FdfBC; FdfABC}];
MAOV3 = [MAOV3 {'df2'; FdfAE; FdfBE; FdfCE; FdfABE; FdfACE; FdfBCE; FdfABCE}];
MAOV3 = [MAOV3 {'F'; FA; FB; FC; FAB; FAC; FBC; FABC}];
MAOV3 = [MAOV3 {'p-value'; pvA; pvB; pvC; pvAB; pvAC; pvBC; pvABC}];
MAOV3 = [MAOV3 {'Conclusion'; CA; CB; CC; CAB; CAC; CBC; CABC}];

end