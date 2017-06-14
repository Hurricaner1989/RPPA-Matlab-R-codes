function [T] = MAOV_ttest (X, E, dfE, alpha)
% MAOV_ttest multivariate analysis of variance based t-test.
% X: data (two-group data)
% E: MANOVA error sum of squares
% dfE: MANOVA error degree of freedom
% alpha: significance level

% Data example:
%
% ------------------------------------------------------------
%         Factor                          Variable
% ------------------------------------------------------------
%           A                         X1     X2   ...
% ------------------------------------------------------------
%           1
%           1
%           .
%           .
%           .
%           1
%           2
%           .
%           .
%           .
%           2
% ------------------------------------------------------------

if nargin < 4
    alpha = 0.05;
end

Y = X;
group = X(:,1);
X = X(:,2:end);
p = size(X,2);
n = [];
ym = [];
for i = 1:2
    n(i) = numel(find(group==i));
    ym(i,:) = mean(X(find(group==i),:),1);
end

stat = (1/((1/n(1)+1/n(2))/dfE))*(ym(1,:)-ym(2,:))*inv(E)*(ym(1,:)-ym(2,:))';
F = (dfE-p+1)*stat/(dfE*p);
df1 = p;
df2 = dfE-p+1;
pv = 1-fcdf(F, df1, df2);
if pv <alpha
    C = 'S';
else
    C = 'NS';
end

T = {'Statistic' 'F' 'df1' 'df2' 'p-value' 'Conclusion'};
T = [T; {stat F df1 df2 pv C}];

end