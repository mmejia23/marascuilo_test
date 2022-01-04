function [dp, dp_SE2, dp_zeta, dp_CI, C] = marascuilo_test(pHpF, OutOfNum)
% Following Marascuilo (1970), calculate 95% CI for each individual subject's d'.
% This function uses the "one signal significance test" (Marascuilo, 1970, pp. 238-239).
% The estimated CI could be used to identify if each d' is statistically different from zero. The idea for this usage was taken from a response in a ResearchGate question: https://www.researchgate.net/post/How-to-determine-whether-a-participants-d-prime-score-differs-significantly-from-what-would-be-expected-by-chance
%
% Arguments:
% pHpF:      is a nx2 matrix where each row has (in values between 0 and 1):
%            [ proportion of hits(pHits), proportion of false alarms (pFA) ]
% OutOfNum:  is a nx2 matrix where each row has:
%            [ number of signal trials (nHits), number of noise trials (nFA) ]
% 
% Output:
% test_dp:     d' for each row.
% test_SE2dp:  SE^2 calculated following Marascuilo (1970).
% zeta_dp:     Z value for each d', Z = d' / SE
% dp_CI:       nx2 matrix where each row has lower and upper 95% CI.
%
% Notes:
% Consider the difference in dp and C calculation with Palamedes toolbox:
% Where R is the ratio between variances of Signal 1 and Signal 2.
% Assuming R=1, we have the simplified formula used in this function.
% If R is not = 1, then:
% k = sqrt(2./(1+R.^2));
% dp = k.*(zHit - R.*zFA);
% C = -R./(1+R).*(zHit+zFA);
% but this is not implemented here.
% 
% Reference:
% Marascuilo, L. A. (1970). Extensions of the significance test for one-parameter signal detection hypotheses. Psychometrika, 35(2), 237–243. https://doi.org/10.1007/BF02291265
% 

pHit = pHpF(:, 1);
pFA = pHpF(:, 2);
nHits = OutOfNum(:, 1);
nFA = OutOfNum(:, 2);

zHit = norminv(pHit);
zFA = norminv(pFA);

dp = zHit - zFA;
dp_SE2 = ( pHit .* (1 - pHit) ) ./ ( nHits .* ( normpdf(zHit) ).^2 ) + ( pFA .* (1 - pFA) ) ./ ( nFA .* ( normpdf(zFA) ).^2 );

dp_zeta = dp ./ sqrt(dp_SE2);
dp_CI = [ dp - 1.96 .* sqrt(dp_SE2), dp + 1.96 .* sqrt(dp_SE2)];

C = -(zHit + zFA)/2;

end