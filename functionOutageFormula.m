function OP = functionOutageFormula(x,alpha1,alpha2,beta1,beta2,b1,b2,c,d)
%Implementation of the outage probability formula in Theorem 1 of the following paper: 
%
%Emil Björnson, Michail Matthaiou, Mérouane Debbah, "A New Look at Dual-Hop
%Relaying: Performance Limits with Hardware Impairments" accepted for
%publication in IEEE Transcations on Communications.
%
%Download article: http://arxiv.org/pdf/1311.2634
%
%This is version 1.2 (Last edited: 2014-03-21).
%
%License: This code is licensed under the GPLv2 license. If you in any way
%use this code for research that results in publications, please cite our
%original article listed above.
%
%The input variables follow the notation from Theorem 1.

if d*x>1
    OP = 1;
else
    factor = zeros(size(b2));
    
    for j=0:alpha1-1
        for n=0:alpha2-1
            for k=0:j
                Cfact = b1.^(alpha2-n-1) .* b2.^(j-k) * beta1^((k-n-1-2*j)/2) * beta2^((n-k+1-2*alpha2)/2) / (factorial(k)*factorial(j-k)*factorial(n)*factorial(alpha2-n-1));
                factor = factor + ( x./(1-d*x) )^(alpha2+j) * Cfact .* (b1.*b2+ c*(1-d*x)/x ).^((n+k+1)/2) .* besselk(n-k+1, 2*sqrt( b1.*b2*x.^2./(beta1*beta2*(1-d*x).^2) + c*x/(beta1*beta2*(1-d*x))));
            end
        end
    end
    
    OP = 1 - 2 * exp(-x./(1-d*x) * (b1/beta2 +b2/beta1) ) .* factor;
    
end

