%This Matlab script can be used to generate Figure 6 of the paper:
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
%Please note that the channels are generated randomly, thus the results
%will not be exactly the same as in the paper.

clear all
close all


%%Simulation parameters

rng('shuffle'); %Initiate the random number generators with a random seed
%%If rng('shuffle'); is not supported by your Matlab version, you can use
%%the following commands instead:
%randn('state',sum(100*clock));


%Select if the upper bound in Theorem 3 should be computed. This is based
%on symbolic toolbox in Matlab, takes a long time to run, and is only
%supported by newer version of Matlab.
computeUpperbound = true;


% Relaying with transceiver hardware impairments

%Number of channel iterations. This must be a very large number to make the
%Monte Carlo simulations fit the analytic curves for small BERs, e.g. 1e-6
nbrOfIterations = 1000000;


%Select fading parameters with Nakagami-m fading
alpha1 = 2; %alpha fading parameter at first hop
beta1 = 1;  %beta fading parameter at first hop

alpha2 = 2; %alpha fading parameter at second hop
beta2 = 1+1e-6;  %beta fading parameter at second hop (with small perturbation for numerical reasons)


%Generate fading realizations. RHO is the squared norm.
RHO = [gamrnd(alpha1,beta1,1,nbrOfIterations); gamrnd(alpha2,beta2,1,nbrOfIterations)];


prelog = 1/2;  %The transmission takes 2 channel uses
R = 5;   %Transmit information rate (per 2 channel uses)
x = 2.^R-1;  %Corresponding outage thresholds

N1 = 1;      %Normalized noise variance at relay
N2 = 1;      %Normalized noise variance at destination

SNR_dB = 0:1:50; %SNR range at first hop
P1_dB = SNR_dB - 10*log10(alpha1*beta1); %Corresponding transmit power range at source (in dB)
P1 = 10.^(P1_dB/10); %Corresponding transmit power range at source

P2 = P1; %Set same transmit power range at relay.

kappaRange = [0.05 0.15];


%Placeholders for Monte-Carlo results with Amplify-and-Forward (AF) with variable gain
rates_ideal_af_v = zeros(length(P1),1); %Ergodic achievable rates, AF variable gain, ideal hardware
rates_nonideal_af_v = zeros(length(P1),length(kappaRange)); %Ergodic achievable rates, AF variable gain, non-ideal hardware

rates_approximation36 = zeros(length(P1),length(kappaRange)); %Ergodic achievable rates, approximation from Eq. (36)
rates_upperbound33 = zeros(length(P1),length(kappaRange)); %Ergodic achievable rates, upper bound from Eq. (33) in Theorem 3


%Define symbolic variable, if the upper bound should be computed
if computeUpperbound==true
    syms t;
end


%Go through all cases of transmit power at the source and assume that the
%power on the relay is varied in the same way. Compute achievable rates
%by Monte-Carlo simulations and using the approximations in Eq. (36) and
%upper bound in Eq. (33) of the paper.
for p = 1:length(P1)
    
    for kk = 1:length(kappaRange)
        
        kappa1 = kappaRange(kk); %Error Vector Magnitude (EVM) at first hop
        kappa2 = kappaRange(kk); %Error Vector Magnitude (EVM) at second hop
        d = kappa1^2+kappa2^2+kappa1^2*kappa2^2; %Recurring function of kappa1 and kappa2
            
        %Monte-Carlo simulations
        for ind = 1:nbrOfIterations
            
            %Extract channel gain for current realization
            rho1 = RHO(1,ind);
            rho2 = RHO(2,ind);

            %Compute achievable rates with ideal hardware only once (since
            %independent of kappa1 and kappa2).
            if kk==1
                SNDR_ideal_af_v = (rho1*rho2)/(rho2*N1/P1(p) + N2*(rho1/P2(p)+N1/P2(p)/P1(p))); %Compute signal-to-noise-and-distortion ratios (SNDRs) with ideal hardware
                rates_ideal_af_v(p,kk) = rates_ideal_af_v(p,kk) + prelog*log2(1+SNDR_ideal_af_v)/nbrOfIterations; %Averaging to get ergodic achievable rates
            end
            
            %Compute achievable rates with non-ideal hardware
            SNDR_nonideal_af_v = (rho1*rho2)/(rho1*rho2*d + rho1*N2*(1+kappa1^2)/P2(p) + rho2*N1*(1+kappa2^2)/P1(p) + N1*N2/P1(p)/P2(p));  %Compute signal-to-noise-and-distortion ratios (SNDRs) with non-ideal hardware
            rates_nonideal_af_v(p,kk) = rates_nonideal_af_v(p,kk) + prelog*log2(1+SNDR_nonideal_af_v)/nbrOfIterations; %Averaging to get ergodic achievable rates
            
            
        end
        
        %Hardware impairment parameters defined in the paper
        b1 = (1+kappa1^2)*N2./P2(p);
        b2 = (1+kappa2^2)*N1./P1(p);
        c = N1*N2./(P1(p).*P2(p));
        a = b1/b2;
        b = c/b2;
        
        %Compute achievable rates according to the approximation in Eq. (33)
        rates_approximation36(p,kk) = prelog*log2(1 + alpha1*alpha2*beta1*beta2 / (alpha1*alpha2*beta1*beta2*d + alpha1*beta1*b1 + alpha2*beta2*b2 + c));
        
        
        %Compute upper bound on the achievable rates using the formula in
        %Theorem 3. This part utilizes the symbolic toolbox in Matlab and
        %is only supported new versions of Matlab. This part can be turned
        %off by setting computeUpperbound=false.
        if computeUpperbound == true
            J = 0;
            for n = 0:alpha1-1
                for k = 0:alpha2-1
                    for m = 0:k
        
                        for q = 0:n+m+2
                            
                            %Define the function inside the curly brackets in Eq. (34)
                            funct = symfun(exp(c*t/(2*b1)) *whittakerW(-(n+m+2)/2, (n-m+1)/2, c/(2*b1) * (t - sqrt(t^2-4*b1/beta1/beta2/b2)) ) * whittakerW(-(n+m+2)/2, (n-m+1)/2, c/(2*b1) * (t + sqrt(t^2-4*b1/beta1/beta2/b2)) ), t);
        
                            derivative = diff( funct, alpha1+k-q+1); %Compute the derivative of the expressions in the curly brackets.

                            %Compute a term of the summation in Eq. (34)
                            J = J + nchoosek(n+m+2,q) * ((n+1)*beta1^((n-m+2-2*alpha1)/2)*beta2^((m-n-2*k)/2))/(factorial(k-m)*factorial(alpha1-n-1)) *  ( (b1/b2)^((n-m+2+2*k)/2) * (c/b1)^q ) / ( c * (-1)^(alpha1+k-q+1) ) *double(derivative(1/beta1+b1/beta2/b2));

                        end
                    end
                end
            end
            rates_upperbound33(p,kk) = prelog*log2(1+ J/(J*d+1));
        end
        
    end
    
end


%Plot the results
figure; hold on; box on;

plot(SNR_dB,rates_ideal_af_v(:,1),'k','LineWidth',1);

for kk = 1:length(kappaRange)
        kappa1 = kappaRange(kk); %Error Vector Magnitude (EVM) at first hop
        kappa2 = kappaRange(kk); %Error Vector Magnitude (EVM) at second hop
        d = kappa1^2+kappa2^2+kappa1^2*kappa2^2; %Recurring function of kappa1 and kappa2
    
    
        plot(SNR_dB,rates_approximation36(:,kk),'b--','LineWidth',1);
        plot(SNR_dB,rates_upperbound33(:,kk),'r-.','LineWidth',1);
        plot([min(SNR_dB) max(SNR_dB)],0.5*log2(1+1/d)*ones(1,2),'k:','LineWidth',1);
        plot(SNR_dB,rates_nonideal_af_v(:,kk),'k','LineWidth',1);
end

xlabel('Average SNR [dB]');
ylabel('Ergodic Capacity [bits/channel use]');

legend('Ergodic Capacity','Approximation in (39)','Upper Bound in (33)','Capacity Ceiling','Location','SouthEast')

axis([0 50 0 5])
