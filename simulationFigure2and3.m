%This Matlab script can be used to generate Figure 2 and Figure 3 of the paper:
%
%Emil Björnson, Michail Matthaiou, Mérouane Debbah, "A New Look at Dual-Hop
%Relaying: Performance Limits with Hardware Impairments" accepted for
%publication in IEEE Transcations on Communications. 
%
%The paper is available for free download: http://arxiv.org/pdf/1311.2634
%
%The code can also generate Figure 2 in "On the impact of transceiver
%impairments on AF relaying" by Emil Björnson, Agisilaos Papadogiannis,
%Michail Matthaiou, Mérouane Debbah, published at IEEE ICASSP 2013.
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


%Relaying with transceiver hardware impairments

%Number of channel iterations. This must be a very large number to make the
%Monte Carlo simulations fit the analytic curves for small BERs, e.g. 1e-6
nbrOfIterations = 1000000;

%Two scenarios: Rayleigh fading (1), Nakagami-m fading (2)
scenario = 2;

%Select fading parameters
if scenario == 1 %Rayleigh fading scenario
    alpha1 = 1; %Not used in Rayleigh fading distribution
    beta1 = 1;  %Variance at first hop
    
    alpha2 = 1; %Not used in Rayleigh fading distribution
    beta2 = 1;  %Variance at second hop
    
    %Generate fading realizations. RHO is the squared norm.
    RHO = [beta1*exprnd(1,1,nbrOfIterations); beta2*exprnd(1,1,nbrOfIterations)];
    
elseif scenario == 2 %Nakagami-m fading scenario
    
    alpha1 = 2; %alpha fading parameter at first hop
    beta1 = 1;  %beta fading parameter at first hop
    
    alpha2 = 2; %alpha fading parameter at second hop
    beta2 = 1;  %beta fading parameter at second hop
    
    %Generate fading realizations. RHO is the squared norm.
    RHO = [gamrnd(alpha1,beta1,1,nbrOfIterations); gamrnd(alpha2,beta2,1,nbrOfIterations)];
    
end


R = [2 5];   %Transmit information rate (per 2 channel uses)
x = 2.^R-1;  %Corresponding outage thresholds

N1 = 1;      %Normalized noise variance at relay
N2 = 1;      %Normalized noise variance at destination

SNR_dB = 0:5:40; %SNR range at first hop
P1_dB = SNR_dB - 10*log10(alpha1*beta1); %Corresponding transmit power range at source (in dB)
P1 = 10.^(P1_dB/10); %Corresponding transmit power range at source

P2 = P1; %Set same transmit power range at relay.

kappa1 = 0.1; %Error Vector Magnitude (EVM) at first hop
kappa2 = 0.1; %Error Vector Magnitude (EVM) at second hop
d = kappa1^2+kappa2^2+kappa1^2*kappa2^2; %Recurring function of kappa1 and kappa2


%Placeholders for Monte-Carlo results: Amplify-and-Forward (AF) case
OP_ideal_af_f = zeros(length(x),length(P1)); %Outage probabilities, AF fixed gain, ideal hardware
OP_ideal_af_v = zeros(length(x),length(P1)); %Outage probabilities, AF variable gain, ideal hardware

OP_nonideal_af_f = zeros(length(x),length(P1)); %Outage probabilities, AF fixed gain, non-ideal hardware
OP_nonideal_af_v = zeros(length(x),length(P1)); %Outage probabilities, AF variable gain, non-ideal hardware

%Placeholders for Monte-Carlo results: Decode-and-Forward (DF) case
OP_ideal_df = zeros(length(x),length(P1));    %Outage probabilities, DF, ideal hardware
OP_nonideal_df = zeros(length(x),length(P1)); %Outage probabilities, DF, non-ideal hardware


%This product will be used in relay gains
exp_rho1 = alpha1*beta1;


%Go through all cases of transmit power at the source and assume that the
%power on the relay is varied in the same way. Compute outage probabilities
%by Monte-Carlo simulations.
for p = 1:length(P1)
    
    G_ideal_f = sqrt(P2(p)/(P1(p)*exp_rho1 + N1));                % Ideal case Fixed Gain
    G_nonideal_f = sqrt(P2(p)/(P1(p)*exp_rho1*(1 + kappa1^2) + N1));    % With impairments Fixed Gain
    
    %Monte-Carlo simulations
    for kk = 1:nbrOfIterations
        
        %Extract channel gain for current realization
        rho1 = RHO(1,kk);
        rho2 = RHO(2,kk);
        
        %Compute signal-to-noise-and-distortion ratios (SNDRs) with AF relaying
        SNDR_ideal_af_f = (rho1*rho2)/(rho2*N1/P1(p) + N2/G_ideal_f^2/P1(p)); %Ideal hardware, fixed gain
        SNDR_ideal_af_v = (rho1*rho2)/(rho2*N1/P1(p) + N2*(rho1/P2(p)+N1/P2(p)/P1(p))); %Ideal hardware, variable gain
        
        SNDR_nonideal_af_f = (rho1*rho2)/(rho1*rho2*d + N1*rho2/P1(p) *(1+kappa2^2) + N2/G_nonideal_f^2/P1(p)); %Non-ideal hardware, fixed gain
        SNDR_nonideal_af_v = (rho1*rho2)/(rho1*rho2*d + rho1*N2*(1+kappa1^2)/P2(p) + rho2*N1*(1+kappa2^2)/P1(p) + N1*N2/P1(p)/P2(p));  %Non-ideal hardware, variable gain
        
        %Compute signal-to-noise-and-distortion ratios (SNDRs) with DF relaying
        SNDR_ideal_df = min([(rho1*P1(p))/N1,(rho2*P2(p))/N2]); %Ideal hardware
        SNDR_nonideal_df = min([(rho1*P1(p))/(P1(p)*rho1*kappa1^2 + N1),(rho2*P2(p))/(P2(p)*rho2*kappa2^2 + N2)]); %Non-ideal hardware
        
        
        %Check if the relaying channel is in outage (if SNDR < x). Store
        %the number of outages in corresponding variables
        for ind = 1:length(x)
            
            if SNDR_ideal_af_f < x(ind) %Fixed gain AF relaying, ideal hardware
                OP_ideal_af_f(ind,p) = OP_ideal_af_f(ind,p) + 1;
            end
            
            if SNDR_nonideal_af_f < x(ind) %Fixed gain AF relaying, non-ideal hardware
                OP_nonideal_af_f(ind,p) = OP_nonideal_af_f(ind,p) + 1;
            end
            
            if SNDR_ideal_af_v < x(ind) %Variable gain AF relaying, ideal hardware
                OP_ideal_af_v(ind,p) = OP_ideal_af_v(ind,p) + 1;
            end
            
            if SNDR_nonideal_af_v < x(ind) %Variable gain AF relaying, non-ideal hardware
                OP_nonideal_af_v(ind,p) = OP_nonideal_af_v(ind,p) + 1;
            end
            
            if SNDR_ideal_df < x(ind) %DF relaying, ideal hardware
                OP_ideal_df(ind,p) = OP_ideal_df(ind,p) + 1;
            end
            
            if SNDR_nonideal_df < x(ind) %DF relaying, non-ideal hardware
                OP_nonideal_df(ind,p) = OP_nonideal_df(ind,p) + 1;
            end
        end
        
    end
end



%Compute outage probabilities using analytic expressions

%Use same transmit power ranges as before but evaluate the performance for
%a denser set of values (every integer in dB).
SNR_dB_dense= min(SNR_dB):max(SNR_dB);
P1_dB_dense = SNR_dB_dense - 10*log10(alpha1*beta1);
P1 = 10.^(P1_dB_dense/10);
P2 = P1;

%Placeholders for analytic results: Amplify-and-Forward (AF) case
Pout_ideal_af_f = zeros(length(x),length(P1)); %Outage probabilities, AF fixed gain, ideal hardware
Pout_ideal_af_v = zeros(length(x),length(P1)); %Outage probabilities, AF variable gain, ideal hardware

Pout_nonideal_af_f = zeros(length(x),length(P1)); %Outage probabilities, AF fixed gain, non-ideal hardware
Pout_nonideal_af_v = zeros(length(x),length(P1)); %Outage probabilities, AF variable gain, non-ideal hardware

%Placeholders for analytic results: Decode-and-Forward (DF) case
Pout_ideal_df = zeros(length(x),length(P1));    %Outage probabilities, DF, ideal hardware
Pout_nonideal_df = zeros(length(x),length(P1)); %Outage probabilities, DF, non-ideal hardware


if scenario == 1
    Omega1 = beta1;    %New notation for channel variance at first hop
    Omega2 = beta2;    %New notation for channel variance at second hop
    
    
    %AF fixed gain relaying with ideal hardware. Outage probability is
    %given by [32]
    G_ideal_f = sqrt(P2./(P1*Omega1+N1));
    c_ideal_f = N2 ./ ( P1 .*(G_ideal_f.^2)*Omega1*Omega2 );
    
    for ind = 1:length(x)
        Pout_ideal_af_f(ind,:) = 1 - 2 * exp(-N1*x(ind)./(P1*Omega1)) .* sqrt(c_ideal_f*x(ind)) .* besselk(1, 2*sqrt(c_ideal_f*x(ind)) );
    end
    
    
    %AF variable gain relaying with ideal hardware. Outage probability is
    %given by [32]
    c_ideal_v = (N1*N2)./(P1.*P2*Omega1*Omega2);
    
    for ind = 1:length(x)
        Pout_ideal_af_v(ind,:) = 1 - 2 * exp(-x(ind)*( N1./(P1*Omega1)+N2./(P2*Omega2))) .* sqrt(c_ideal_v*(x(ind)+x(ind)^2)) .* besselk(1, 2*sqrt(c_ideal_v*(x(ind)+x(ind)^2)) );
    end
    
    
    %AF fixed gain relaying with non-ideal hardware. Outage probability is
    %given by Theorem 1
    G_nonideal_f = sqrt(P2./(P1*Omega1*(1+kappa1^2)+N1));
    A_nonideal_f = N2./(P1.* (G_nonideal_f.^2) *Omega1*Omega2);
    B_nonideal_f = N1*(1+kappa2^2)./(Omega1*P1);
    d = (kappa1^2+kappa2^2+kappa1^2*kappa2^2);
    
    for ind = 1:length(x)
        if d*x(ind)>1
            Pout_nonideal_af_f(ind,:) = ones(size(B_nonideal_f));
        else
            Pout_nonideal_af_f(ind,:) = 1 - 2 * exp(-B_nonideal_f*x(ind)/(1-d*x(ind))) .* sqrt(A_nonideal_f*x(ind)/(1-d*x(ind))) .* besselk(1, 2*sqrt(A_nonideal_f*x(ind)/(1-d*x(ind))) );
        end
    end
    
    
    %AF variable gain relaying with non-ideal hardware. Outage probability
    %is given by Theorem 1
    A_nonideal_v = (N1*N2)./(P1.* P2 *Omega1*Omega2);
    B_nonideal_v = N1*(1+kappa2^2)./(Omega1*P1)+N2*(1+kappa1^2)./(Omega2*P2);
    
    for ind = 1:length(x)
        if d*x(ind)>1
            Pout_nonideal_af_v(ind,:) = ones(size(B_nonideal_v));
        else
            Pout_nonideal_af_v(ind,:) = 1 - 2 * exp(-B_nonideal_v*x(ind)/(1-d*x(ind))) .* sqrt(A_nonideal_v*(x(ind)+x(ind)^2))/(1-d*x(ind)) .* besselk(1, 2*sqrt(A_nonideal_v*(x(ind)+x(ind)^2))/(1-d*x(ind)) );
        end
    end
    
    
    %DF relaying with ideal hardware. Outage probability is given in [33, Eq. (21)]
    for ind = 1:length(x)
        prefactor1 = (N1*x(ind)/beta1);
        prefactor2 = (N2*x(ind)/beta2);
        Pout_ideal_df(ind,:) = 1 - exp(-prefactor1./P1 - prefactor2./P2);
    end
    
    
    %DF relaying with non-ideal hardware. Outage probability is given by Theorem 2
    for ind = 1:length(x)
        prefactor1 = (N1*x(ind)/(1-kappa1^2*x(ind))/beta1);
        prefactor2 = (N2*x(ind)/(1-kappa2^2*x(ind))/beta2);
        
        delta = max([kappa1^2 kappa2^2]);
        
        if delta*x(ind)>1
            Pout_nonideal_df(ind,:) = ones(size(B_nonideal_v));
        else
            Pout_nonideal_df(ind,:) = 1 - exp(-prefactor1./P1 - prefactor2./P2);
        end
    end
    
elseif scenario == 2
    
    %AF fixed gain relaying with ideal hardware. Outage probability is
    %given by [32]
    G_ideal_f = sqrt(P2./(P1*alpha1*beta1+N1));
    b1 = 0;
    b2 = N1./P1;
    c = N2./(P1.*G_ideal_f.^2);
    d = 0;
    
    for ind = 1:length(x)
        Pout_ideal_af_f(ind,:) = functionOutageFormula(x(ind),alpha1,alpha2,beta1,beta2,b1,b2,c,d);
    end
    
    
    %AF variable gain relaying with ideal hardware. Outage probability is
    %given by [32]
    b1 = N2./P2;
    b2 = N1./P1;
    c = N1*N2./(P1.*P2);
    d = 0;
    
    for ind = 1:length(x)
        Pout_ideal_af_v(ind,:) = functionOutageFormula(x(ind),alpha1,alpha2,beta1,beta2,b1,b2,c,d);
    end
    
    
    %AF fixed gain relaying with non-ideal hardware. Outage probability is
    %given by Theorem 1
    G_nonideal_f = sqrt(P2./(P1*alpha1*beta1*(1+kappa1^2)+N1));
    b1 = 0;
    b2 = (1+kappa2^2)*N1./P1;
    c = N2./(P1.*G_nonideal_f.^2);
    d = kappa1^2+kappa2^2+kappa1^2*kappa2^2;
    
    for ind = 1:length(x)
        Pout_nonideal_af_f(ind,:) = functionOutageFormula(x(ind),alpha1,alpha2,beta1,beta2,b1,b2,c,d);
    end
    
    
    %AF variable gain relaying with non-ideal hardware. Outage probability
    %is given by Theorem 1
    b1 = (1+kappa1^2)*N2./P2;
    b2 = (1+kappa2^2)*N1./P1;
    c = N1*N2./(P1.*P2);
    d = kappa1^2+kappa2^2+kappa1^2*kappa2^2;
    
    for ind = 1:length(x)
        Pout_nonideal_af_v(ind,:) = functionOutageFormula(x(ind),alpha1,alpha2,beta1,beta2,b1,b2,c,d);
    end
    
    
    %DF relaying with ideal hardware. Outage probability is given in [33, Eq. (21)]
    for ind = 1:length(x)
        prefactor1 = (N1*x(ind)/beta1);
        prefactor2 = (N2*x(ind)/beta2);
        part1 = zeros(size(P1));
        part2 = zeros(size(P1));
        
        for j = 0:alpha1-1
            part1 = part1 + exp(-prefactor1./P1) .* (prefactor1./P1).^j / factorial(j);
        end
        
        for j = 0:alpha2-1
            part2 = part2 + exp(-prefactor2./P2) .* (prefactor2./P2).^j / factorial(j);
        end
        
        Pout_ideal_df(ind,:) = 1 - part1 .* part2;
    end
    
    
    %DF relaying with non-ideal hardware. Outage probability is given by Theorem 2
    for ind = 1:length(x)
        prefactor1 = (N1*x(ind)/(1-kappa1^2*x(ind))/beta1);
        prefactor2 = (N2*x(ind)/(1-kappa2^2*x(ind))/beta2);
        
        part1 = zeros(size(P1));
        part2 = zeros(size(P1));
        
        for j = 0:alpha1-1
            part1 = part1 + exp(-prefactor1./P1) .* (prefactor1./P1).^j / factorial(j);
        end
        
        for j = 0:alpha2-1
            part2 = part2 + exp(-prefactor2./P2) .* (prefactor2./P2).^j / factorial(j);
        end
        
        delta = max([kappa1^2 kappa2^2]);
        
        if delta*x(ind) > 1
            Pout_nonideal_df(ind,:) = ones(size(B_nonideal_v));
        else
            Pout_nonideal_df(ind,:) = 1 - part1 .* part2;
        end
    end
    
end


%Plot the results for AF relaying
figure; hold on; box on;

for ind = 1:length(x)
    plot(SNR_dB(1),OP_nonideal_af_f(ind,1)./nbrOfIterations,'k+-');
    plot(SNR_dB(1),OP_ideal_af_f(ind,1)./nbrOfIterations,'rd-.');
    plot(SNR_dB(1),OP_nonideal_af_v(ind,1)./nbrOfIterations,'bs--');
    plot(SNR_dB(1),OP_ideal_af_v(ind,1)./nbrOfIterations,'go:');
    
    plot(SNR_dB,OP_nonideal_af_f(ind,:)./nbrOfIterations,'k+');
    plot(SNR_dB,OP_ideal_af_f(ind,:)./nbrOfIterations,'rd');
    plot(SNR_dB,OP_nonideal_af_v(ind,:)./nbrOfIterations,'bs');
    plot(SNR_dB,OP_ideal_af_v(ind,:)./nbrOfIterations,'go');
    
    plot(SNR_dB_dense,Pout_nonideal_af_f(ind,:),'k');
    plot(SNR_dB_dense,Pout_ideal_af_f(ind,:),'r-.');
    plot(SNR_dB_dense,Pout_nonideal_af_v(ind,:),'b--');
    plot(SNR_dB_dense,Pout_ideal_af_v(ind,:),'g:');
end

set(gca,'yscale','log')

legend('Impairments (Fixed Gain)','Ideal Hardware (Fixed Gain)','Impairments (Variable Gain)','Ideal Hardware (Variable Gain)','Location','SouthWest')

xlabel('Average SNR [dB]');
ylabel('Outage Probability (OP)');
axis([0 40 1e-6 1]);



%Plot the results for DF relaying
figure; hold on; box on;

for ind = 1:length(x)
    plot(SNR_dB(1),OP_nonideal_df(ind,1)./nbrOfIterations,'rs-.');
    plot(SNR_dB(1),OP_ideal_df(ind,1)./nbrOfIterations,'ko-');
    
    plot(SNR_dB,OP_nonideal_df(ind,:)./nbrOfIterations,'rs');
    plot(SNR_dB,OP_ideal_df(ind,:)./nbrOfIterations,'ko');
    
    plot(SNR_dB_dense,Pout_nonideal_df(ind,:),'r-.');
    plot(SNR_dB_dense,Pout_ideal_df(ind,:),'k');
end

set(gca,'yscale','log')

legend('Impairments (DF)','Ideal Hardware (DF)','Location','SouthWest')

xlabel('Average SNR [dB]');
ylabel('Outage Probability (OP)');
axis([0 40 1e-6 1]);
