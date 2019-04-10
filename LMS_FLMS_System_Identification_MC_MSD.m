% LMS and fractional LMS for System Identification

clear
close all
clc

tic

num = 2000; % number of measurements
dim = 8; % dimension of the parameters
gain = 10; % gain
w = rand(dim,1)*1*gain; % generate weight
 
len = num + dim - 1; % number of input signals

% num2 = num;
% len2 = num2 + dim - 1; % number of input signals

% mu1 = 0.0005;
% mu2 = 0.0005;
% fp = 0.5; % fractional power
mu = 0.01; % LMS
% mu3 = 0.005;

% ----------------------------------------
% muEq = mu1 + mu2*mean(w.^(1-fp)/gamma(2-fp));
%

Wn = ones(dim,1)*1*gain*2; % LMS
e = zeros(1,num);
Err = zeros(1,num);
Wni = Wn;

% Wn2 =  ones(dim,1)*1*gain*2; % F-LMS
% e2 = zeros(1,num);
% Err2 = zeros(1,num);
% Wn2i = Wn2;

% Wn3 =  ones(dim,1)*1*gain*2; % MCC
% e3 = zeros(1,num);
% Err3 = zeros(1,num);
% Wn3i = Wn3;

    % ----------------------------------------------------------
    % generate signal
    input = randn(len,1); % input signal
    x = zeros(dim,num); 
    for it = 1:num
        x(:,it) = input(it:it+dim-1,1); % signal in matrix form
    end
    noise = randn(1,num);
    snr = 20; % in dB scale
    d = w.'*x + 10^(-snr/20)*noise;
    % ----------------------------------------------------------
    Rxx = x*x'/num;
    wo = mu*10^(-snr/10)*trace(Rxx)/2;
%     woEq = muEq*10^(-snr/10)*trace(Rxx)/2; % Fractional LMS
 

runs = 1000;

for itr = 1:runs
    disp(['runs: ', num2str(itr)])
    % ----------------------------------------------------------
    % generate signal
    input = randn(len,1); % input signal
    x = zeros(dim,num); 
    for it = 1:num
        x(:,it) = input(it:it+dim-1,1); % signal in matrix form
    end
    noise = randn(1,num);
    snr = 20; % in dB scale
    d = w.'*x + 10^(-snr/20)*noise;
    % ----------------------------------------------------------
%     Rxx = x*x'/num;
%     wo = mu*10^(-snr/10)*trace(Rxx)/2;

    for itn = 1:num
        % LMS algorithm
        e(1,itn) = d(1,itn) - Wn(:,itn)'*x(:,itn);
        Wn(:,itn+1) = Wn(:,itn) + mu*e(1,itn)*x(:,itn);
        Err(itr,itn) = norm(Wn(:,itn+1) - w)^2;
    
%         % Fractional LMS algorithm - REAL
% %             Wn2(:,itn) = real(Wn2(:,itn));
%         e2(1,itn) = d(1,itn) - Wn2(:,itn)'*x(:,itn);
%         temF = Wn2(:,itn).^(1-fp);
%         Wn2(:,itn+1) = Wn2(:,itn) + mu1*e2(1,itn)*x(:,itn) + mu2*e2(1,itn)*(x(:,itn).*temF)/gamma(2-fp); 
%         Err2(itr,itn) = norm(real(Wn2(:,itn+1)) - w(:,1))^2;
        
%         % 3. MCC
% %         gkw = 0.1; % Gaussian kernel width
%         e3(1,itn) = d(1,itn) - Wn3(:,itn)'*x(:,itn);
%         Wn3(:,itn+1) = Wn3(:,itn) + mu3* sign(e3(1,itn))*x(:,itn);
%         Err3(itr,itn) = norm(Wn3(:,itn+1) - w)^2;
        
    end
 
      
    Wn = Wni;
%     Wn2 = Wn2i;
%     Wn3 = Wn3i;
end
 
mErr = mean(Err); % 1. LMS
% mErr2 = mean(Err2); % 2. Fractional-LMS
% mErr3 = mean(Err3); % 3. MCC
mEopt(1,1:num) = wo;
% mEopt2(1,1:num) = woEq;

figure
px = 1:num;
semilogy(px,mErr,'r-',px,mEopt,'k:','linewidth',2) 
legend('LMS','LMS (theoretical)')
xlabel('iteration','fontsize',12)
ylabel('Mean squared deviation (MSD)','fontsize',12)
grid
%title(['w = [',num2str(w.'),' ], \mu = \mu_1 = \mu_2 = ',num2str(mu1),', v = ',num2str(fp)],'fontsize',12)
% title(['w \sim U(', num2str(gain),',',num2str(2*gain),'), ',num2str(dim),...
%     ' taps, \mu = ',num2str(mu),', \mu_1 = ',num2str(mu1),', \mu_2 = ',num2str(mu2),', v = ',num2str(fp),', ', num2str(runs), ' runs.'])
toc