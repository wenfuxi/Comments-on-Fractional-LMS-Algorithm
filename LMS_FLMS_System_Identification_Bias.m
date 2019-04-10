% LMS and fractional LMS for System Identification

clear all
close all
clc

tic

num = 20000; % number of measurements
dim = 16; % dimension of the parameters
gain = 10; % gain
w = gain + ones(dim,1)*1*gain; % generate weight
 
len = num + dim - 1; % number of input signals

num2 = num;
len2 = num2 + dim - 1; % number of input signals

mu1 = 0.0005;
mu2 = 0.0005;
fp = 0.5; % fractional power
mu = 0.0005; % LMS

% ----------------------------------------
muEq = mu1 + mu2*mean(w.^(1-fp)/gamma(2-fp));
%
    snrL = 10;
    snrS = 5;
    snrH = 30;

Wn = ones(dim,1)*1*gain; % LMS
e = zeros(1,num);
Err = zeros(1,length(snrL:snrS:snrH));
Wni = Wn;

Wn2 =  ones(dim,1)*1*gain; % F-LMS
e2 = zeros(1,num);
Err2 = zeros(1,length(snrL:snrS:snrH));
Wn2i = Wn2;
% 
%     % ----------------------------------------------------------
%     % generate signal
%     input = randn(len,1); % input signal
%     x = zeros(dim,num); 
%     for it = 1:num
%         x(:,it) = input(it:it+dim-1,1); % signal in matrix form
%     end
%     noise = randn(1,num);
%     snr = 20; % in dB scale
%     d = w.'*x + 10^(-snr/20)*noise;
%     % ----------------------------------------------------------
%     Rxx = x*x'/num;
%     wo = mu*10^(-snr/10)*trace(Rxx)/2;
%     woEq = muEq*10^(-snr/10)*trace(Rxx)/2;
    % ----------------------------------------------------------
    % generate signal
%     input2 = randn(len2,1); % input signal
%     x = zeros(dim,num2); 
%     for it = 1:num2
%         x2(:,it) = input2(it:it+dim-1,1); % signal in matrix form
%     end
%     noise2 = randn(1,num2);
%     snr = 20; % in dB scale
%     d2 = w.'*x2 + 10^(-snr/20)*noise2;
%     % ----------------------------------------------------------
%     t1 = w.^0.5;
%     t2 = w.^1.5;
%     Rdx = mean(ones(dim,1)*d2.*x2,2);
%     Rxx = x2*x2'/num2;
%     wot = Rxx\Rdx;
%     
%     wo2 = (wot.*(1+t1/gamma(1.5)) - t2/gamma(1.5));
    % ----------------------------------------------------------
    


runs = 100;

itsnr = 0;

for snr = snrL:snrS:snrH
    
    itsnr = itsnr + 1;
    disp(['snr: ', num2str(snr)])
    for itr = 1:runs
        
        % ----------------------------------------------------------
        % generate signal
        input = randn(len,1); % input signal
        x = zeros(dim,num); 
        for it = 1:num
            x(:,it) = input(it:it+dim-1,1); % signal in matrix form
        end
        noise = randn(1,num);
    %     snr = 20; % in dB scale
        d = w.'*x + 10^(-snr/10)*noise;
    % ----------------------------------------------------------
%     Rxx = x*x'/num;
%     wo = mu*10^(-snr/10)*trace(Rxx)/2;

        for itn = 1:num
            % LMS algorithm
            e(1,itn) = d(1,itn) - Wn(:,itn)'*x(:,itn);
            Wn(:,itn+1) = Wn(:,itn) + mu*e(1,itn)*x(:,itn);
%             Err(itr,itn) = norm(Wn(:,itn+1) - w)^2;
    
        % Fractional LMS algorithm - REAL
%             Wn2(:,itn) = real(Wn2(:,itn));
            e2(1,itn) = d(1,itn) - Wn2(:,itn)'*x(:,itn);
            temF = Wn2(:,itn).^(1-fp);
            Wn2(:,itn+1) = Wn2(:,itn) + mu1*e2(1,itn)*x(:,itn) + mu2*e2(1,itn)*(x(:,itn).*temF)/gamma(2-fp); 
%             Err2(itr,itn) = norm(real(Wn2(:,itn+1)) - w(:,1))^2;
        
%         Eopt(itr,it) = wo;
        end
    
%     t1 = w.^0.5;
%     t2 = w.^1.5;
    t1 = Wn2(:,end).^0.5;
    t2 = Wn2(:,end).^1.5;
    Rdx = mean(ones(dim,1)*d.*x,2);
    Rxx = x*x'/num;
    wot = Rxx\Rdx;
    wo2 = (wot.*(1+t1/gamma(1.5)) - t2/gamma(1.5));
    wo(itsnr,itr) = mean(wo2 - w); 
    
%     Wn = Wni;
%     Wn2 = Wn2i;
Err(itsnr,itr) = mean(Wn(:,end) - w);
Err2(itsnr,itr) = mean(Wn2(:,end) - w);
end
end

mErr = mean(Err,2);
mErr2 = mean(Err2,2);
mEopt = mean(wo,2); 
% mEopt2(1,1:num) = norm(wo2 - w)^2;
% mEopt2(1,1:num) = mean(Erro);
% mEopt2(1,1:length(snrL:snrS:snrH)) = woEq;

figure
px = snrL:snrS:snrH;
plot(px,abs(mErr2),'-.',px,abs(mErr),'r-',px,mEopt,'s-','linewidth',2),grid 
legend('Fractional LMS','LMS','Theoretical')
xlabel('SNR (dB)','fontsize',12)
ylabel('Mean squared deviation (MSD)','fontsize',12)
%title(['w = [',num2str(w.'),' ], \mu = \mu_1 = \mu_2 = ',num2str(mu1),', v = ',num2str(fp)],'fontsize',12)
% title(['w \sim U(', num2str(gain),',',num2str(2*gain),'), ',num2str(dim),...
%     ' taps, \mu = ',num2str(mu),', \mu_1 = ',num2str(mu1),', \mu_2 = ',num2str(mu2),', v = ',num2str(fp),', ', num2str(runs), ' runs.'])
toc