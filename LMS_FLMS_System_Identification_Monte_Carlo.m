% LMS and fractional LMS for System Identification

clear all
close all
clc

tic

num = 6000; % number of measurements
dim = 26; % dimension of the parameters
gain = 10; % gain
w = (-10:1:dim-11).';
% w = (-1:0.1:1.5).';
% w = (1.1:0.1:0.1*(dim-1)+1.1).';
% w = gain - 2*gain*rand(dim,1);
% w = rand(dim,1);

% dim = 1;
% w = 0.1;
 
len = num + dim - 1; % number of input signals

num2 = num;
len2 = num2 + dim - 1; % number of input signals

mu1 = 0.0018;
mu2 = 0.0018;
fp = 0.5; % fractional power
mu = 0.005; % LMS

runs = 100;
% ----------------------------------------
muEq = mu1 + mu2*mean(w.^(1-fp)/gamma(2-fp)) 
%

Wn = ones(dim,1)*2*gain; % LMS
e = zeros(1,num);
Err = zeros(1,num);
Wni = Wn;

Wn2 =  ones(dim,1)*2*gain; % F-LMS with ABS
e2 = zeros(1,num);
Err2 = zeros(1,num);
Wn2i = Wn2;

Wn3 =  ones(dim,1)*2*gain; % F-LMS with REAL
e3 = zeros(1,num);
Err3 = zeros(1,num);
Wn3i = Wn3;

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
    woEq = muEq*10^(-snr/10)*trace(Rxx)/2;
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
        
        wn(:,itn,itr) = Wn(:,itn+1);
        % Fractional LMS algorithm - REAL
%             Wn2(:,itn) = real(Wn2(:,itn));
            
        Wn2(:,itn) = abs(Wn2(:,itn));
        e2(1,itn) = d(1,itn) - Wn2(:,itn)'*x(:,itn);
        temF = Wn2(:,itn).^(1-fp);
        Wn2(:,itn+1) = Wn2(:,itn) + mu1*e2(1,itn)*x(:,itn) + mu2*e2(1,itn)*(x(:,itn).*temF)/gamma(2-fp); 
        Err2(itr,itn) = norm(real(Wn2(:,itn+1)) - w(:,1))^2;
        wn2(:,itn,itr) = Wn2(:,itn+1);
        
        Wn3(:,itn) = real(Wn3(:,itn));
        e3(1,itn) = d(1,itn) - Wn3(:,itn)'*x(:,itn);
        temF3 = Wn3(:,itn).^(1-fp);
        Wn3(:,itn+1) = Wn3(:,itn) + mu1*e3(1,itn)*x(:,itn) + mu2*e3(1,itn)*(x(:,itn).*temF3)/gamma(2-fp); 
        Err3(itr,itn) = norm(real(Wn3(:,itn+1)) - w(:,1))^2;
        wn3(:,itn,itr) = Wn3(:,itn+1);
%         Eopt(itr,it) = wo;
    end
    

    

    
%       Wnit(itr,:) = Wn;
%       Wn2it(itr,:)= Wn2;
      
%       t1 = Wn2it(itr,:).^(1-fp);
%       t2 = Wn2it(itr,:).^(2-fp);
%       Rdx = mean(ones(dim,1)*d.*x,2);
%       Rxx = x*x'/num;
%       wot = Rxx\Rdx;
%       wo2(itr,:) = (wot.*(1+t1/gamma(2-fp)) - t2/gamma(2-fp));
      
%           Erro(itr,:) = norm(wo2(1,itr) - w)^2; 
      
    Wn = Wni;
    Wn2 = Wn2i;
    Wn3 = Wn3i;
end

iwn = mean(wn,3);
bwn = iwn-w*ones(1,num);

iwn2 = mean(wn2,3);
bwn2 = iwn2-w*ones(1,num);

iwn3 = mean(wn3,3);
bwn3 = iwn3-w*ones(1,num);

% figure
% subplot(121)
% plot(iwn.','linewidth',2),grid
% xlabel('iteration','fontsize',12)
% ylabel('Instantaneous weights','fontsize',12)
% % axis([1 num 8 22])
% title('LMS','fontsize',11)
% 
% subplot(122)
% plot(bwn.','linewidth',2),grid
% xlabel('iteration','fontsize',12)
% ylabel('Bias','fontsize',12)
% % axis([1 num 8 22])
% title('LMS','fontsize',11)

% figure
% subplot(121)
% plot(iwn2.','linewidth',2),grid
% xlabel('iteration','fontsize',12)
% ylabel('Instantaneous weights','fontsize',12)
% % axis([1 num 8 22])
% title('Fractional LMS','fontsize',11)
% 
% subplot(122)
% plot(bwn2.','linewidth',2),grid
% xlabel('iteration','fontsize',12)
% ylabel('Bias','fontsize',12)
% % axis([1 num 8 22])
% title('Fractional LMS','fontsize',11)

% ----------------------------------------------
figure

px = 1:num;
% temp = (wn-wn2).*(wn-wn2);
% temp2 = mean(temp,3);
% xy = (diag(temp2.'*temp2));

temp = (wn-wn2).*(wn-wn2);
temp2 = mean(temp,3);

if dim~=1
    xy = (mean(temp2));
else
    xy = (temp2);
end

tempb = (wn-wn3).*(wn-wn3);
temp2b = mean(tempb,3);

if dim~=1
    xyb = (mean(temp2b));
else
    xyb = (temp2b);
end

semilogy(px,real(xy),'-.',px,real(xyb),'--','linewidth',2) 
xlabel('number of iterations','fontsize',12)
ylabel('similarity measure','fontsize',12)
legend('FLMS (ABS)','FLMS (REAL)')
% axis([0 num 10^(-5) 10^3])
% ----------------------------------------------

% figure
% % subplot(212)
% plot(bwn2.','linewidth',2),grid
% xlabel('iteration','fontsize',12)
% ylabel('Bias','fontsize',12)
% % axis([1 num 8 22])
% title('Fractional LMS','fontsize',11)
% figure
% mWn = mean(Wnit);
% mWn2 = mean(Wn2it);
% mwo2 = 2*gain-mean(wo2);
% 
% px = 1:num;
% plot(px,mWn(1:num),'--',px,mWn2(1:num),'-',px,mwo2(1:num),':','linewidth',2),grid
% xlabel('iteration','fontsize',12)
% ylabel('weight','fontsize',12)
% title(['w =', num2str(gain),', ',num2str(dim),...
%     ' taps, \mu = ',num2str(mu),', \mu_1 = ',num2str(mu1),', \mu_2 = ',num2str(mu2),', v = ',num2str(fp),', ', num2str(runs), ' runs.'],'fontsize',11)
% legend('LMS','Fractional LMS','Fractional LMS (theoretical)')
% 
mErr = mean(Err);
mErr2 = mean(Err2);
mErr3 = mean(Err3);
mEopt(1,1:num) = wo;
% mEopt2(1,1:num) = norm(wo2 - w)^2;
% mEopt2(1,1:num) = mean(Erro);
mEopt2(1,1:num) = woEq;

figure
px = 1:num;
% semilogy(px,mErr2,'-.',px,real(mEopt2),'--',px,mErr,'r-',px,mEopt,'k:','linewidth',2) 
% legend('Fractional LMS','Fractional LMS (theoretical)','LMS','LMS (theoretical)')
semilogy(px,mErr2,'-.',px,mErr3,'--',px,mErr,'r-','linewidth',2) 
legend('FLMS (ABS)','FLMS (REAL)','LMS')
xlabel('number of iterations','fontsize',12)
ylabel('MSD','fontsize',12)
 
% title(['w \sim U(', num2str(0*gain),',',num2str(1*gain),'), ',num2str(dim),...
%     ' taps, \mu_1 = ',num2str(mu1),', \mu_2 = ',num2str(mu2),', \mu = ',num2str(mu),', v = ',num2str(fp),', ', num2str(runs), ' runs.'])

% title(['w = ', num2str(gain),',',num2str(dim),...
%     ' tap, \mu_1 = ',num2str(mu1),', \mu_2 = ',num2str(mu2),', \mu = ',num2str(mu),', v = ',num2str(fp),', ', num2str(runs), ' runs.'])

toc