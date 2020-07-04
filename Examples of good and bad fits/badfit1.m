function play
%% TO SPECIFY 1:
% select column containing FI data and scale the data to fit the histogram
% provide file names to load data files, log10 transformation of data
% and fitting into histogram bins.
SC=3;
maxval_new_FC = 10;


X=load('1.csv'); [Xh(1,:),t] = hist(256*log10((X(:,SC))/maxval_new_FC),0:1023);
X=load('2.csv'); [Xh(2,:),t] = hist(256*log10((X(:,SC))/maxval_new_FC),0:1023);
X=load('3.csv'); [Xh(3,:),t] = hist(256*log10((X(:,SC))/maxval_new_FC),0:1023);
X=load('4.csv'); [Xh(4,:),t] = hist(256*log10((X(:,SC))/maxval_new_FC),0:1023);
X=load('5.csv'); [Xh(5,:),t] = hist(256*log10((X(:,SC))/maxval_new_FC),0:1023);
X=load('6.csv'); [Xh(6,:),t] = hist(256*log10((X(:,SC))/maxval_new_FC),0:1023);
X=load('7.csv'); [Xh(7,:),t] = hist(256*log10((X(:,SC))/maxval_new_FC),0:1023);
X=load('8.csv'); [Xh(8,:),t] = hist(256*log10((X(:,SC))/maxval_new_FC),0:1023);
X=load('9.csv'); [Xh(9,:),t] = hist(256*log10((X(:,SC))/maxval_new_FC),0:1023);
X=load('10.csv'); [Xh(10,:),t] = hist(256*log10((X(:,SC))/maxval_new_FC),0:1023);
X=load('11.csv'); [Xh(11,:),t] = hist(256*log10((X(:,SC))/maxval_new_FC),0:1023);
X=load('12.csv'); [Xh(12,:),t] = hist(256*log10((X(:,SC))/maxval_new_FC),0:1023);
X=load('13.csv'); [Xh(13,:),t] = hist(256*log10((X(:,SC))/maxval_new_FC),0:1023);
X=load('14.csv'); [Xh(14,:),t] = hist(256*log10((X(:,SC))/maxval_new_FC),0:1023);
X=load('15.csv'); [Xh(15,:),t] = hist(256*log10((X(:,SC))/maxval_new_FC),0:1023);

number_of_files=size(Xh,1);

%% TO SPECIFY 2:
% enter input data for cell concentration, number of data
% points in each file and relative time at which the sample was meaured
cc=[
18781500
20363000
19621000
19086000
20569500
21124500
22486500
24128000
26328500
27831000
34287500
38629000
50150000
63900000
72675000
];

g=[
37563
40726
39242
38172
41139
42249
44973
48256
52657
55662
68575
77258
10030
12780
14535
]';

tt=[
0
16.08333333
16.83333333
17.58333333
18.33333333
19.08333333
19.83333333
20.58333333
21.33333333
22.08333333
22.83333333
23.58333333
24.33333333
25.08333333
25.83333333
];

all_cc=cc;
all_tt=tt;
f=cc./g';

%% TO SPECIFY 3
% specify which samples should be used to perform the fit
indx=[2 : 15]; % window for bigaussian fit
cc_indx=[2 : 15]; % window for extended fit to cellcounts
li= length(indx);

%% TO SPECIFY 4
% specify the way data is pre-treated and some initial parameters

method='moving'; no=30; % smoothing of input data
bg_cutoff=1; % background threshold cutoff; for machines that can export gated data like Accuri C6 it is not needed
noisemag=1; noisewidth=1; % background triangle subtraction; for machines that can export gated data like Accuri C6 it is not needed
weight_on_cc = 1; % weight on CC fit; weight on bigaussian fit is 1
dye_degr=0.015; % 1/h - dye degradation rate
evaporation_rate=0.117; % ml/h - water evaporation rate in the flask 
ng_mu_upperbound=0.05; % 1/h - upper bound of non-growing population growth rate
ng_mu_lowerbound=0.0001; % 1/h - lower bound of non-growing population growth rate

%% Cell count correction due to evaporation
V0=50; %ml
V=V0-evaporation_rate*tt; %ml
cc=cc.*V./V0;
V=V0-evaporation_rate*all_tt; %ml
all_cc=all_cc.*V./V0;

%% determine mu (growth rate) from cell counts alone
tmp = fitoptions('Method','NonlinearLeastSquares',...
'Lower',[0,0],...
'Upper',[Inf,Inf],...
'Startpoint',[cc(indx(1)) 0.2]);
tmp2 = fittype('a*exp(b*x)','options',tmp);
[results,stats] = fit(tt(cc_indx),cc(cc_indx),tmp2)
mu_from_cc = results.b;

%% Smooth FI data and cut off background below threshold
for i=1:number_of_files,
    Xs(i,:)=smooth(Xh(i,bg_cutoff:end-1),no,method);
end
t = t(bg_cutoff:end-1);
I = 10.^(t/256); % convert t to I

%% Subtract background triangle from FI data above background cutoff
figure,hold on
title('Background correction check')
plot (Xs(1,:),'k:')
plot (Xs(indx(1),:),'b:')
noisex=1:noisewidth;
noisey=noisemag-noisemag/noisewidth*noisex;
for i=1:size(Xs,1),
    Xs(i,noisex)=Xs(i,noisex)-noisey;
    Xs(i,noisex)=max(0,Xs(i,noisex));
end
plot (Xs(1,:),'k')
plot (Xs(indx(1),:),'b')
plot (noisey,'r')
legend('t=0h','first gaussian t','t=0h','first gaussian t')

%% Scale FI data with measured cell number
for i=1:size(Xs,1),
    Xs(i,:)=Xs(i,:).*f(i);
end

%% Find I_0 value
disp('Attempting I_0 localization...')
options = optimset('Display','final','TolFun',1e-5,'TolX',1e-5,'MaxFunEvals',Inf);
tmp=fmincon (@gaussmodel,[750 100 1000*f(1)],[],[],[],[],[200 30 0],[1000 200 1E11],[],options);
    function out=gaussmodel (in)
        model = normpdf(t(bg_cutoff:end-1),in(1),in(2))*in(3);
        exp = Xs(1,bg_cutoff:end-1);
        out=(model-exp)*(model-exp)';
    end
t0=tmp(1);
I0 = 10.^(t0/256);
figure,h=subplot(1,1,1);hold on
title('I_0 determination check')
plot (I,Xs(1,:),'b')
plot (I,normpdf(t,tmp(1),tmp(2))*tmp(3),'r','Linewidth',2)
set (h,'XScale','log')
xlabel('FI')
ylabel('Cell count')
%% Bigaussian fit
disp ('Attempting Bigaussian fit...')

%% SPECIFY 5
% initial guesses for parameters
IG=[
    5000    % nongrowing population initial FI mean
    100     % nongrowing population FI stdev
    2e+04  % nongrowing population cell concentration
    0.01        % growth rate of non-growing population
    50          % growing population FI stdev
    2e+05  % initial growing population cell concentration 
    0.5         % growth rate of the population that starts to grow normally
    1000         % growing population initial FI mean
    10          % unstained cell background autofluorescence
    ];
% bounds for parameters (lower / upper)
bounds= [
            5000   15000     % nongrowing population initial FI mean
            10       500      % nongrowing population FI stdev
            1E3       1E5   % nongrowing population cell concentration 
            ng_mu_lowerbound       ng_mu_upperbound % growth rate of non-growing population
            15      200     % growing population FI stdev
            1E3     1E9    % initial growing population cell concentration
            0.3      0.6 % growth rate of the population that starts to grow normally
            800       2000    % growing population initial FI mean
            1        250     % cell background autofluorescence
            ];

options = optimset('Display','final','TolFun',1e-5,'TolX',1e-5,'MaxFunEvals',Inf);
sol=fmincon(@SSQ,IG,[],[],[],[],bounds(:,1),bounds(:,2),[],options);
function ssq=SSQ(in)
    ng_mean0 = in(1);
    ng_stdev = in(2);
    ng_scale0= in(3);
    ng_mu    = in(4);
    gr_stdev = in(5);
    gr_scale0= in(6);
    gr_mu    = in(7);
    gr_mean0 = in(8);
    EI       = in(9);
    ssq=0;nr=0;
    for i=1:li, % fit bigaussians to FI
        nr=nr+1;
        ng_mean=EI+(ng_mean0-EI).*exp(-(ng_mu+dye_degr).*(tt(indx(i))-tt(indx(1))));
        ng_mean=256*log10(ng_mean); % convert from I to t
        model_non(nr,:)=normpdf(t,ng_mean,ng_stdev)*ng_scale0*exp(ng_mu*(tt(indx(i))-tt(indx(1))));
        gr_mean=EI+(gr_mean0-EI).*exp(-(gr_mu+dye_degr).*(tt(indx(i))-tt(indx(1))));
        gr_mean=256*log10(gr_mean); % convert from I to t
        model_g(nr,:)=normpdf(t,gr_mean,gr_stdev)*gr_scale0*exp(gr_mu*(tt(indx(i))-tt(indx(1))));
        data=Xs(indx(i),:);
        model=model_non(nr,:)+model_g(nr,:);
        ssq=ssq+1/cc(indx(i))*(model-data)*(model-data)';
    end
    fit_CC_exp = all_cc(cc_indx); % choose exp. CC data
    fit_CC_t   = all_tt(cc_indx); % give timepoints of exp. CC data
    fit_CC_model_non =  ng_scale0*exp(ng_mu*(fit_CC_t-tt(indx(1)))); % NOT cc_indx
    fit_CC_model_gr = gr_scale0*exp(gr_mu*(fit_CC_t-tt(indx(1)))); % NOT cc_indx
    fit_CC_model_sum = fit_CC_model_gr + fit_CC_model_non;
    fit_ssq = 1/all_cc(cc_indx).*(fit_CC_model_sum-fit_CC_exp)'*(fit_CC_model_sum-fit_CC_exp);
    ssq=ssq+weight_on_cc*fit_ssq;
end

in=sol;
ng_mean0 = in(1);
ng_stdev = in(2);
ng_scale0= in(3);
ng_mu    = in(4);
gr_stdev = in(5);
gr_scale0= in(6);
gr_mu    = in(7);
gr_mean0 = in(8);
EI       = in(9);
nr=0;
for i=1:length(indx),
    nr=nr+1;
    ng_mean=EI+(ng_mean0-EI).*exp(-(ng_mu+dye_degr).*(tt(indx(i))-tt(indx(1))));
    ng_mean=256*log10(ng_mean); % convert from I to t
    model_non(nr,:)=normpdf(t,ng_mean,ng_stdev)*ng_scale0*exp(ng_mu*(tt(indx(i))-tt(indx(1))));
    gr_mean=EI+(gr_mean0-EI).*exp(-(gr_mu+dye_degr).*(tt(indx(i))-tt(indx(1))));
    gr_mean=256*log10(gr_mean); % convert from I to t
    model_g(nr,:)=normpdf(t,gr_mean,gr_stdev)*gr_scale0*exp(gr_mu*(tt(indx(i))-tt(indx(1))));
    data=Xs(indx(i),:);
    model=model_non(nr,:)+model_g(nr,:);
end

%% Calculate alpha
xx = tt(indx);
It = EI+(gr_mean0-EI).*exp(-(gr_mu+dye_degr).*(xx-xx(1)));
areaG = gr_scale0.*exp(gr_mu.*(tt(indx)-tt(indx(1))));
celldiv = log((I0-EI)./(It-EI))./log(2);
alfa = areaG./2.^celldiv./cc(1);

%% Display output
disp ('-------------------------------------------------')
disp (['alpha     = ' num2str(alfa(1))])
disp (['mu_growin = ' num2str(gr_mu)])
disp (['mu_nongrw = ' num2str(ng_mu)])
disp (['sgma_grow = ' num2str(gr_stdev)])
disp (['sgma_nong = ' num2str(ng_stdev)])
disp (['I_0       = ' num2str(I0)])
disp (['EI        = ' num2str(EI)])
disp (['dye_degr  = ' num2str(dye_degr)])
disp (['BG cutoff = ' num2str(bg_cutoff)])
disp (['BG magnit = ' num2str(noisemag)])
disp (['BG width  = ' num2str(noisewidth)])
disp (['CC weight = ' num2str(weight_on_cc)])
disp (['Bigaus pt = ' num2str(tt(indx)')])
disp (['CC pts    = ' num2str(all_tt(cc_indx)')])
disp ('-------------------------------------------------')

%% Plot all bigaussian fits in one figure
figure,h=subplot(1,1,1);hold on, box on, grid off
title ('Bigaussian fit to FI data')
for i=indx, % exp data
    plot (I,Xs(i,:),'b')
end
nr=0;
for i=indx, % sum of gaussians
    nr=nr+1;
    plot (I,model_g(nr,:)+model_non(nr,:),'k','LineWidth',3)
end
nr=0;
for i=indx, % non-growing gaussians
    nr=nr+1;
    plot (I,model_non(nr,:),'r','LineWidth',2)
end
nr=0;
for i=indx, % growing gaussians
    nr=nr+1;
    plot (I,model_g(nr,:),'g','LineWidth',2)
end
set(h,'XScale','log')
set(h,'YScale','lin')
xlim([5 10000])
xlabel('FI')
ylabel('Cell count')

%% Plot bigaussian fits, one per subplot
figure, nr=0; 
sgtitle('Bigaussian fit for each time point')

for i=1:length(tt),
    map=[1 4 7 10 13 2 5 8 11 14 3 6 9 12 15];
    h=subplot(5,3,map(i));hold on, box on, grid off
    plot (I,Xs(i,:),'b')      % plot experimental data
    if any(i==indx),        % plot model fit where it exists
        nr=nr+1;
        plot (I,model_g(nr,:)+model_non(nr,:),'k','LineWidth',3)
        plot (I,model_non(nr,:),'r','LineWidth',2)
        plot (I,model_g(nr,:),'g','LineWidth',2)
    end
    set(h,'XScale','log')
    set(h,'YScale','lin')
    title(h, strcat('t = ', num2str(tt(i))))
    xlim([5 10000])
    xlabel('FI')
    ylabel('Cell count')
end

%% Plot all time points in one figure, without bigaussian fit
figure, h=subplot(1,1,1);hold on, box on, grid off
for i=1:length(tt),
    if i==1,
        color='b';
    else if any(i==indx),
            color='r';
        else 
            color='m';
        end
    end
    plot (I,Xs(i,:),color)
end
set(h,'XScale','log')
set(h,'YScale','lin')
xlim([5 I(end)])
title ('Fluorescence Data')
xlabel('FI')
ylabel('Cell count')

%% Plot cellcount curve
xx = tt(indx(1)):0.1:tt(indx(end));
yyg = gr_scale0.*exp(gr_mu.*(xx-tt(indx(1))));
yyn = ng_scale0.*exp(ng_mu.*(xx-tt(indx(1))));
xx_ext = all_tt(cc_indx(1)):0.1:all_tt(cc_indx(end));
yyg_ext = gr_scale0.*exp(gr_mu.*(xx_ext-tt(indx(1))));
yyn_ext = ng_scale0.*exp(ng_mu.*(xx_ext-tt(indx(1))));
figure, h=subplot(1,1,1); hold on
title('Cellcount curve fit check')
plot (xx_ext,yyg_ext,'m')
plot (xx_ext,yyn_ext,'c')
plot (xx_ext,yyg_ext+yyn_ext,'r')
plot (xx,yyg,'m','LineWidth',2)
plot (xx,yyn,'c','LineWidth',2)
plot (xx,yyg+yyn,'r','LineWidth',2)
scatter(all_tt,all_cc,'ko')
scatter(all_tt(cc_indx),all_cc(cc_indx),'ko','filled')
scatter(tt(indx),cc(indx),'ro','filled')

%scatter(3.5,2.248284837060616e+07,'ko') % Wert 2 und 3 von cc.m 
%scatter(10,2.812223684132589e+07,'ko') % Wert 2 und 3 von cc.m

set(h,'YScale','log')
xlabel('time')
ylabel('Cell count')

end


 


