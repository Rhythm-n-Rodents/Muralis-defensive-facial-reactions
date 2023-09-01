%% plot figure S2C

 Dataall=load('....Data_S2_ABC.mat')

figure
fs = 2e4;

period = [4*fs:50*fs];
for i=2:5

    Y = Dataall{i} ;
    minmaxnorm = 1-Y.Meanresonse.SmoothBreathing;
    minmaxnorm=(minmaxnorm-min(minmaxnorm(period)))./(max(minmaxnorm(period))-min(minmaxnorm(period)));
%     subplot (2,1,1)
%     plot(minmaxnorm(period))
    NormBreathing(i,:) = minmaxnorm;
%     hold on

    minmaxnorm = Y.Meanresonse.Smoothspikes;
    minmaxnorm=(minmaxnorm-min(minmaxnorm(period)))./(max(minmaxnorm(period))-min(minmaxnorm(period)));
%     subplot (2,1,2)
%     plot(minmaxnorm(period))
%     hold on

    NormSpike(i,:) = minmaxnorm;


    if isfield(Y.Meanresonse,'SmoothBreathing2')
        minmaxnorm = 1-Y.Meanresonse.SmoothBreathing2;
        minmaxnorm=(minmaxnorm-min(minmaxnorm(period)))./(max(minmaxnorm(period))-min(minmaxnorm(period)));
        NormBreathing(i+5,:) = minmaxnorm;
        minmaxnorm = Y.Meanresonse.Smoothspikes2;
        minmaxnorm=(minmaxnorm-min(minmaxnorm(period)))./(max(minmaxnorm(period))-min(minmaxnorm(period)));
        NormSpike(i+5,:) = minmaxnorm;


    end

end


figure
allbreathis.thisx = [];
allbreathis.thisy = [];
allbreathis.thisysmooth=[];
fs = 2e4;

yyaxis left
y = mean(NormSpike(:,period));
y = mean(NormSpike(:,:));

x= (1:numel(y))./fs*100-3/fs*100;
% x= (1:numel(y))./fs*100-15;
x= (1:numel(y))./fs-(numel(y)./fs)/2;
x= (1:numel(y))./fs-18;
allbreathis.tsmooth=[];

plot(x,y,'-b' )
yyaxis right
maxi=0;
for i=[2:5]

    Y = Dataall{i} ;
    fs = 2e4;
    Stim_time = Y.AmoniaStim.onset;
    Breathing_IBI = [0 diff(Y.Breathing.Locations) ];
    Breathingsig=Y.Breathing;

    for j=1:numel(Stim_time)

        theseperiod= Stim_time(j)./fs;
        breathing = Breathingsig.Locations(Breathingsig.Locations>theseperiod-15&Breathingsig.Locations<theseperiod+40);

        % subplot 211
        % plot(breathing-theseperiod,ones(size(theseperiod))+i,'k.')

        plot(breathing-theseperiod,Breathing_IBI((Breathingsig.Locations>theseperiod-15&Breathingsig.Locations<theseperiod+40)),'r.','markersize',15)

        allbreathis.thisx = [allbreathis.thisx breathing-theseperiod];
        allbreathis.thisy = [allbreathis.thisy Breathing_IBI((Breathingsig.Locations>theseperiod-15&Breathingsig.Locations<theseperiod+40))];
        % [y, Ty] = resample( Breathing_IBI((Breathingsig.Locations>theseperiod-15&Breathingsig.Locations<theseperiod+40)), breathing-theseperiod ,fs);
        Ty = linspace (-15,40,20000);
        xSpline = abs(interp1(breathing-theseperiod,Breathing_IBI((Breathingsig.Locations>theseperiod-15&Breathingsig.Locations<theseperiod+40)),Ty,'nearest'));

        allbreathis.thisysmooth(j+maxi,:) = [ xSpline];
        allbreathis.tsmooth(j+maxi,:) = [ Ty];

        hold on

    end
    maxi =size(allbreathis.thisysmooth,1);
end
xlim([-10 40])
ylim([0 20])
figurecorrectionOPTO(gca,[0.05, 0.05],12)
plot(Ty,medfilt1(nanmedian(allbreathis.thisysmooth),1000),'r')


%% plot figure S2F


 Dataall=load('....Data_S2_ABC.mat')

figure
fs = 2e4;
clear NormSpike
period = [4*fs:50*fs];
for i=2:5

    Y = Dataall{i} ;
    minmaxnorm = 1-Y.Meanresonse.SmoothBreathing;
    minmaxnorm=(minmaxnorm-min(minmaxnorm(period)))./(max(minmaxnorm(period))-min(minmaxnorm(period)));
    subplot (2,1,1)
    plot(minmaxnorm(period))
    NormBreathing(i,:) = minmaxnorm;
    hold on

    minmaxnorm = Y.Meanresonse.Smoothspikes;
    minmaxnorm=(minmaxnorm-min(minmaxnorm(period)))./(max(minmaxnorm(period))-min(minmaxnorm(period)));
    subplot (2,1,2)
    plot(minmaxnorm(period))
    hold on

    NormSpike(i,:) = minmaxnorm;


    if isfield(Y.Meanresonse,'SmoothBreathing2')
        minmaxnorm = 1-Y.Meanresonse.SmoothBreathing2;
        minmaxnorm=(minmaxnorm-min(minmaxnorm(period)))./(max(minmaxnorm(period))-min(minmaxnorm(period)));
        NormBreathing(i+5,:) = minmaxnorm;
        minmaxnorm = Y.Meanresonse.Smoothspikes2;
        minmaxnorm=(minmaxnorm-min(minmaxnorm(period)))./(max(minmaxnorm(period))-min(minmaxnorm(period)));
        NormSpike(i+5,:) = minmaxnorm;


    end

end


figure
allbreathis.thisx = [];
allbreathis.thisy = [];
allbreathis.thisysmooth=[];
fs = 2e4;

yyaxis left
y = mean(NormSpike(:,period));
y = mean(NormSpike(:,:));

x= (1:numel(y))./fs*100-3/fs*100;
% x= (1:numel(y))./fs*100-15;
x= (1:numel(y))./fs-(numel(y)./fs)/2;
x= (1:numel(y))./fs-18;
allbreathis.tsmooth=[];

plot(x,y,'-b' )
yyaxis right
maxi=0;
for i=[1:5]

    Y = Dataall{i} ;
    fs = 2e4;
    Stim_time = Y.AmoniaStim.onset;
    Breathing_IBI = [0 diff(Y.Breathing.Locations) ];
    Breathingsig=Y.Breathing;

    for j=1:numel(Stim_time)

        theseperiod= Stim_time(j)./fs;
        breathing = Breathingsig.Locations(Breathingsig.Locations>theseperiod-15&Breathingsig.Locations<theseperiod+40);

        % subplot 211
        % plot(breathing-theseperiod,ones(size(theseperiod))+i,'k.')

        plot(breathing-theseperiod,Breathing_IBI((Breathingsig.Locations>theseperiod-15&Breathingsig.Locations<theseperiod+40)),'r.','markersize',15)

        allbreathis.thisx = [allbreathis.thisx breathing-theseperiod];
        allbreathis.thisy = [allbreathis.thisy Breathing_IBI((Breathingsig.Locations>theseperiod-15&Breathingsig.Locations<theseperiod+40))];
        % [y, Ty] = resample( Breathing_IBI((Breathingsig.Locations>theseperiod-15&Breathingsig.Locations<theseperiod+40)), breathing-theseperiod ,fs);
        Ty = linspace (-15,40,20000);
        xSpline = abs(interp1(breathing-theseperiod,Breathing_IBI((Breathingsig.Locations>theseperiod-15&Breathingsig.Locations<theseperiod+40)),Ty,'nearest'));

        allbreathis.thisysmooth(j+maxi,:) = [ xSpline];
        allbreathis.tsmooth(j+maxi,:) = [ Ty];

        hold on

    end
    maxi =size(allbreathis.thisysmooth,1);
end
xlim([-10 40])
ylim([0 20])
figurecorrectionOPTO(gca,[0.05, 0.05],12)
plot(Ty,medfilt1(nanmedian(allbreathis.thisysmooth),1000),'r')



