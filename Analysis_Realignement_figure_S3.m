%% 
% list of days:
ploteachtrial =1;
Fs =1e4; % sampling frequency

filename{1} = 'sept_22_2022_Rat332';
filename{2} = 'sept_22_2022_Rat335';
filename{3} = 'sept_23_2022_Rat333';
filename{4} = 'sept_23_2022_Rat336';
filename{5} = 'sept_23_2022_Rat337';
filename{6} = 'sept_26_2022_Rat332';
filename{7} = 'sept_26_2022_Rat335';
filename{8} = 'sept_26_2022_Rat336';
filename{9} = 'sept_26_2022_Rat337';



filename2{1} = '22sept22_Rat_332';
filename2{2} = '22sept22_Rat_335';
filename2{3} = '23sept22_Rat_333';
filename2{4} = '23sept22_Rat_336';
filename2{5} = '23sept22_Rat_337';
filename2{6} = '26sept22_Rat_332';
filename2{7} = '26sept22_Rat_335';
filename2{8} = '26sept22_Rat_336';
filename2{9} = '26sept22_Rat_337';

% timing of the stimulation each day
Stimtime{1}=[1445 2271 2684 4336 5162 5988];
Stimtime{2}=[619 1445 2271 3097 3923 4749 5575];
Stimtime{3}=[1032 3510];
Stimtime{4}=[1032 1858 2684 3510 4336 5162 5988 ];
Stimtime{5}=[3097 3923 4749 6401 nan];
Stimtime{6}=[206 1032 1858 2684 3510 4336 5162 5988];
Stimtime{7}=[nan 4336 3510 2684 1858 1032 206];
Stimtime{8}=[1858 2684 3510 5162 5988 ];
Stimtime{9}=[206 1032 1858 2684 3510 5162 5575 5988];





Air{1}=[1032 1858 3087 4749 5575 nan];
Air{2}=[1032 1858 2684 3510 4336 5162 5988 6401];
Air{3}=[1445 3923];
Air{4}=[1445 2271 3097 3923 4749 5575 6401 ];
Air{5}=[3510 4336 5162 nan];
Air{6}=[619 1445 2271 3097 3923 4749 5575 6401 ];
Air{7}=[nan nan 3923 3097 2271 1445 619 ];
Air{8}=[2271 3097 3923 5575 nan];
Air{9}=[619 1445 2271 3097 3923 nan nan 6401 ];



















%% digital filter designe ;
creatfilter = 1; % place 0 if your own digital filter 
if creatfilter==1
    Fpass = 50;
    Fstop = 70;
    Ap = 1;
    Ast = 30;
    Fs = 1e4;
    [A,B,C,D] = butter(10,[10 50]/1e4);
    d = designfilt('bandstopiir','FilterOrder',2, ...
        'HalfPowerFrequency1',59,'HalfPowerFrequency2',61, ...
        'DesignMethod','butter','SampleRate',Fs);

    lpFilt = designfilt('lowpassiir','FilterOrder',8, ...
        'PassbandFrequency',1,'PassbandRipple',0.2, ...
        'SampleRate',1e4);

    sos = ss2sos(A,B,C,D);
    fvt = fvtool(lpFilt,d,'Fs',1500);
    legend(fvt,'butter','designfilt')

end

%%


Sessioncounter=0;


for session=1:9

    for i=1:numel(Stimtime{session})
        Sessioncounter=Sessioncounter+1;
        switch i
            case 1
                Matdata =load(['\\dk-server.dk.ucsd.edu\data\mdeschenesvis\Grimace ammonia\Figure 3 Awake-rat measurement\Matlab segments corresponding to movie\' filename{session} '_s1.mat']);
            case 2
                Matdata =load(['\\dk-server.dk.ucsd.edu\data\mdeschenesvis\Grimace ammonia\Figure 3 Awake-rat measurement\Matlab segments corresponding to movie\' filename{session}  '_s2.mat']);
            case 3
                Matdata =load(['\\dk-server.dk.ucsd.edu\data\mdeschenesvis\Grimace ammonia\Figure 3 Awake-rat measurement\Matlab segments corresponding to movie\' filename{session} '_s3.mat']);
            case 4
                Matdata =load(['\\dk-server.dk.ucsd.edu\data\mdeschenesvis\Grimace ammonia\Figure 3 Awake-rat measurement\Matlab segments corresponding to movie\' filename{session} '_s4.mat']);
            case 5
                Matdata =load(['\\dk-server.dk.ucsd.edu\data\mdeschenesvis\Grimace ammonia\Figure 3 Awake-rat measurement\Matlab segments corresponding to movie\' filename{session} '_s5.mat']);
            case 6
                Matdata =load(['\\dk-server.dk.ucsd.edu\data\mdeschenesvis\Grimace ammonia\Figure 3 Awake-rat measurement\Matlab segments corresponding to movie\' filename{session} '_s6.mat']);
            case 7
                Matdata =load(['\\dk-server.dk.ucsd.edu\data\mdeschenesvis\Grimace ammonia\Figure 3 Awake-rat measurement\Matlab segments corresponding to movie\' filename{session} '_s7.mat']);
            case 8
                Matdata =load(['\\dk-server.dk.ucsd.edu\data\mdeschenesvis\Grimace ammonia\Figure 3 Awake-rat measurement\Matlab segments corresponding to movie\' filename{session} '_s8.mat']);
        end

        Data(Sessioncounter).thermistor2=Matdata.data(Matdata.datastart(1):Matdata.dataend(1));
        Data(Sessioncounter).whisker=Matdata.data(Matdata.datastart(2):Matdata.dataend(2));
        Data(Sessioncounter).Trigger=Matdata.data(Matdata.datastart(3):Matdata.dataend(3));
        Data(Sessioncounter).AirNh3=Matdata.data(Matdata.datastart(4):Matdata.dataend(4));
        Data(Sessioncounter).Trigger2=Matdata.data(Matdata.datastart(5):Matdata.dataend(5));
        Data(Sessioncounter).ratname=filename{session};
        [~, Data(Sessioncounter).onset]=findpeaks(Data(Sessioncounter).Trigger2,1e4,'MinPeakProminence',.5,'MinPeakDistance',5);
    end
end
%%

Sessioncounter=0;


for session=1:9

    %%
    % genfolder= 'D:\Dropbox\Expression of defensive facial reactions\Grimace ammonia\Figure 3 Awake-rat measurement\Movies\26september_rat335_20220926_161253\';
    % genfolder= 'D:\Matlab segments corresponding to movie\26september_rat335_20220926_161253\';
    vidname2 ='26september_rat335_20220926_131556_20220926_161253';
    Nose = load(['\\dk-server.dk.ucsd.edu\data\mdeschenesvis\Grimace ammonia\Figure 3 Awake-rat measurement\Nose tracking coordonnates\' filename2{session} '_Coord_nez.mat']);
    Eye=load(['\\dk-server.dk.ucsd.edu\data\mdeschenesvis\Grimace ammonia\Figure 3 Awake-rat measurement\Eye aperture tracking\' filename2{session} '_OuvertureOeil.mat']);
    % Ear =load(['\\dk-server.dk.ucsd.edu\data\mdeschenesvis\Grimace ammonia\Figure 3 Awake-rat measurement\Ears tracking\' filename2{session} '_Oreille.mat']);
    % vMain = VideoReader([genfolder   vidname2 '.avi']);
    % video = read(vMain);
    Coord_Nez=Nose.Coord_Nez;

    %%
    % figure
    % buttLoop = filtfilt(d,Data(1).thermistor2);
    % y = filtfilt(lpFilt,[zeros(1,1e5) buttLoop]);
    % filtered = ([zeros(1,1e5) buttLoop]-y);
    % filtered(1:1e5)=[];

    %%
    figure
    Lastframe=800;
    framestart=100;

    thisStimtime = Stimtime{session};
    thisAir= Air{session};

    % Airstim = [619 1445 2271 3097 3923 4749 5575 6401];
    Fcounter=0;

    for i=1:numel(thisStimtime)
        Sessioncounter=Sessioncounter+1;

        if isnan(thisStimtime(i))
            continue
        end
        subplot(3,1,1)
        plot(Coord_Nez(thisStimtime(i)-200:thisStimtime(i)+200,2))
        AllMoviemeasures.noseX(Sessioncounter,:)=Coord_Nez(thisStimtime(i)-200:thisStimtime(i)+200,2)';
        AllMoviemeasures.noseY(Sessioncounter,:)=Coord_Nez(thisStimtime(i)-200:thisStimtime(i)+200,3)';
        AllMoviemeasures.eyenorm_opening(Sessioncounter,:)=Eye.norm_opening(thisStimtime(i)-200:thisStimtime(i)+200);
        if ~isnan(thisAir(i))
            AllMoviemeasures.ControlnoseX(Sessioncounter,:)=Coord_Nez(thisAir(i)-200:thisAir(i)+200,2)';
            AllMoviemeasures.ControlnoseY(Sessioncounter,:)=Coord_Nez(thisAir(i)-200:thisAir(i)+200,3)';
            AllMoviemeasures.Controleyenorm_opening(Sessioncounter,:)=Eye.norm_opening(thisAir(i)-200:thisAir(i)+200);
        else

            AllMoviemeasures.ControlnoseX(Sessioncounter,:)=nan(1,401);
            AllMoviemeasures.ControlnoseY(Sessioncounter,:)=nan(1,401);
            AllMoviemeasures.Controleyenorm_opening(Sessioncounter,:)=nan(1,401);

        end


        % AllMoviemeasures.ear(Sessioncounter,:)=Ear.PPLateralSmooth_norm(thisStimtime(i)-200:thisStimtime(i)+200);

        hold on
        subplot(3,1,2)
        thisstart =round(Data(Sessioncounter).onset(1)*1e4);
        buttLoop = filtfilt(d,Data(Sessioncounter).thermistor2);
        y = filtfilt(lpFilt,[zeros(1,1e5) buttLoop]);
        filtered = ([zeros(1,1e5) buttLoop]-y);
        filtered(1:1e5)=[];
        X=(thisstart-5e4:thisstart+5e4)./1e4;
        plot(X-thisstart/1e4,Data(Sessioncounter).thermistor2(thisstart-5e4:thisstart+5e4))
        hold on
        subplot(3,1,3)
        [peaks, indexes] = findpeaks(filtered,1e4,'MinPeakProminence',.2,'MinPeakDistance',0.15,'MinPeakHeight',0.1);
        thesepeaks = indexes(indexes>(thisstart-5e4)/1e4&indexes<(thisstart+5e4)/1e4);
        plot(thesepeaks-(thisstart)/1e4,ones(size(thesepeaks))+i,'.k','markersize',12)
        hold on
        plot(X-((thisstart)/1e4),(filtered(thisstart-5e4:thisstart+5e4)./5)+i+1,'color',[0.5 0.5 0.5])
        AllBreathing.Breathingtimes{Sessioncounter} =thesepeaks-(thisstart)/1e4;
        AllBreathing.Breathingsignal(Sessioncounter,:)=filtered(thisstart-5e4:thisstart+5e4);
        AllBreathing.Breathingsignalraw(Sessioncounter,:)=Data(Sessioncounter).thermistor2(thisstart-5e4:thisstart+5e4);
        AllBreathing.whiskerraw(Sessioncounter,:)=Data(Sessioncounter).whisker(thisstart-5e4:thisstart+5e4);




        if numel(Data(Sessioncounter).onset)>1
            thisstart =round(Data(Sessioncounter).onset(2)*1e4);

            X=(thisstart-5e4:thisstart+5e4)./1e4;
            plot(X-thisstart/1e4,Data(Sessioncounter).thermistor2(thisstart-5e4:thisstart+5e4))
            hold on
            subplot(3,1,3)
            thesepeaks = indexes(indexes>(thisstart-5e4)/1e4&indexes<(thisstart+5e4)/1e4);
            if ploteachtrial
            plot(thesepeaks-(thisstart)/1e4,ones(size(thesepeaks))+i,'.k','markersize',12)
            hold on
            plot(X-((thisstart)/1e4),(filtered(thisstart-5e4:thisstart+5e4)./5)+i+1,'color',[0.5 0.5 0.5])
            end
            ControlBreathing.Breathingtimes{Sessioncounter} =thesepeaks-(thisstart)/1e4;
            ControlBreathing.Breathingsignal(Sessioncounter,:)=filtered(thisstart-5e4:thisstart+5e4);
            ControlBreathing.Breathingsignalraw(Sessioncounter,:)=Data(Sessioncounter).thermistor2(thisstart-5e4:thisstart+5e4);
            ControlBreathing.whiskerraw(Sessioncounter,:)=Data(Sessioncounter).whisker(thisstart-5e4:thisstart+5e4);
        end

    end

end
%%
figure
for i=1:55


    plot(AllBreathing.Breathingtimes{i},ones(size(AllBreathing.Breathingtimes{i}))+i,'.k')
    hold on
end



%%
figure
for i=1:55
    % plot(AllBreathing.whiskerraw(i,:)-nanmean(AllBreathing.whiskerraw(i,1:5e4)))
    hold on
    thisfilter=AllBreathing.whiskerraw(i,:)-nanmean(AllBreathing.whiskerraw(i,4e4:5e4));
    % TF = isoutlier(thisfilter);
    thisfilter(thisfilter>2.5)=NaN;
    thisfilter = medfilt1((thisfilter),1e3);
    plot(thisfilter)
        plot(AllBreathing.whiskerraw(i,:))

    AllBreathing.whiskerfilter(i,:)=thisfilter;
end


%%

figure
for i=1:55
    % plot(AllBreathing.whiskerraw(i,:)-nanmean(AllBreathing.whiskerraw(i,1:5e4)))
    hold on
    thisfilter=AllBreathing.whiskerraw(i,:)-nanmean(AllBreathing.whiskerraw(i,4e4:5e4));
    firstbreath=find(AllBreathing.Breathingtimes{i}>0,1,'first');
    firsttime=round(AllBreathing.Breathingtimes{i}(firstbreath)*1e4+5e4);
    %
    %               Airfirsttime= round(AllBreathing.Breathingtimes{i}(firstbreath)*1e4+5e4);  ControlBreathing.Breathingtimes{Sessioncounter} =thesepeaks-(thisstart)/1e4;

    lastbreathi=find(AllBreathing.Breathingtimes{i}<0,1,'last');

    lastbreathi2=find(AllBreathing.Breathingtimes{i}<0,2,'last');


    Lastbefore=round(AllBreathing.Breathingtimes{i}(lastbreathi)*1e4+5e4);
    if any(lastbreathi2)
        if numel(lastbreathi2)>1
            IBI1(i)=diff(round(AllBreathing.Breathingtimes{i}(lastbreathi2)*1e4+5e4));
        end
    end
    if firstbreath>1
        IBI2(i)=diff(round(AllBreathing.Breathingtimes{i}(firstbreath-1:firstbreath)*1e4+5e4));
    end
    if firstbreath<numel(AllBreathing.Breathingtimes{i})
        IBI3(i)=diff(round(AllBreathing.Breathingtimes{i}(firstbreath:firstbreath+1)*1e4+5e4));
    end
    if firstbreath<numel(AllBreathing.Breathingtimes{i})-1
        IBI4(i)=diff(round(AllBreathing.Breathingtimes{i}(firstbreath+1:firstbreath+2)*1e4+5e4));
    end

    if isempty(firsttime)
        continue
    end

    Breathtimeafter(i,:)= (firsttime-5e4)./1e4;

    if isempty(Lastbefore)
        Breathtimebefore(i,:)=NaN;
        Moviestarindbefore(i,:)= NaN;

    else
        Breathtimebefore(i,:)= (Lastbefore-5e4)./1e4;
        Moviestarindbefore(i,:)= round((Lastbefore-5e4)./1e4*200);

    end
    % TF = isoutlier(thisfilter);
%     thisfilter(thisfilter>2.5)=NaN;
%     thisfilter = medfilt1((thisfilter),1e3);
%     plot(thisfilter)
%     AllBreathing.whiskerfilter(i,:)=thisfilter;
%     Breathaliged(i,:) = thisfilter(firsttime-4e4:firsttime+2e4);
    Moviestarind(i,:)= round((firsttime-5e4)./1e4*200+200);

end

%%
figure
bins = [0:0.1:3]
thishist1=hist(Breathtimeafter,bins)
thishist2=hist(-Breathtimebefore,bins)
h=bar(bins,[thishist1 ;thishist2]');
xlim([-0.1 1.5])
figurecorrectionOPTO(gca,[0.05, 0.05],12)

%% plot nose as similar to Martins

AllMoviemeasures.noseX(AllMoviemeasures.noseX==0)=NaN;
AllMoviemeasures.noseY(AllMoviemeasures.noseY==0)=NaN;
badtrials = nanmean(AllMoviemeasures.eyenorm_opening')==0;
AllMoviemeasures.eyenorm_opening(badtrials,:)=NaN;
for i=1:55

    AllMoviemeasures.noseXnorm(i,:)=AllMoviemeasures.noseX(i,:)-nanmedian(AllMoviemeasures.noseX(i,1:180));
    AllMoviemeasures.noseYnorm(i,:)=AllMoviemeasures.noseY(i,:)-nanmedian(AllMoviemeasures.noseY(i,1:180));
    AllMoviemeasures.CnoseXnorm(i,:)=AllMoviemeasures.ControlnoseX(i,:)-nanmedian(AllMoviemeasures.ControlnoseX(i,1:180));
    AllMoviemeasures.CnoseYnorm(i,:)=AllMoviemeasures.ControlnoseY(i,:)-nanmedian(AllMoviemeasures.ControlnoseY(i,1:180));


end
subplot(4,1,1:2)
hold off
B = localcontrast(x26september_rat335_20220926_131556_20220926_155929);
imagesc(B*1.1)

hold on

plot(AllMoviemeasures.noseXnorm(:,1:180)+mean(nanmedian(AllMoviemeasures.noseX(:,1:180))),AllMoviemeasures.noseYnorm(:,1:180)+mean(nanmedian(AllMoviemeasures.noseY(:,1:180)))-5,'.y')

hold on
plot(AllMoviemeasures.noseXnorm(:,210:225)+mean(nanmedian(AllMoviemeasures.noseX(:,1:180))),AllMoviemeasures.noseYnorm(:,210:225)+mean(nanmedian(AllMoviemeasures.noseY(:,1:180)))-5,'.r')
axis image

% plot(Aaxis image
xlim([180 450])
ylim([220 400])


% llMoviemeasures.noseXnorm(:,210:250),AllMoviemeasures.noseYnorm(:,210:250),'.r')
subplot(4,1,3)

hold off

% yyaxis left

plot(nanmean(AllMoviemeasures.CnoseXnorm))
hold on

plot(nanmean(AllMoviemeasures.noseXnorm))
% yyaxis right




subplot(4,1,4)
hold off
% yyaxis right
plot(-nanmean(AllMoviemeasures.CnoseYnorm))
hold on

plot(-nanmean(AllMoviemeasures.noseYnorm))





%%



D = sqrt([AllMoviemeasures.noseXnorm.^2+ AllMoviemeasures.noseYnorm.^2]);

[~,id]=max(D(:,200:400)')


[RHO,PVAL] = corr(id(~isnan(Moviestarind))',Moviestarind(~isnan(Moviestarind)),'Type','Spearman');
Moviestarind(Moviestarind==0)=NaN;



%%

figure
hist(IBI2(IBI2>0)./1e4,50)

IBI2(IBI2==0)=NaN;
IBI1(IBI1==0)=NaN;
IBI3(IBI3==0)=NaN;

x=[IBI1 ;IBI2; IBI3]
boxplot(x')

for i=1:55
    % plot(AllBreathing.whiskerraw(i,:)-nanmean(AllBreathing.whiskerraw(i,1:5e4)))
    hold on
    thisfilter=ControlBreathing.whiskerraw(i,:)-nanmean(ControlBreathing.whiskerraw(i,4e4:5e4));
    firstbreath=find(ControlBreathing.Breathingtimes{i}>0,1,'first');
    firsttime=round(ControlBreathing.Breathingtimes{i}(firstbreath)*1e4+5e4);
    lastbreathi=find(ControlBreathing.Breathingtimes{i}<0,1,'last');

    lastbreathi2=find(ControlBreathing.Breathingtimes{i}<0,2,'last');
    if isempty(firsttime)
            ControlBreathtimeafter(i)= NaN;

        continue
    end

    Lastbefore=round(ControlBreathing.Breathingtimes{i}(lastbreathi)*1e4+5e4);
    if any(lastbreathi2)
        if numel(lastbreathi2)>1
            ControlIBI1(i)=diff(round(ControlBreathing.Breathingtimes{i}(lastbreathi2)*1e4+5e4));
        end
    end
    if firstbreath>1
        ControlIBI2(i)=diff(round(ControlBreathing.Breathingtimes{i}(firstbreath-1:firstbreath)*1e4+5e4));
    end
    if firstbreath<numel(ControlBreathing.Breathingtimes{i})
        ControlIBI3(i)=diff(round(ControlBreathing.Breathingtimes{i}(firstbreath:firstbreath+1)*1e4+5e4));
    end

    if firstbreath<numel(ControlBreathing.Breathingtimes{i})-1
        ControlIBI4(i)=diff(round(ControlBreathing.Breathingtimes{i}(firstbreath+1:firstbreath+2)*1e4+5e4));
    end



    ControlBreathtimeafter(i)= (firsttime-5e4)./1e4;

    if isempty(Lastbefore)
        ControlBreathtimebefore(i,:)=NaN;
        ControlMoviestarindbefore(i,:)= NaN;

    else
        ControlBreathtimebefore(i,:)= (Lastbefore-5e4)./1e4;
        ControlMoviestarindbefore(i,:)= round((Lastbefore-5e4)./1e4*200);

    end
    airMoviestarind(i,:)= round((firsttime-5e4)./1e4*200+200);

end

figure
ControlIBI4(ControlIBI4==0)=NaN;

ControlIBI3(ControlIBI3==0)=NaN;
ControlIBI1(ControlIBI1==0)=NaN;
ControlIBI2(ControlIBI2==0)=NaN;

% x=[ControlIBI1 ;ControlIBI2; ControlIBI3]
% boxplot(x')
%%
%% plot realigend
 
clearvars thisrealingedx thisrealingedy airthisrealingedx airthisrealingedy thisrealingedx
for i=1:55
    if Moviestarind(i)<200||isnan(Moviestarind(i))

 airthisrealingedx.Nose(i,:) = NaN(1,401);
        airthisrealingedy.Nose(i,:) = NaN(1,401);
        airthisrealingedx.Eye(i,:) = NaN(1,401);
         thisrealingedy.Nose(i,:) = NaN(1,401);
        thisrealingedx.Nose(i,:) = NaN(1,401);
        thisrealingedx.Eye(i,:) = NaN(1,401);

        continue
    end
    Y = [AllMoviemeasures.noseXnorm(i,:) NaN(1,200)];
    thisrealingedx.Nose(i,:)=  Y(Moviestarind(i)-200:Moviestarind(i)+200);
    Y = [AllMoviemeasures.noseYnorm(i,:) NaN(1,200)];
    thisrealingedy.Nose(i,:)=  Y(Moviestarind(i)-200:Moviestarind(i)+200);
     Y = [AllMoviemeasures.eyenorm_opening(i,:) NaN(1,200)];
    thisrealingedx.Eye (i,:)=  Y(Moviestarind(i)-200:Moviestarind(i)+200);

%     thisrealingedy(i,:)= AllMoviemeasures.noseYnorm(i,Moviestarind-100:Moviestarind+100);


    if airMoviestarind(i)>=200 && airMoviestarind(i)<400

        Y = [AllMoviemeasures.CnoseXnorm(i,:) NaN(1,200)];
        airthisrealingedx.Nose(i,:)= Y(airMoviestarind(i)-200:airMoviestarind(i)+200); 
                Y = [AllMoviemeasures.CnoseYnorm(i,:) NaN(1,200)];
        airthisrealingedy.Nose(i,:)= Y(airMoviestarind(i)-200:airMoviestarind(i)+200); 
        Y = [AllMoviemeasures.Controleyenorm_opening(i,:) NaN(1,200)];
        airthisrealingedx.Eye(i,:)= Y(airMoviestarind(i)-200:airMoviestarind(i)+200);


    else
        airthisrealingedx.Nose(i,:) = NaN(1,401);
        airthisrealingedy.Nose(i,:) = NaN(1,401);
        airthisrealingedx.Eye(i,:) = NaN(1,401);


    end

end





%%  replot Nose X 
figure
subplot(2,1,1)

x=(1:401)/200-1;
y = nanmean(AllMoviemeasures.CnoseXnorm);
errBar = nanstd(AllMoviemeasures.CnoseXnorm)./sqrt(size(AllMoviemeasures.CnoseXnorm,1));
h2 = shadedErrorBar(x,y,errBar,'lineprops',{'b-','markerfacecolor','b'})
hold on
y = nanmean(AllMoviemeasures.noseXnorm);
errBar = nanstd(AllMoviemeasures.noseXnorm)./sqrt(size(AllMoviemeasures.noseXnorm,1));
h2 = shadedErrorBar(x,y,errBar,'lineprops',{'r-','markerfacecolor','r'})

xlim([-1 1])
subplot(2,1,2)
% x=(1:201)/200-0.5;
y = nanmean(airthisrealingedx.Nose);
errBar = nanstd(airthisrealingedx.Nose)./sqrt(size(airthisrealingedx.Nose,1));
h2 = shadedErrorBar(x,y,errBar,'lineprops',{'b-','markerfacecolor','b'})
hold on
y = nanmean(thisrealingedx.Nose);
errBar = nanstd(thisrealingedx.Nose)./sqrt(size(thisrealingedx.Nose,1));
h2 = shadedErrorBar(x,y,errBar,'lineprops',{'r-','markerfacecolor','r'})
xlim([-1 1])


%%  Nose Y 
figure
subplot(2,1,1)

x=(1:401)/200-1;
y = nanmean(-AllMoviemeasures.CnoseYnorm);
errBar = nanstd(-AllMoviemeasures.CnoseYnorm)./sqrt(size(AllMoviemeasures.CnoseYnorm,1));
% Mean.Contorl.
h2 = shadedErrorBar(x,y,errBar,'lineprops',{'b-','markerfacecolor','b'})
hold on
y = nanmean(-AllMoviemeasures.noseYnorm);
errBar = nanstd(-AllMoviemeasures.noseYnorm)./sqrt(size(AllMoviemeasures.noseYnorm,1));
h2 = shadedErrorBar(x,y,errBar,'lineprops',{'r-','markerfacecolor','r'})

xlim([-1 1])
subplot(2,1,2)
% x=(1:201)/200-0.5;
y = nanmean(-airthisrealingedy.Nose);
errBar = nanstd(-airthisrealingedy.Nose)./sqrt(size(airthisrealingedy.Nose,1));
h2 = shadedErrorBar(x,y,errBar,'lineprops',{'b-','markerfacecolor','b'})
hold on
y = nanmean(-thisrealingedy.Nose);
errBar = nanstd(-thisrealingedy.Nose)./sqrt(size(thisrealingedy.Nose,1));
h2 = shadedErrorBar(x,y,errBar,'lineprops',{'r-','markerfacecolor','r'})
xlim([-1 1])

%%  Eye opening
figure
subplot(2,1,1)

x=(1:401)/200-1;
y = nanmean(AllMoviemeasures.Controleyenorm_opening);
errBar = nanstd(AllMoviemeasures.Controleyenorm_opening)./sqrt(size(AllMoviemeasures.Controleyenorm_opening,1));
h2 = shadedErrorBar(x,y,errBar,'lineprops',{'b-','markerfacecolor','b'})
hold on
y = nanmean(AllMoviemeasures.eyenorm_opening);
errBar = nanstd(AllMoviemeasures.eyenorm_opening)./sqrt(size(AllMoviemeasures.eyenorm_opening,1));
h2 = shadedErrorBar(x,y,errBar,'lineprops',{'r-','markerfacecolor','r'})

xlim([-1 1])
subplot(2,1,2)
% x=(1:201)/200-0.5;
y = nanmean(airthisrealingedx.Eye );
y(370:401)=NaN;

errBar = nanstd(airthisrealingedx.Eye)./sqrt(size(airthisrealingedx.Eye,1));
h2 = shadedErrorBar(x,y,errBar,'lineprops',{'b-','markerfacecolor','b'})
hold on
y = nanmean(thisrealingedx.Eye);
y(370:401)=NaN;
errBar = nanstd(thisrealingedx.Eye)./sqrt(size(thisrealingedx.Eye,1));
h2 = shadedErrorBar(x,y,errBar,'lineprops',{'r-','markerfacecolor','r'})
xlim([-1 1])

%%  viberissa opening
figure
% Breathtimeafter
% ControlBreathtimeafter
 airthisrealingedx.Whisker=NaN(55,80001);
 thisrealingedx.Whisker=NaN(55,80001);
 ControlBreathing.whiskerraw(ControlBreathing.whiskerraw==0)=NaN;
  ControlBreathing.whiskerraw(ControlBreathing.whiskerraw>7)=NaN;
  AllBreathing.whiskerraw(AllBreathing.whiskerraw>7)=NaN;

  AllBreathing.whiskerraw(AllBreathing.whiskerraw==0)=NaN;

for i=1:55
 thisrealingedx.Whisker(i,:)= medfilt1(AllBreathing.whiskerraw(i,Breathtimeafter(i)*1e4+5e4-4e4:Breathtimeafter(i)*1e4+5e4+4e4),2e3);
  
 AllBreathing.whiskerfilter(i,:)= medfilt1(AllBreathing.whiskerraw(i,:),2e3);
  ControlBreathing.whiskerfilter(i,:)= medfilt1(ControlBreathing.whiskerraw(i,:),2e3);

%     thisrealingedy.Nose(i,:)= AllMoviemeasures.whiskerfilter(i,Moviestarind(i)-100:Moviestarind(i)+100);
if ControlBreathtimeafter(i)~=0&&~isnan(ControlBreathtimeafter(i))
 airthisrealingedx.Whisker(i,:)=medfilt1(ControlBreathing.whiskerraw(i,ControlBreathtimeafter(i)*1e4+5e4-4e4:ControlBreathtimeafter(i)*1e4+5e4+4e4),2e3);
end
end
figure
subplot(2,1,1)
x=(1:100001)/1e4-5;
thesedownsample = round(linspace(1,numel(x),1000))

y = medfilt1(nanmean(-AllBreathing.whiskerfilter),1e3)./0.35*40;
errBar = nanstd(AllBreathing.whiskerfilter./0.35*40)./sqrt(size(AllBreathing.whiskerraw,1));
h2 = shadedErrorBar(x(thesedownsample),y(thesedownsample)+500,errBar(thesedownsample),'lineprops',{'r-','markerfacecolor','r'})
hold on
y = medfilt1(nanmean(-ControlBreathing.whiskerfilter./0.35*40),1e3);
errBar = nanstd(ControlBreathing.whiskerfilter./0.35*40)./sqrt(size(ControlBreathing.whiskerraw,1));
h2 = shadedErrorBar(x(thesedownsample),y(thesedownsample)+500,errBar(thesedownsample),'lineprops',{'b-','markerfacecolor','b'})

xlim([-1 1])
subplot(2,1,2)
x=(1:80001)/1e4-4;
thesedownsample = round(linspace(1,numel(x),1000))

y = medfilt1(nanmean(-airthisrealingedx.Whisker ./0.35*40),1e3);

errBar = nanstd(airthisrealingedx.Whisker./0.35*40)./sqrt(size(airthisrealingedx.Whisker,1));
h2 = shadedErrorBar(x(thesedownsample),y(thesedownsample)+500,errBar(thesedownsample),'lineprops',{'b-','markerfacecolor','b'})
hold on
y =medfilt1(nanmean(-thisrealingedx.Whisker./0.35*40),1e3);
errBar = nanstd(thisrealingedx.Whisker./0.35*40)./sqrt(size(thisrealingedx.Whisker,1));
h2 = shadedErrorBar(x(thesedownsample),y(thesedownsample)+500,errBar(thesedownsample),'lineprops',{'r-','markerfacecolor','r'})
xlim([-5 5])




%%
%%  viberissa opening normalized
figure
% Breathtimeafter
% ControlBreathtimeafter
 airthisrealingedx.Whisker=NaN(55,80001);
 thisrealingedx.Whisker=NaN(55,80001);
 ControlBreathing.whiskerraw(ControlBreathing.whiskerraw==0)=NaN;
  ControlBreathing.whiskerraw(ControlBreathing.whiskerraw>7)=NaN;
  AllBreathing.whiskerraw(AllBreathing.whiskerraw>7)=NaN;

  AllBreathing.whiskerraw(AllBreathing.whiskerraw==0)=NaN;

for i=1:55
 thisrealingedx.Whisker(i,:)= medfilt1(AllBreathing.whiskerraw(i,Breathtimeafter(i)*1e4+5e4-4e4:Breathtimeafter(i)*1e4+5e4+4e4),2e3);
   AllBreathing.whiskerfilter(i,:)= medfilt1(AllBreathing.whiskerraw(i,:),2e3);
  ControlBreathing.whiskerfilter(i,:)= medfilt1(ControlBreathing.whiskerraw(i,:),2e3);

 thisrealingedx.Whisker(i,:)= thisrealingedx.Whisker(i,:)-nanmean(thisrealingedx.Whisker(i,1e4:2e4));
   AllBreathing.whiskerfilter(i,:)= medfilt1(AllBreathing.whiskerfilter(i,:),2e3)-nanmean(AllBreathing.whiskerfilter(i,1e4:2e4));
  ControlBreathing.whiskerfilter(i,:)= medfilt1(ControlBreathing.whiskerfilter(i,:),2e3)-nanmean(ControlBreathing.whiskerfilter(i,1e4:2e4));
if ControlBreathtimeafter(i)~=0&&~isnan(ControlBreathtimeafter(i))


 airthisrealingedx.Whisker(i,:)=medfilt1(ControlBreathing.whiskerraw(i,ControlBreathtimeafter(i)*1e4+5e4-4e4:ControlBreathtimeafter(i)*1e4+5e4+4e4),2e3);
 airthisrealingedx.Whisker(i,:)=  airthisrealingedx.Whisker(i,:)-nanmean( airthisrealingedx.Whisker(i,1e4:2e4));

end
end
figure
Ratio =1./0.35*40;
Ratio = 1;
subplot(2,1,1)
x=(1:100001)/1e4-5;
thesedownsample = round(linspace(1,numel(x),1000))

y = medfilt1(nanmean(-AllBreathing.whiskerfilter),1e3)*Ratio;
errBar = nanstd(AllBreathing.whiskerfilter*Ratio)./sqrt(size(AllBreathing.whiskerraw,1));
h2 = shadedErrorBar(x(thesedownsample),y(thesedownsample),errBar(thesedownsample),'lineprops',{'r-','markerfacecolor','r'})
hold on
y = medfilt1(nanmean(-ControlBreathing.whiskerfilter*Ratio),1e3);
errBar = nanstd(ControlBreathing.whiskerfilter*Ratio)./sqrt(size(ControlBreathing.whiskerraw,1));
h2 = shadedErrorBar(x(thesedownsample),y(thesedownsample),errBar(thesedownsample),'lineprops',{'b-','markerfacecolor','b'})

xlim([-1 1])
subplot(2,1,2)
x=(1:80001)/1e4-4;
thesedownsample = round(linspace(1,numel(x),1000))

y = medfilt1(nanmean(-airthisrealingedx.Whisker*Ratio),1e3);

errBar = nanstd(airthisrealingedx.Whisker*Ratio)./sqrt(size(airthisrealingedx.Whisker,1));
h2 = shadedErrorBar(x(thesedownsample),y(thesedownsample),errBar(thesedownsample),'lineprops',{'b-','markerfacecolor','b'})
hold on
y =medfilt1(nanmean(-thisrealingedx.Whisker*Ratio),1e3);
errBar = nanstd(thisrealingedx.Whisker*Ratio)./sqrt(size(thisrealingedx.Whisker,1));
h2 = shadedErrorBar(x(thesedownsample),y(thesedownsample),errBar(thesedownsample),'lineprops',{'r-','markerfacecolor','r'})
xlim([-5 5])



%%
thishist1 = hist(ControlBreathtimeafter,[0:0.09:3]);
thishist2 = hist(-ControlBreathtimebefore,[0:0.09:3]);

figure

h = bar([0:0.09:3],[thishist1 ;thishist2]');
%%
figure
tbl.IBI=[ControlIBI1 IBI1 ControlIBI2 IBI2 ControlIBI3 IBI3  ControlIBI4 IBI4];
tbl.Order = [ones(size(ControlIBI1)) ones(size(ControlIBI1)) ones(size(ControlIBI1))*3 ones(size(ControlIBI1))*3 ones(size(ControlIBI1))*6 ones(size(ControlIBI1))*6  ones(size(ControlIBI1))*9 ones(size(ControlIBI1))*9];
tbl.ControlvsStim=[ones(size(ControlIBI1)) ones(size(ControlIBI1))*2 ones(size(ControlIBI1)) ones(size(ControlIBI1))*2 ones(size(ControlIBI1)) ones(size(ControlIBI1))*2  ones(size(ControlIBI1)) ones(size(ControlIBI1))*2];



h = boxchart(tbl.Order,(tbl.IBI),'GroupByColor',tbl.ControlvsStim);

tbl2=struct2table(tbl);
% tblstats1 = grpstats(tbl,"ControlvsStim")

%%  Eye opening

clearvars thisrealingedx thisrealingedy airthisrealingedx airthisrealingedy thisrealingedx
for i=1:55
    if Moviestarind(i)<200||isnan(Moviestarind(i))

 airthisrealingedx.Nose(i,:) = NaN(1,401);
        airthisrealingedy.Nose(i,:) = NaN(1,401);
        airthisrealingedx.Eye(i,:) = NaN(1,401);
         thisrealingedy.Nose(i,:) = NaN(1,401);
        thisrealingedx.Nose(i,:) = NaN(1,401);
        thisrealingedx.Eye(i,:) = NaN(1,401);

        continue
    end
    Y = [AllMoviemeasures.noseXnorm(i,:) NaN(1,200)];
    thisrealingedx.Nose(i,:)=  Y(Moviestarind(i)-200:Moviestarind(i)+200);
    Y = [AllMoviemeasures.noseYnorm(i,:) NaN(1,200)];
    thisrealingedy.Nose(i,:)=  Y(Moviestarind(i)-200:Moviestarind(i)+200);
     Y = [AllMoviemeasures.eyenorm_opening(i,:) NaN(1,200)];
    thisrealingedx.Eye (i,:)=  Y(Moviestarind(i)-200:Moviestarind(i)+200);

%     thisrealingedy(i,:)= AllMoviemeasures.noseYnorm(i,Moviestarind-100:Moviestarind+100);


    if airMoviestarind(i)>=200 && airMoviestarind(i)<400

        Y = [AllMoviemeasures.CnoseXnorm(i,:) NaN(1,200)];
        airthisrealingedx.Nose(i,:)= Y(airMoviestarind(i)-200:airMoviestarind(i)+200); 
                Y = [AllMoviemeasures.CnoseYnorm(i,:) NaN(1,200)];
        airthisrealingedy.Nose(i,:)= Y(airMoviestarind(i)-200:airMoviestarind(i)+200); 
        Y = [AllMoviemeasures.Controleyenorm_opening(i,:) NaN(1,200)];
        airthisrealingedx.Eye(i,:)= Y(airMoviestarind(i)-200:airMoviestarind(i)+200);


    else
        airthisrealingedx.Nose(i,:) = NaN(1,401);
        airthisrealingedy.Nose(i,:) = NaN(1,401);
        airthisrealingedx.Eye(i,:) = NaN(1,401);


    end

end

figure
subplot(2,1,1)

x=(1:401)/200-1;
y = nanmean(AllMoviemeasures.Controleyenorm_opening);
errBar = nanstd(AllMoviemeasures.Controleyenorm_opening)./sqrt(size(AllMoviemeasures.Controleyenorm_opening,1));
h2 = shadedErrorBar(x,y,errBar,'lineprops',{'b-','markerfacecolor','b'})
hold on
y = nanmean(AllMoviemeasures.eyenorm_opening);
errBar = nanstd(AllMoviemeasures.eyenorm_opening)./sqrt(size(AllMoviemeasures.eyenorm_opening,1));
h2 = shadedErrorBar(x,y,errBar,'lineprops',{'r-','markerfacecolor','r'})

xlim([-1 1])
subplot(2,1,2)
% x=(1:201)/200-0.5;
y = nanmean(airthisrealingedx.Eye );
y(370:401)=NaN;

errBar = nanstd(airthisrealingedx.Eye)./sqrt(size(airthisrealingedx.Eye,1));
h2 = shadedErrorBar(x,y,errBar,'lineprops',{'b-','markerfacecolor','b'})
hold on
y = nanmean(thisrealingedx.Eye);
y(370:401)=NaN;
errBar = nanstd(thisrealingedx.Eye)./sqrt(size(thisrealingedx.Eye,1));
h2 = shadedErrorBar(x,y,errBar,'lineprops',{'r-','markerfacecolor','r'})
xlim([-1 1])

