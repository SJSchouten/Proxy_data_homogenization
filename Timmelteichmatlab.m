clear all
close all

% --- Defining the cores ---
%% 

dat = xlsread('Excel Timmelteich.xlsx','XRF');
dat = sortrows(dat,3);

names = {'K','Ca','Ti','Mn','Fe','Rb','Sr','Zr','Ba'}
indexing = [11,13,15,19,21,29,31,33,39]

loi_2020 = xlsread('Excel Timmelteich.xlsx','TT-20 Stan');
loi_2011 = xlsread('Excel Timmelteich.xlsx','TT-11 Renee');
loi_2013 = sortrows(xlsread('Excel Timmelteich.xlsx','TT-13 Renee'),6);
gisp = xlsread('GISPII_oxygen_isotopes.xls');
icecores = xlsread('groenland comb');
ieccores = sortrows(icecores,1)

figure, subplot(4,1,1),stairs(loi_2020(:,6),loi_2020(:,16));
hold on, stairs(loi_2011(:,6),loi_2011(:,16));
hold on, stairs(loi_2013(:,6)-13,loi_2013(:,16));
xlabel('depth'); ylabel('LOI(%)'); legend('2020','2011','2013'); grid; xlim([365,470]);

subplot(4,1,2), stairs(loi_2020(:,6),(loi_2020(:,14)./loi_2020(:,1)))
hold on, yyaxis right,stairs(loi_2020(:,6),(loi_2020(:,13)./loi_2020(:,1)))
legend('Klastische dichtheid', 'Organische dichtheid');xlim([365,470]);title('2020 resampling')

subplot(4,1,3), stairs(loi_2011(:,6),(loi_2011(:,14)./loi_2011(:,1)))
hold on, yyaxis right,stairs(loi_2011(:,6),(loi_2011(:,13)./loi_2011(:,1)))
title('2011 core');xlim([365,470]);

subplot(4,1,4),stairs(loi_2013(:,6)-13,(loi_2013(:,14)./loi_2013(:,1)));
hold on, yyaxis right,stairs(loi_2013(:,6)-13,(loi_2013(:,13)./loi_2013(:,1)));title('2013 core');xlim([365,470]);

% The amount of carbon

figure,
for i = 1:length(names)
    subplot(length(names)+1,1,i),plot(dat(:,3),dat(:,(indexing(i)))*100,'x-k');grid;xlabel('Depth (cm)');ylabel('ppm')
%     hold on, scatter(data_xrf(:,3),(data_xrf(:,6+i)-data_xrf(7+i)),10,'--r')
%     hold on, plot(data_xrf(:,3),(data_xrf(:,6+i)+data_xrf(7+i)),10,'--r')
    xlim([min(dat(:,3)),max(dat(:,3))]);
    title(names{i});
end

%% Calling the age depth model

inputfilename = 'Excel Timmelteich.xlsx'
inputfiletabname = 'TT-11 composite'
exponent = [0.3:0.1:2.5]

depthfix = [506.5,469,423,131]
timefix = [14642,12846,8500,1150]
fixnames = {'GI-1 start','GS-1 start','Alnus','Guess core top'}

tiepnts_polx = [441.9,441.9,433.4,433.4];
tiepnts_poly = [10300,10800,9700,10300];
polnames = {'Corylus','Quercus & Ulmus'}
choice = 0.9
method = 2

[out, fig, f] = Agedepthmodelling(inputfilename,inputfiletabname,exponent,timefix,depthfix,fixnames,tiepnts_polx,tiepnts_poly,polnames,choice);

%% --- Looping detrend and loading operations over all images ---

cores = {'1'};
c_depths = sort([370,470] ,'descend');

s = 1;

for i = 1:length(cores)
    I{i} = imread('TT_Finals.png');
    
    d_scale = [c_depths((i*2)-1) c_depths((i*2))];
    Image = I{i};

    % figure,
    Im = rgb2gray(Image);
    % imshow(Im),improfile

    sz = size(Im);
    som = zeros(sz(1),1);

    % Finding the scanline domain in the vertical

    for j = 1:8
        xh = [(j/8)*sz(2),(j/8)*sz(2)];
        yh = [1,sz(1)];
        hor = improfile(Im,xh,yh);
        som = som + hor(:,:,1);
    end
    
    som = som/8;
    % figure,plot(som)
    ind = find(((diff(som)) > 50)==1);
    ind = ind(ind > 100 & ind < (length(som)-100));

    ranges = [ind(find(diff(ind)>1)) ind(find(diff(ind)>1) +1)];

    lines = ((ranges(1)+100):100:(ranges(2)-100));

    % --- Making scanlines ---
    %%

    sommetje = zeros(sz(2),1);
    
    for z = 1:length(lines)
        y = [lines(z),lines(z)];
        x = [1,sz(2)];
        n = improfile(Im,x,y);
        vect = n(:,:,1);
        sommetje = sommetje + vect;
    end
    
    vec = flip(sommetje/length(lines));

    % Manual decracking
    
%     if i == 1
%     delpos = [6009:1:6770];    
%     vec(delpos) = [];
%     end
    
    stretch = length(vec);
    ims{i} = [d_scale(1):(d_scale(1)-d_scale(2))/-stretch:d_scale(2)];
    
    inc = find(vec > 252); 
    bots = inc(find(diff(inc) > 1));
    tops = inc(find(diff(inc) > 1)+1);

%     % General detrending of holocene cores
%     
%     if i >= 3
%         vec = detrend(vec) + 100;
%     elseif i == 1
%         vec = vec+5;
%     end

%      if i == 2
%          vec = vec-65;
%      end

    %% --- Scanline average block detrending ---
 
    for comp =1:15 
        for n = 1:length(bots)-1 
            differ{i}(comp,n) = mean(vec((tops(n)-(10*comp)):(tops(n)-5)))-mean(vec((bots(n+1)+5):(bots(n+1)+(10*comp)))); 
        end
    end

    richtingcoff = mean(mean(differ{i}(:,(11:13)),1))/mean(tops-bots);
    blocks{i} = (tops-bots)*richtingcoff;

    for j = 1:length(blocks{i})
%         detvec{j} = vec(bots(j):tops(j))'-(-(blocks{i}(j)/2):richtingcoff:(blocks{i}(j)/2));
        detvec{j} = vec(bots(j):tops(j))'
        xvec{j} = (bots(j):1:tops(j));
    end
    
    figure,
    for p = 1:length(detvec)
        subplot(4,1,1),plot(xvec{p},detvec{p},':r')
        hold on, plot(xvec{p},movmean(detvec{p},60),'k');grid;set(gca, 'Xdir', 'reverse');set(gca, 'Ydir', 'reverse');
        
        allvec{p,i} = detvec{p};

        subplot(4,1,2),hold on, findpeaks(movmean(detvec{p},30),xvec{p},'MinPeakDistance',5,'MinPeakProminence',15,'Annotate','extents');grid
        [pks{i},locs{i}] =  findpeaks(movmean(detvec{p},30),xvec{p},'MinPeakDistance',5,'MinPeakProminence',15);
        text(locs{i}',pks{i}+10',num2str((s:(s-1)+(numel(pks{i})))'),'FontSize',8,'HorizontalAlignment','center')

        s = s+(numel(pks{i}));

        allvec{p,i} = detvec{p};
        allpks{p,i} = pks{i};
        allpkslocs{p,i} = locs{i};
        alldpkslocs{p,i} = []
        
    end
      
    grid;set(gca, 'Ydir', 'reverse');title('core '+string(cores{i})+' Greyscale scan'); set(gca, 'Xdir', 'reverse')
    subplot(4,1,3),plot(ims{i}(1:length(vec)),vec);grid;set(gca, 'Ydir', 'reverse');xlim([min(d_scale) max(d_scale)]) 
    subplot(4,1,4),imshow(I{i}); title('core '+string(cores{i}));
    
    o(1) = 1;
    for y = 2:length(tops)+1
        o(y) = length(bots(y-1):tops(y-1));
    end
    coff = (d_scale(1)-d_scale(2))/sum(o);
    depth_scale{i} = flip((d_scale(2)+coff):coff:d_scale(1));
    for y = 1:length(tops)
        alldvec{y,i} = depth_scale{i}(o(y):(o(y)+o(y+1))-1);
        if length(allpkslocs{y,i}) > 0
            for b = 1:length(allpkslocs{y,i})
                intermezzo(b) = alldvec{y,i}(find(xvec{y} == allpkslocs{y,i}(b)));
            end
                alldpkslocs{y,i} = intermezzo;
                clear intermezzo;
        end
        
        o(y+1) = o(y)+o(y+1);
    end
    
    clear o;    
    clear xvec;
    clear detvec;
    clear o;
    clear delposit;
    clear vec;
    clear xvect;
    clear count;
end

s = 1

    counter = 0
    
    for h = 1:length(alldvec)
        counter = counter + length(alldvec{h}); 
    end
    counter/(c_depths(1)-c_depths(2));
    temp = [find(out.Depth == [c_depths(1)-1]),find(out.Depth == [c_depths(2)+5])]
    ageselect = out.Age(temp(1):temp(2))
    depthselect = out.Depth(temp(1):temp(2))
    d_deler = c_depths(1)-1 - (c_depths(2)+5)

    [age_res,depth_res] = resample(ageselect,depthselect,(counter)/d_deler,1,1)
    allagevec = flip(age_res(1:length(age_res)-1))
    
    tel = 1
    te = 1
    
    for w = 1:length(allvec)
        longvector(tel:(tel-1+length(allvec{w}))) = allvec{w}
        longdepthvector(tel:(tel-1+length(alldvec{w}))) = alldvec{w}
        pksdlocs(te:(te-1+length(alldpkslocs{w}))) = alldpkslocs{w}
        allpeaks(te:(te-1+length(allpks{w}))) = allpks{w}
        tel = tel+length(allvec{w})
        te = te+length(alldpkslocs{w})
    end
    
%     figure, plot(age_res,depth_res,'x-',ageselect,depthselect,'o')
%     hold on, yyaxis right, plot(flip(age_res(1:length(age_res)-1)),longvector)
%     
%     figure, plot(longdepthvector',flip(age_res((1:length(depth_res)-1))))

%% plotting the event strat on depth

figure,
for p = 1:(14)
    subplot(6,1,1),plot(alldvec{p},allvec{p},':r')
    hold on, plot(alldvec{p},movmean(allvec{p},60),'k')
    hold on, plot(alldvec{p},movmean(allvec{p},300),'b')
    hold on,
    text(alldpkslocs{p},allpks{p}+(160-allpks{p}),num2str((s:(s-1)+(numel(allpks{p})))'),'FontSize',8,'HorizontalAlignment','center')
    s = s+(numel(allpks{p}));
    
%     counter = 0
%     for h=1:length(alldvec); counter = counter + length(alldvec{h}); 
%     end 
%     subplot(6,1,2),plot(alldve
end

set(gca, 'Ydir', 'reverse'); grid; xlim([min(c_depths) max(c_depths)]);xlim([365,470]);  

subplot(6,1,3),stairs(out.Depth,out.LOI);

xlim([365,470]);

subplot(6,1,4),plot(dat(:,1),dat(:,13)*100,'x-b');grid;xlabel('Depth (mm)');ylabel('% of cps');
title(names{5}); xlim([365,470]);
hold on, yyaxis right,plot(dat(:,1),dat(:,15)*100,'.-k');grid;xlabel('Depth (mm)');ylabel('% of cps');
title(names{5}); xlim([365,470]);
legend('Ca','Ti')

subplot(6,1,5),plot(dat(:,1),dat(:,19)*100,'^-b');grid;xlabel('Depth (mm)');ylabel('% of cps')
title(names{5}); xlim([365,470]);
hold on, yyaxis right,plot(dat(:,1),dat(:,33)*100,'o-r');grid;xlabel('Depth (mm)');ylabel('% of cps')
title(names{5}); xlim([365,470]);
legend('Mn','Zr')

subplot(6,1,6),imshow(I{i}); title('Timmelteich 2011 V');

figure,
for i = 1:length(names)
    subplot(length(names),1,i),plot(dat(:,3),dat(:,(indexing(i)))*100,'x-k');grid;xlabel('Depth (cm)');ylabel('ppm')
%     hold on, scatter(data_xrf(:,3),(data_xrf(:,6+i)-data_xrf(7+i)),10,'--r')
%     hold on, plot(data_xrf(:,3),(data_xrf(:,6+i)+data_xrf(7+i)),10,'--r')
    xlim([min(dat(:,3)),max(dat(:,3))]);
    title(names{i});
end

%% plotting on age

s = 1

figure,

subplot(6,1,1),plot(allagevec,longvector,':r')
hold on, plot(allagevec,movmean(longvector,60),'k')
hold on, plot(allagevec,movmean(longvector,300),'b')
hold on,
text(pksdlocs,allpeaks,num2str((s:(s-1)+(numel(allpeaks)))'),'FontSize',8,'HorizontalAlignment','center')
s = s+(numel(allpeaks));

set(gca, 'Ydir', 'reverse'); grid; xlim([min(c_depths) max(c_depths)]);xlim([min(allagevec),max(allagevec)]);  

subplot(6,1,2),stairs(out.Age,out.LOI);
xlim([min(allagevec),max(allagevec)]);

subplot(6,1,3), stairs(gisp(:,1),gisp(:,2)+2,'b');
hold on, stairs(icecores(:,1)-50,icecores(:,10),'k');
% hold on, stairs(icecores(:,1),icecores(:,3),'g')
% stairs(icecores(:,1),icecores(:,5),'b')
% stairs(icecores(:,1),icecores(:,7),'r')
legend('GISPII','Ice core average','Ngrip','grip','dye3');
xlim([min(allagevec),max(allagevec)]);;grid;

subplot(6,1,4),plot(dat(:,13),dat(:,1)*100,'x-b');grid;xlabel('Depth (mm)');ylabel('% of cps');
title('Calcium'); ylim([365,470]);
legend('Ca')

subplot(6,1,5),plot(dat(:,33),dat(:,1)*100,'o-r');grid;xlabel('Depth (mm)');ylabel('% of cps')
title('Zirconium'); ylim([365,470]);


subplot(6,1,6),plot(ageselect,depthselect);xlim([min(allagevec),max(allagevec)]);;grid;


figure,
for i = 1:length(names)
    subplot(length(names),1,i),plot(dat(:,3),dat(:,(indexing(i)))*100,'x-k');grid;xlabel('Depth (cm)');ylabel('ppm')
%     hold on, scatter(data_xrf(:,3),(data_xrf(:,6+i)-data_xrf(7+i)),10,'--r')
%     hold on, plot(data_xrf(:,3),(data_xrf(:,6+i)+data_xrf(7+i)),10,'--r')
    xlim([min(dat(:,3)),max(dat(:,3))]);
    title(names{i});
end



