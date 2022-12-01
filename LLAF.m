clear
close all

% --- Defining the cores ---
%%

dat = xlsread('XRF LLAS14 BII.xlsx','data for graph');

names = {'Al/total cps','Si/total cps','S/total cps','K/total cps','Ca/total cps','Ti/total cps','Mn/total cps','Fe/total cps','Zn/total cps','Br/total cps','Rb/total cps','Sr/total cps','Zr/total cps','Si/Al','Ca/Sr','K/Ti','Fe/Ca','Rb/Zn','Ti/Zn','Ca/Al','Ti/Al'};

lla = xlsread('Llangorse South Stan.xls','LOI composite');
% lla(:,1) = (lla(:,2)+lla(:,3))/2

datas = xlsread('Excel Timmelteich.xlsx','Timmelteich_2020_LOI_composite');
gisp = xlsread('GISPII_oxygen_isotopes.xls');
icecores = xlsread('groenland comb');
ieccores = sortrows(icecores,1)

figure, subplot(3,1,1),stairs(lla(:,6),lla(:,19)); xlim([880,1670]); xlabel('Depth (cm)'); grid; ylabel('LOI(%)');title('Loss on ignition percentage')
subplot(3,1,3), stairs(lla(:,6),lla(:,16));xlabel('Depth (cm)'); ylabel('\rho_{organic}'); grid; xlim([880,1670]);
hold on ,yyaxis right, stairs(lla(:,6),lla(:,17));xlabel('Depth (cm)'); ylabel('\rho_{clastic}'); grid; xlim([880,1670]); grid; title('Individual density components');legend('\rho_{organic}','\rho_{clastic}')
subplot(3,1,2), stairs(lla(:,6),(lla(:,17)./lla(:,16))); xlabel('Depth(cm)'); ylabel('\rho_{clastic}/\rho_{organic}'); grid; xlim([880,1670]); set(gca ,'Ydir','reverse');title('\rho_{clastic} / \rho_{organic} ratio')

figure, subplot(3,1,1),stairs(datas(:,3),datas(:,14));
hold on, stairs(datas(:,22),datas(:,23));
hold on, stairs(datas(:,25),datas(:,26));
xlabel('depth'); ylabel('LOI(%)'); legend('2020','2011','2013'); grid; xlim([0,470]);
subplot(3,1,2), stairs(lla(:,1),lla(:,12));
xlabel('depth'); ylabel('LOI(%)'); grid; xlim([1200,1465]);
subplot(3,1,3), stairs(gisp(:,1),gisp(:,2),'k');
xlim([0,13000]);grid

% The amount of carbon

figure,
for i = 1:length(names)
    subplot(ceil(length(names)/2),2,i),plot(dat(:,3),dat(:,3+(i))*100,'x-k');grid;xlabel('Depth (mm)');ylabel('% of cps')
%     hold on, scatter(data_xrf(:,3),(data_xrf(:,6+i)-data_xrf(7+i)),10,'--r')
%     hold on, plot(data_xrf(:,3),(data_xrf(:,6+i)+data_xrf(7+i)),10,'--r')
    xlim([min(dat(:,3)),max(dat(:,3))]);
    title(names{i});
end

% cores = {'3','2','1'};
% c_depths = ([1465,1365,1365,1265,1265,1207]);
cores = {'1','2','3'};
c_depths = sort([884.5,1265,1265,1510.1,1529.7,1715] ,'descend');

% --- Looping detrend and loading operations over all images ---
%% 

s = 1;

for i= 1:length(cores)
    I{i} = imread('LLAF_'+ string(cores{i})+'.png');
    
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

    richtingcoff = mean(mean(differ{i},1))/mean(tops-bots);
    blocks{i} = (tops-bots)*richtingcoff;

    for j = 1:length(blocks{i})
        detvec{j} = vec(bots(j):tops(j))'-(-(blocks{i}(j)/2):richtingcoff:(blocks{i}(j)/2));
%         detvec{j} = vec(bots(j):tops(j))'
        xvec{j} = (bots(j):1:tops(j));
    end
    
    
    figure,
    for p = 1:length(detvec)
        subplot(4,1,1),plot(xvec{p},detvec{p},':r')
        hold on, plot(xvec{p},movmean(detvec{p},60),'k');grid;set(gca, 'Xdir', 'reverse');set(gca, 'Ydir', 'reverse');
        
        allvec{p,i} = detvec{p};

        subplot(4,1,2),hold on, findpeaks(movmean(detvec{p},30),xvec{p},'MinPeakDistance',100,'MinPeakProminence',15,'Annotate','extents');grid
        [pks{i},locs{i}] =  findpeaks(movmean(detvec{p},30),xvec{p},'MinPeakDistance',100,'MinPeakProminence',15);
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

figure,
for p = 1:(3*17)
    subplot(6,1,[1:2]),plot(alldvec{p},allvec{p},':r')
    hold on, plot(alldvec{p},movmean(allvec{p},60),'k')
    hold on, plot(alldvec{p},movmean(allvec{p},300),'b')
    hold on,
    
    text(alldpkslocs{p},allpks{p}+(160-allpks{p}),num2str((s:(s-1)+(numel(allpks{p})))'),'FontSize',8,'HorizontalAlignment','center')
    s = s+(numel(allpks{p}));
    
end

set(gca, 'Ydir', 'reverse'); grid; xlim([min(c_depths) max(c_depths)]);xlim([850,1720]);  

subplot(6,1,3), stairs(lla(:,1),lla(:,12));
xlabel('depth'); ylabel('LOI(%)'); grid; xlim([850,1720]);

subplot(6,1,4),plot(dat(:,27),dat(:,3+(5))*100,'x-b');grid;xlabel('Depth (mm)');ylabel('% of cps');
title(names{5}); xlim([850,1720]);
hold on, yyaxis right,plot(dat(:,27),dat(:,3+(6))*100,'.-k');grid;xlabel('Depth (mm)');ylabel('% of cps');
title(names{5}); xlim([850,1720]);
legend('Ca','Ti')

subplot(6,1,5),plot(dat(:,27),dat(:,3+(7))*100,'^-b');grid;xlabel('Depth (mm)');ylabel('% of cps')
title(names{5}); xlim([850,1720]);
hold on, yyaxis right,plot(dat(:,27),dat(:,3+(13))*100,'o-r');grid;xlabel('Depth (mm)');ylabel('% of cps')
title(names{5}); xlim([850,1720]);
legend('Mn','Zr')

subplot(6,1,6), stairs(gisp(:,1),gisp(:,2)+2,'b');
hold on, stairs(icecores(:,1),icecores(:,10),'k')
% hold on, stairs(icecores(:,1),icecores(:,3),'g')
% stairs(icecores(:,1),icecores(:,5),'b')
% stairs(icecores(:,1),icecores(:,7),'r')
legend('GISPII','Ice core average','Ngrip','grip','dye3')
xlim([5000,13100]);grid


