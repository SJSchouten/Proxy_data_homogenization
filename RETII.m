clear
close all

Bathymetry = 0
XRFdata = 1
Stratdata = 1
Picturedata = 1

%% Bathymetry

for bathymetry = 1 
    if bathymetry == 1
    bath = xlsread('Retournemer_bathymetry.xlsx','Retournemer_bathymetry');
    figure,scatter3(bath(:,1),bath(:,2),bath(:,3))

    x = bath(:,1);
    y = bath(:,2);
    z = bath(:,3);

    [xq,yq] = meshgrid(((min(x)):0.0001:(max(x))), (min(y):0.0001:(max(y))));

    vq = griddata(x,y,z,xq,yq,'natural');

    figure,
    surf(xq,yq,vq); hold on, plot3(x,y,z,'x');
    end
end

%% General data loading

for XRFdata = 1
   if XRFdata == 1
       [~,~, dat] = xlsread('Excel Retournemer.xlsx','XRF select');
       dat = dat((8:223),:);
       names = {'S','Cl','K','Ca','Ti','Cr','Mn','Fe','Ni','Cu','Zn','As','Rb','Sr','Zr','Mo','Sn','Ba','Pb'};
       dat = sortrows(dat,1);
       xrf.datorg = dat;
       halfsize = floor((size(dat))/2);
       dats = zeros([216,20]);
       for i = 1:halfsize(2)-2; 
           xrf.ranges(i) = range(cell2mat(dat(:,(4+(2*i)))));
           dats(:,1) = cell2mat(dat(:,1));
           cores(:,1) = (dat(:,3));
           dats(:,(i+1)) = cell2mat(dat(:,(4+(2*i))));
       end
       
       k = 1;z = 1;b = 1;l = 1;g = 1;
       for m = 1:length(dats); 
           if dat{m,5} == 'pla'; xrf.plastic = dats(m,:); 
           elseif dat{m,5} == 'dou'; xrf.double(b,:) = dats(m,:);b = b+1; 
           elseif dat{m,5} == 'lyr'; xrf.layer(l,:) = dats(m,:);l = l+1;
           elseif dat{m,5} == 'res'; xrf.hires(z,:) = dats(m,:);z = z+1;
           elseif dat{m,5} == 'com'; xrf.compare(k,:) = dats(m,:); xrf.compcorename(k,:) = cores(m,:); k = k+1;
           else xrf.regular(g,:) = dats(m,:); g = g+1; end;
       end
       
       xrf.signpla = (xrf.plastic(2:20)./xrf.ranges)*100;
       relevant.double = xrf.double(:,[6,8,9,13,16]);
       relevant.layer = xrf.layer(:,[6,8,9,13,16]);
       relevant.hires = xrf.hires(:,[6,8,9,13,16]);
       relevant.compare = xrf.compare(:,[6,8,9,13,16]);
       relevant.regular = xrf.regular(:,[6,8,9,13,16]);
       
       sz = size(xrf.double);
       for a = 1:sz(1)/2; 
           xrf.dblper(a, 1) = xrf.double((a*2)-1,1);
           xrf.dblper(a,(2:sz(2))) = abs(((xrf.double((a*2)-1,(2:20))-xrf.double(a*2,(2:20)))./xrf.ranges)*100);
       end
       
       indice = [find(nanmean(xrf.dblper(:,(2:20))) > 10)+1,find(xrf.signpla > 10)+1,3];
       xrf.compare(:,indice) = [];
       xrf.hires(:,indice) = [];
       xrf.layer(:,indice) = [];
       xrf.regular(:,indice) = [];
       dats(:,indice) = [];
       names(:,indice-1) = [];
         
       k = 1;z = 1;b = 1
       for m = 1:length(xrf.compcorename)
           if mean(xrf.compcorename{m} == 'Qref') == 1; comp.hov(b,:) = xrf.compare(m,:); b = b+1;
           elseif mean(xrf.compcorename{m} == 'Qrnp') == 1; comp.hovnoplastic(z,:) = xrf.compare(m,:); z = z+1;
           elseif mean(xrf.compcorename{m} == 'Qsam') == 1; comp.closecont(k,:) = xrf.compare(m,:); k = k+1; end
       end
       
       figure,
       for v = 1:length(comp.closecont)-1
       subplot(ceil(length(comp.closecont)/2),2,v),plot(comp.closecont(:,1), comp.closecont(:,(v+1)),'-or');
       hold on, plot(comp.hov(:,1), comp.hov(:,(v+1)),'-+b');
       hold on,plot(comp.hovnoplastic(:,1), comp.hovnoplastic(:,(v+1)),'-xb');title(names{v});legend('Close contact with foil','Hover with foil','Hover without foil');xlim([2750,2830]);
       grid;xlabel('Depth(cm)');ylabel('ppm');
       end
       
       tr = size(xrf.hires);
       tr_1 = size(xrf.layer);
       figure,
       for v = 1:(tr(2)-1)
       subplot(ceil((tr(2)-1)/2),2,v),scatter(xrf.hires(:,1),xrf.hires(:,(v+1)),50,xrf.hires(:,(v+1)),'filled');title(names{v});xlabel('Depth (cm)');ylabel('ppm')
       hold on, scatter(xrf.layer(:,1),xrf.layer(:,(v+1)),50,xrf.layer(:,(v+1)),'x');grid
       hold on, plot(xrf.compare(:,1),xrf.compare(:,(v+1)),'--k')
       hold on, scatter(xrf.regular(:,1),xrf.regular(:,(v+1)),30,xrf.regular(:,(v+1)),'+');
       hold on, plot(xrf.regular(:,1),xrf.regular(:,(v+1)),'-k');
       end
       legend('Highresolution layer scans','Layer spot scans','Comparisson scans core Q','Other systematical data points','Other data line')
       
       tr_2 = size(relevant.layer)
       select = [3,5,6,8,11]
       figure,
       for v = 1:tr_2(2)
       subplot(tr_2(2),1,v),scatter(xrf.hires(:,1),relevant.hires(:,v),50,relevant.hires(:,v),'filled');title(names{select(v)});xlabel('Depth (cm)');ylabel('ppm')
       hold on, scatter(xrf.layer(:,1),relevant.layer(:,v),50,relevant.layer(:,v),'square','filled');grid
       hold on, plot(xrf.compare(:,1),relevant.compare(:,v),'--k')
       hold on, scatter(xrf.regular(:,1),relevant.regular(:,v),30,relevant.regular(:,v),'+');
       hold on, plot(xrf.regular(:,1),relevant.regular(:,v),'-k');
       end
       legend('Highresolution layer scans','Layer spot scans','Comparisson scans core Q','Other systematical data points','Other data line')
   end
end

ret = xlsread('Excel Retournemer 2.xlsx','Relevant data');
gisp = xlsread('GISPII_oxygen_isotopes.xls');
datas = xlsread('Excel Timmelteich.xlsx','Timmelteich_2020_LOI_composite');
icecores = xlsread('groenland comb');
icecores = sortrows(icecores,1);

figure, subplot(3,1,1),stairs(ret(:,2),ret(:,10)); xlim([2350,2930]); xlabel('Depth (cm)'); grid; ylabel('LOI(%)');title('Loss on ignition percentage');
hold on, plot(ret(:,2),movmean(ret(:,10),70,'omitnan'))
subplot(3,1,3), stairs(ret(:,2),ret(:,6));xlabel('Depth (cm)'); ylabel('\rho_{organic}'); grid; xlim([2350,2930]);
hold on ,yyaxis right, stairs(ret(:,2),ret(:,7));xlabel('Depth (cm)'); ylabel('\rho_{clastic}'); grid; xlim([2350,2930]); grid; title('Individual density components');legend('\rho_{organic}','\rho_{clastic}')
% subplot(3,1,2), stairs(ret(:,2),(ret(:,7)./ret(:,6))); xlabel('Depth(cm)'); ylabel('\rho_{clastic}/\rho_{organic}'); grid; xlim([2350,2930]); set(gca ,'Ydir','reverse');title('\rho_{clastic} / \rho_{organic} ratio')
movmeananom = ret(:,10)-movmean(ret(:,10),70,'omitnan')
bits = movmeananom > 0
low = movmeananom(bits)
high = movmeananom(bits == 0)
subplot(3,1,2), bar(ret(bits,2),low,'r');
hold on, bar(ret(bits == 0,2),high,'b');
xlabel('Depth (cm)'); grid; ylabel('LOI anomaly (%)')

figure,subplot(4,1,1), stairs(ret(:,4),ret(:,6));
xlabel('depth'); ylabel('CD'); grid; xlim([7000,14000]);
subplot(4,1,2), stairs(ret(:,4),ret(:,7));
xlabel('depth'); ylabel('KlasD'); grid; xlim([7000,14000]);
subplot(4,1,3), stairs(ret(:,4),ret(:,10));
xlabel('depth'); ylabel('LOI(%)'); grid; xlim([7000,14000]);
subplot(4,1,4), stairs(icecores(:,1),icecores(:,10),'k')
hold on, stairs(icecores(:,1),icecores(:,3),'g')
stairs(icecores(:,1),icecores(:,5),'b')
stairs(icecores(:,1),icecores(:,7),'r')
legend('Ice core average','Ngrip','grip','dye3')
xlim([7000,14000]);grid

norm7d = diff(ret(:,7))/(max(diff(ret(:,7)))-min(diff(ret(:,7))));
norm6d = diff(ret(:,6))/(max(diff(ret(:,6)))-min(diff(ret(:,6))));
norm7 = ret(:,7)/(max(ret(:,7))-min(ret(:,7)));
norm6 = ret(:,6)/(max(ret(:,6))-min(ret(:,6)));

figure,subplot(4,1,1),plot(ret(:,4),(ret(:,7)/(max(ret(:,7)))-ret(:,6)/(max(ret(:,6)))))
subplot(4,1,2),plot(ret(:,4),norm7)
hold on, yyaxis right, plot(ret(:,4),norm6)
subplot(4,1,3), plot(ret((1:506),4),norm7d)
hold on,yyaxis right, plot(ret((1:506),4),norm6d)
subplot(4,1,4), plot(ret((1:506),4),norm6d-norm7d)

% 
% figure,
% for i = 1:19
%     subplot(19,1,i),plot(dat(:,1),dat(:,3+(i*2)),'x-k');grid;xlabel('depth');ylabel('ppm')
% %     hold on, scatter(data_xrf(:,3),(data_xrf(:,6+i)-data_xrf(7+i)),10,'--r')
% %     hold on, plot(data_xrf(:,3),(data_xrf(:,6+i)+data_xrf(7+i)),10,'--r')
%     xlim([2390,2900]);
%     title(names{i});
% end

%% --- Running the Age depth model for retournemer core ---

inputfilename = 'Excel Retournemer 2.xlsx'
inputfiletabname = 'Retournemer_2020_LOI_composite'
exponent = [0.6:0.05:0.9];

depthfix = [2905.5,2859.5,2854.5,2596.5]; % depthfix = [2905.5,2895.5,2859.5,2854.5,2801.5,2596.5];
timefix = [15842,12950,12900,9250]; % timefix = [15842,14542,12950,12900,11650,9250];
fixnames = {'Pleniglacial end','LST','LST','Tilia'};

tiepnts_polx = [2762.5,2762.5,2595,2595];
tiepnts_poly = [10800,11000,9400,9200];
polnames = {'Corylus','Tilia'};
endpoint = 2392.5;
choice = 0.8;
method = 1;

[out, fig, f] = Agedepthmodelling_ret(inputfilename,inputfiletabname,exponent,timefix,depthfix,fixnames,tiepnts_polx,tiepnts_poly,polnames,choice,method,endpoint);

% --- Defining the cores ---
%%

cores = {'R','Q','P','O','N','M','L','K','J'};
c_depths = ([2940,2840,2840,2740,2740,2725,2690,2590,2590,2490,2490,2392.5,2390,2290,2290,2190,2190,2106]);

% --- Looping detrend and loading operations over all images ---
%% 

s = 1;

for i= 1:length(cores)
%% --- Preparing the depth scale and scanline positions --- 
    
    I{i} = imread('scan_'+ string(cores{i})+'_final.png');
    
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

%% --- Making scanlines ---

    sommetje = zeros(sz(2),1);
    
    for z = 1:length(lines)
        y = [lines(z),lines(z)];
        x = [1,sz(2)];
        n = improfile(Im,x,y);
        vect = n(:,:,1);
        sommetje = sommetje + vect;
    end
    
    vec = flip(sommetje/length(lines));
    xvect = [1:length(vec)]';
    count = 0;
    
%% --- Decracking the core ---
    
    % Manual decracking
    
        delpositie{3} = [{5104:1:5232}];
        delpositie{4} = [{9885:1:10050},{17630:1:17720},{28770:1:28890}];
        delpositie{5} = [{6518:1:6587},{7249:1:7313},{11810:1:11920},{12800:1:12880},{13960:1:14110},{17000:1:17140},{18920:1:19040},{19990:1:20140},{21030:1:21070},{23300:1:23500},{23730:1:23800},{24230:1:24390},{25640:1:25720}];
        delpositie{6} = [{6144:1:6273},{8782:1:8923},{9845:1:9947},{15120:1:15210},{21190:1:21510},{29600:1:29640}];
        delpositie{7} = [{8767:1:8954},{13810:1:14060},{14560:1:14660},{17020:1:17200},{18640:1:18880},{20470:1:20830},{21720:1:22100},{23380:1:23530},{24300:1:24600},{25760:1:26130},{27210:1:27620},{29060:1:29530}];
        delpositie{8} = [{2601:1:2684},{4869:1:5246},{5594:1:5665},{5929:1:6056},{11050:1:11220},{11280:1:11400},{12450:1:12640},{14570:1:14670},{16790:1:16980},{17880:1:18110},{18440:1:18550},{19690:1:19840},{20250:1:20350},{20510:1:20680},{21180:1:21320},{22520:1:22770},{24280:1:24660},{27000:1:27300},{29030:1:29170},{29200:1:29370}]; 
        delpositie{9} = [{13880:1:14170},{14950:1:15060},{17960:1:18300},{18690:1:18900},{21210:1:21660},{21990:1:22050},{23058:1:23350},{24840:1:24970},{25870:1:26470},{27260:1:27730}];
        
        if i == 3 
            vec(delpositie{i}{1}) = vec(delpositie{i}{1})+25; 
        else
            for a = 1:length(delpositie{i})
                delpos = delpositie{i}{a}-count;
                vec(delpos) = [];
                xvect = [1:length(vec)];
                count = count+length(delpos);
            end
        end
        
    % Automatic decracking
    
%         if i > 3
%             figure,
%             [cracksdepth{i},crackloc{i},crackwidth{i}] = findpeaks(movmean((-vec),30),xvect,'MinPeakHeight',-80,'MinPeakDistance',100,'MinPeakProminence',30);
%             [a,b,crackwidthref{i}] = findpeaks(movmean((-vec),10),xvect,'MinPeakDistance',100,'MinPeakProminence',30);
%             findpeaks(movmean((-vec),30),xvect,'MinPeakHeight',-80,'MinPeakDistance',100,'MinPeakProminence',30,'Annotate','extents')
%             
%             for a = 1:length(crackwidth{i})
%                 if crackwidth{i}(a) > 300 
%                     crackwidth{i}(a) = 200;
%                 end
%                 delposit{a} = [crackloc{i}(a)-count - floor((crackwidth{i}(a))/2), floor((crackwidth{i}(a))/2) + crackloc{i}(a)-count];
%                 delpos{i,a} = [(find(xvect == delposit{a}(1))):1:(find(xvect == delposit{a}(2)))];
%                 vec(delpos{i,a}) = [];
%                 xvect = [1:length(vec)]';
%                 count = count+length(delpos{i,a})
%             end
%             hold on,plot(movmean((-vec-50),30))
%         end
        
%% --- Full Detrending alligning and reorganizing data ---
        
    stretch = length(vec);
    ims{i} = [d_scale(2):(d_scale(1)-d_scale(2))/stretch:d_scale(1)];
    
   
    inc = find(vec > 252); 
    bots = inc(find(diff(inc) > 1));
    tops = inc(find(diff(inc) > 1)+1);

    % General detrending of holocene cores
    
    if i >= 3
        vec = detrend(vec) + 100;
    elseif i == 1
        vec = vec+5;
    end

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
        xvec{j} = (bots(j):1:tops(j));  
    end
    
%% --- Figure plotting and peak finding ---
    
    figure,
    for p = 1:length(detvec)
        subplot(5,1,1),plot(xvec{p},detvec{p},':r')
        hold on, plot(xvec{p},movmean(detvec{p},60),'k');grid
        hold on,
        
        grid;set(gca, 'Ydir', 'reverse');title('core '+string(cores{i})+' Greyscale scan'); set(gca, 'Xdir', 'reverse')

        if i == 1
            allxvec{p,i} = xvec{p};
        else
            allxvec{p,i} = xvec{p} + max(allxvec{length(blocks{i-1}),i-1});
        end
        
        % Creating event based stratigraphy
        
        if i == 5
            subplot(5,1,2),findpeaks(movmean(detvec{p},30),xvec{p},'MinPeakDistance',100,'MinPeakProminence',20,'Annotate','extents')
            [pks{i},locs{i}] =  findpeaks(movmean(detvec{p},30),xvec{p},'MinPeakDistance',100,'MinPeakProminence',20);
            text(locs{i}',pks{i}+10',num2str((s:(s-1)+(numel(pks{i})))'),'FontSize',8,'HorizontalAlignment','center')
        elseif i < 3 
            subplot(5,1,2),findpeaks(movmean(detvec{p},30),xvec{p},'MinPeakDistance',80,'MinPeakProminence',12,'Annotate','extents')
            [pks{i},locs{i}] =  findpeaks(movmean(detvec{p},30),xvec{p},'MinPeakDistance',80,'MinPeakProminence',12);
            text(locs{i}',pks{i}+10',num2str((s:(s-1)+(numel(pks{i})))'),'FontSize',8,'HorizontalAlignment','center')
        else    
            subplot(5,1,2),findpeaks(movmean(detvec{p},30),xvec{p},'MinPeakDistance',100,'MinPeakProminence',15,'Annotate','extents')
            [pks{i},locs{i}] =  findpeaks(movmean(detvec{p},30),xvec{p},'MinPeakDistance',100,'MinPeakProminence',15);
            text(locs{i}',pks{i}+10',num2str((s:(s-1)+(numel(pks{i})))'),'FontSize',8,'HorizontalAlignment','center')
        end
        hold on,
        s = s+(numel(pks{i}));
        
        allvec{p,i} = detvec{p};
        allpks{p,i} = pks{i};
        allpkslocs{p,i} = locs{i};
        
    end
    
    grid;set(gca, 'Ydir', 'reverse');title('core '+string(cores{i})+' Greyscale scan'); set(gca, 'Xdir', 'reverse')
    subplot(5,1,3), stairs(ret(:,2),ret(:,10));
    xlabel('depth')
    ylabel('LOI(%)');xlim([d_scale(2) d_scale(1)])
    subplot(5,1,4),plot(ims{i}(1:length(vec)),vec);grid;set(gca, 'Ydir', 'reverse');set(gca, 'Xdir', 'reverse')
    subplot(5,1,5),imshow(I{i}); title('core '+string(cores{i}));

%% --- Generating block based depthvectors ---
    
    o(1) = 1;
    for y = 2:length(tops)+1
        o(y) = length(bots(y-1):tops(y-1));
    end
    sum(o);
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
    
    clear xvec;
    clear detvec;
    clear o;
    clear delposit;
    clear vec;
    clear xvect;
    clear count;
end

s = 1;

figure,
for p = 1:(19*9)
    subplot(10,1,[1:3]),plot(alldvec{p},allvec{p},':r')
    hold on, plot(alldvec{p},movmean(allvec{p},60),'k')
    hold on, plot(alldvec{p},movmean(allvec{p},300),'b')
    hold on,
    text(alldpkslocs{p},allpks{p}+(160-allpks{p}),num2str((s:(s-1)+(numel(allpks{p})))'),'FontSize',8,'HorizontalAlignment','center')
    s = s+(numel(allpks{p}));
end
ylim([50, 190]); set(gca, 'Ydir', 'reverse'); grid; xlim([min(c_depths) max(c_depths)]); 
subplot(10,1,4), stairs(ret(:,2),ret(:,10));
xlabel('depth')
ylabel('LOI(%)');xlim([min(c_depths) max(c_depths)]);grid

% for p = [1:35 39:41 58:73 77:88 96:114]
%     ser = movmean(allvec{p},50) - movmean(allvec{p},500);
%     sigm = (max(ser)-min(ser))*0.10;
%     sigma = [(max(ser)-sigm),(min(ser)+sigm)];
%     binma = ser > sigma(1);
%     binmi = ser < sigma(2);
% 
%     subplot(4,1,3), yyaxis left, plot(alldvec{p},movmean(allvec{p},50),'k-')
%     hold on, plot(alldvec{p},movmean(allvec{p},500),'b-')
%     ylim([60, 185]);set(gca, 'Ydir', 'reverse'); grid; xlim([min(c_depths) max(c_depths)])   
% %     yyaxis right, bar(alldvec{p},binmi,'g')
%     yyaxis right, bar(alldvec{p},binma,'r')
%     
%     hold on,
% end

tr_2 = size(relevant.layer);
select = [3,5,6,8,11];
for v = 1:tr_2(2)
   subplot(tr_2(2)+5,1,v+4),scatter(xrf.hires(:,1),relevant.hires(:,v),50,relevant.hires(:,v),'filled');title(names{select(v)});ylabel('ppm');
   hold on, scatter(xrf.layer(:,1),relevant.layer(:,v),50,relevant.layer(:,v),'square','filled');grid; xlim([min(c_depths) max(c_depths)]);
   hold on, plot(xrf.compare(:,1),relevant.compare(:,v),'--k');
   hold on, scatter(xrf.regular(:,1),relevant.regular(:,v),30,relevant.regular(:,v),'+');
   hold on, plot(xrf.regular(:,1),relevant.regular(:,v),'-k');
end
legend('Highresolution layer scans','Layer spot scans','Comparisson scans core Q','Other systematical data points','Other data line')

subplot(10,1,10), stairs(gisp(:,1),gisp(:,2)+2,'b');
hold on, stairs(icecores(:,1),icecores(:,10),'k');
% hold on, stairs(icecores(:,1),icecores(:,3),'g')
% stairs(icecores(:,1),icecores(:,5),'b')
% stairs(icecores(:,1),icecores(:,7),'r')
legend('GISPII','Ice core average','Ngrip','grip','dye3')
xlim([8000,13000]);grid;

% detvec{j}
% figure,plot(tops(1:18),diff)
% for k = 1:length(vec)
% detr(k,:) = vec(k,:)-(difference/length(vec)*k);
% end

% figure,
% for p = 1:(19*6)
%     subplot(2,1,1),plot(alldvec{p},allvec{p},':r')
%     hold on, plot(alldvec{p},movmean(allvec{p},60),'k')
%     hold on, plot(alldvec{p},movmean(allvec{p},300),'b')
%     hold on,
%     text(alldpkslocs{p},allpks{p}+(160-allpks{p}),num2str((s:(s-1)+(numel(allpks{p})))'),'FontSize',8,'HorizontalAlignment','center')
%     s = s+(numel(allpks{p}));
% end
% ylim([50, 190]); set(gca, 'Ydir', 'reverse'); grid; xlim([min(c_depths) max(c_depths)]); 
% subplot(2,1,2), stairs(ret(:,2),ret(:,10));
% xlabel('depth')
% ylabel('LOI(%)');xlim([min(c_depths) max(c_depths)]);grid;

%% Wholescale age devision

%     counter = 0
    
%     for h = 1:(sz(1)*sz(2))
%         counter = counter + length(alldvec{h});
%         boundaries = [min(alldvec{h}),max(alldvec{h})]
%     end
%     counter/(c_depths(1)-c_depths(2));
%     temp = [1,find(out.Depth == (c_depths(length(c_depths))))]
%     ageselect = out.Age(temp(1):temp(2))
%     find(alldvec{1} == out.Depth(1)) 
%     depthselect = out.Depth(temp(1):temp(2))
%     d_deler = out.Depth(1) - (c_depths(length(c_depths)))
%     numerix = find((alldvec{7} > 2905.503 | alldvec{7} < 2905.497) == 0)
%     
%     counter2 = 0
%     for h = 1:6
%         counter2 = counter2 + length(alldvec{h})
%     end
%     compfact = counter2 + numerix
% 
%     [age_res,depth_res] = resample(ageselect,depthselect,(counter-compfact)/d_deler,1,1)
%     allagevec = flip(age_res(1:length(age_res)-1))
%     
%     figure, plot(age_res,depth_res,'-x')
%     counter-compfact
     
%% Core specific scaled age determination

sz = size(allvec);

% Creating wholecorescale greyscale vectors

for h = 1:sz(2)
    tel = 1;
    te = 1;
    for w = 1:sz(1)
        c.longvector{h}(tel:(tel-1+length(allvec{w,h}))) = allvec{w,h};
        c.longdepthvector{h}(tel:(tel-1+length(alldvec{w,h}))) = alldvec{w,h};
        c.pksdlocs{h}(te:(te-1+length(alldpkslocs{w,h}))) = alldpkslocs{w,h};
        c.allpeaks{h}(te:(te-1+length(allpks{w,h}))) = allpks{w,h};
        tel = tel+length(allvec{w,h});
        te = te+length(alldpkslocs{w,h});
    end
end

for h = 1:9
    for t = 1:length(c.pksdlocs{h}); c.pksdlocs_ind{h}(t) = find(c.longdepthvector{h} == c.pksdlocs{h}(t)); end;
    if h == 1; tempnew = c.pksdlocs_ind{1} - 10000; c.pksdlocs_ind{1} = tempnew(tempnew > 0); c.allpeaks{1} = c.allpeaks{1}(tempnew > 0); end
end 

% Start position correction
c.longdepthvector{1} = c.longdepthvector{1}(10000:end);
c.longvector{1} = c.longvector{1}(10000:end);

% New core depths with first core start substituted
new_c_depths = c_depths;
new_c_depths(1) = 2905.5;

% Resampling age and depth vectors comming from the Age depth mdoel
[re.age,re.de] = resample(out.Age,out.Depth,2,1,1);
gevonden = find(re.de == new_c_depths);
for j = 1:length(gevonden); gevondenecht(j) = gevonden(j) - (j-1)*length(re.de);end
plusje = 2;
for i = 1:length(gevondenecht)/2; c.a{i} = re.age(gevondenecht(14-(plusje)):gevondenecht(14-(plusje+1))); c.d{i} = re.de(gevondenecht(14-plusje):gevondenecht(14-(plusje+1))); plusje  = plusje + 2; end

% positions
    
%     c.a{1} = out.Age(1:68)
%     c.d{1}  = out.Depth(1:68)
%     c.a{2} = out.Age(69:168)
%     c.d{2}  = out.Depth(69:168)
%     c.a{3} = out.Age(169:183)
%     c.d{3}  = out.Depth(169:183)
%     c.a{4} = out.Age(190:289)
%     c.d{4}  = out.Depth(190:289)
%     c.a{5} = out.Age(290:389)
%     c.d{5} = out.Depth(290:389)
%     c.d{6}  = out.Depth(390:488)
%     c.a{6} = out.Age(390:488)
    
% Resampling and plotting the new resampled age depth model

figure,
for h = 1:9; 
    c.longdepth{h} = [max(c.longdepthvector{h}),min(c.longdepthvector{h})]; 
    c.lengths{h} = length(c.longdepthvector{h});
    c.increment{h} = c.longdepthvector{h}(1)-c.longdepthvector{h}(2);
    [c.ageres{h},c.depthres{h}] = resample(c.a{7-h},c.d{7-h},1/c.increment{h},1,1,'spline');
    subplot(4,1,1),plot(c.ageres{h},c.depthres{h},'.-','LineWidth',2); hold on,
end
xlabel('Depth(cm)');ylabel('Age (cal. yr BP)');grid;title('Age depth model segemented and resampled');xlim([7000,16000]);
hold on, scatter(timefix,depthfix,40,'o','filled'); hold on, text(timefix+40,depthfix,fixnames)

c.shift{1} = c.depthres{1}(2:end-1) - sort(c.longdepthvector{1}','ascend');
c.shift{2} = c.depthres{2}(2:end) - sort(c.longdepthvector{2}','ascend');
c.shift{3} = c.depthres{3}(2:end) - sort(c.longdepthvector{3}','ascend');
c.shift{4} = c.depthres{4}(2:end-1) - sort(c.longdepthvector{4}','ascend');
c.shift{5} = c.depthres{5}(2:end) - sort(c.longdepthvector{5}','ascend');
c.shift{6} = c.depthres{6}(2:end-1) - sort(c.longdepthvector{6}','ascend');

c.ageres{1} = flip(c.ageres{1}(2:end-1));
c.ageres{2} = flip(c.ageres{2}(2:end));
c.ageres{3} = flip(c.ageres{3}(2:end));
c.ageres{4} = flip(c.ageres{4}(2:end-1)); 
c.ageres{5} = flip(c.ageres{5}(2:end));
c.ageres{6} = flip(c.ageres{6}(2:end-1));

s = 1;

subplot(4,1,2),
for n = 1:length(cores)
    c.pksages{n} = c.ageres{n}(c.pksdlocs_ind{n});
    plot(c.ageres{n},c.longvector{n},'k'); hold on,
    plot(c.ageres{n},movmean(c.longvector{n},60),'r','LineWidth',2), hold on,
    text(c.pksages{n},c.allpeaks{n}+(160-c.allpeaks{n}),num2str((s:(s-1)+(numel(c.allpeaks{n})))'+5),'FontSize',8,'HorizontalAlignment','center')
    s = s+(numel(c.allpeaks{n})); hold on,
end

text(timefix+40,[50,50,50,50],fixnames); hold on,
scatter(timefix,[50,50,50,50],40,'r','o','filled');
set(gca,'Ydir','reverse');xlim([7000,16000]);ylim([40,180]);grid

subplot(4,1,3), stairs(out.Age,out.LOI);
xlabel('Age (cal. yr BP)')
ylabel('LOI(%)');xlim([7000,16000]);grid

subplot(4,1,4), stairs(gisp(:,1),gisp(:,2)+2,'b');
hold on, stairs(icecores(:,1)-50,icecores(:,10),'k');
% hold on, stairs(icecores(:,1),icecores(:,3),'g')
% stairs(icecores(:,1),icecores(:,5),'b')
% stairs(icecores(:,1),icecores(:,7),'r')
legend('GISPII','Ice core average','Ngrip','grip','dye3');
xlim([7000,16000]);grid;
    
%    (c.lengths{h}/length(c.a{7-h}))*2
%     sz = size(allvec)
%     increment = 0
%     for g = 1:sz(2) 
%         counter2 = 0 
%     for h = 1:(sz(1))
%         counter2 = counter2 + length(alldvec{h,g});
%         boundaries{h,g} = [min(alldvec{h,g}),max(alldvec{h,g})];
%     end
%     c{g} = counter2
%     depths{g} = [c_depths(g+increment),c_depths(g+increment+1)]
%     d{g} = c_depths(g+increment)-c_depths(g+increment+1)
%     increment = increment+1
%     [minval{g},mindex{g}] = min(abs(depth_res-depths{g}(1)))
%     [maxval{g},maxdex{g}] = min(abs(depth_res-depths{g}(2)))
%     counts{g} = mindex{g}-maxdex{g}
%     [agecoreres{g},depthcoreres{g}] = resample(age_res,depth_res,(counts{g}/c{g}),1,1)
%     end
%     d{3} = d{3} + 35


%     figure, plot(allagevec,longdepthvector(compfact:length(longvector)-1))
%     length(longvector)-compfact
    

%     figure, plot(age_res,depth_res,'x-',ageselect,depthselect,'o')
%     hold on, yyaxis right, plot(flip(age_res(1:length(age_res)-1)),longvector)
%     
%     figure, plot(longdepthvector',flip(age_res((1:length(depth_res)-1))))

% figure,subplot(2,1,1),plot((detr-movmean(detr,100)))
% array = (detr-movmean(detr,100))
% normarray = array./(max(array)-min(array))/100
% subplot(2,1,2), plot(normarray)
% 
% ma = 0
% mi = 1
% stretch = length(average)-mi-ma
% 
% i = [2724:216/stretch:2940];
% length(i)
% figure,subplot(3,1,1),plot(dat(:,1),dat(:,13),'x-k');grid;xlabel('depth');ylabel('ppm');xlim([2724,2940]);set(gca, 'Ydir', 'reverse')
% subplot(3,1,2),plot(i,average(mi:mi+stretch))
% hold on, plot(i,movmean(average(mi:mi+stretch),60));grid
% set(gca, 'Ydir', 'reverse');;xlim([2724,2940])
% subplot(3,1,3), stairs(ret(:,3),ret(:,14));
% xlabel('depth')
% ylabel('LOI(%)')
% grid; xlim([2724,2940]);
% 
% figure, plot(detr)
