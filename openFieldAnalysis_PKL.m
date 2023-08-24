[fileName, folderName] = uigetfile('*.xlsx');
xlsfile = fullfile(folderName,fileName);
[~, sheetNames] = xlsfinfo(xlsfile);


for sheetIdx = 1:length(sheetNames)
    sheetData{sheetIdx}= xlsread(xlsfile,sheetNames{sheetIdx});
end


[~,txtString,~] = xlsread(xlsfile,sheetNames{sheetIdx});
txtString = txtString(1,:);


for sheetIdx = 1:length(sheetNames)
    temp = sheetData{sheetIdx};
    
    for colIdx = 1
        timeStamps = temp(:,colIdx);
        RangeIdx = timeStamps<300;
        timeStamps = timeStamps(RangeIdx);
        xCoor =  temp(RangeIdx,colIdx+1);
        yCoor =  temp(RangeIdx,colIdx+2);
        diff(timeStamps);
        velo = sqrt(xCoor(1:end-1).^2 + yCoor(1:end-1).^2)./diff(timeStamps);
        txtString{colIdx};
        
        veloIdx = velo>16000;
        sum(veloIdx)*mean(diff(timeStamps));
        mean(velo);
        mean(velo(veloIdx));
        velo (veloIdx);
        mean(velo(~veloIdx));
        finalData.(sheetNames{sheetIdx}).(txtString{colIdx}).meanVelo=mean(velo);
        finalData.(sheetNames{sheetIdx}).(txtString{colIdx}).activemeanVelo=mean(velo);
        finalData.(sheetNames{sheetIdx}).(txtString{colIdx}).xCoor=xCoor;
        finalData.(sheetNames{sheetIdx}).(txtString{colIdx}).yCoor=yCoor;
%          finalDC{1} =txtString{colIdx}
%          finalDC{2} =mean(velo);
% 
%          figure
%          plot(xCoor,yCoor)
%          title(sprintf('%s%s',sheetNames{sheetIdx},txtString{colIdx})) 
%          saveas(gcf,sprintf('Openfield_%s%s',sheetNames{sheetIdx},txtString{colIdx}),'pdf');
%     
        %Calculate total distance traveled
        DistanceTraveled=zeros(1, 8999); % Aaron
        for i = 2:8999
        xCoorDist = xCoor(i,1) - xCoor(i-1,1);
        yCoorDist = yCoor(i,1) - yCoor(i-1,1);
        DistanceTraveled(i) = sqrt(xCoorDist^2+yCoorDist^2);
        end
        
        finalData.(sheetNames{sheetIdx}).(txtString{colIdx}).TotalDistanceTraveled=sum(DistanceTraveled);
    end
end
%%
DistanceWT = [finalData.PKL20.Timestamp.TotalDistanceTraveled ...
    finalData.PKL23.Timestamp.TotalDistanceTraveled ...
    finalData.PKL24.Timestamp.TotalDistanceTraveled ...
    finalData.PKL26.Timestamp.TotalDistanceTraveled ...
    finalData.PKL33.Timestamp.TotalDistanceTraveled finalData.PKL37.Timestamp.TotalDistanceTraveled finalData.PKL39.Timestamp.TotalDistanceTraveled ...
    finalData.PKL41.Timestamp.TotalDistanceTraveled finalData.PKL42.Timestamp.TotalDistanceTraveled];
DistanceHT = [finalData.PKL1.Timestamp.TotalDistanceTraveled finalData.PKL2.Timestamp.TotalDistanceTraveled finalData.PKL3.Timestamp.TotalDistanceTraveled finalData.PKL5.Timestamp.TotalDistanceTraveled finalData.PKL6.Timestamp.TotalDistanceTraveled finalData.PKL7.Timestamp.TotalDistanceTraveled finalData.PKL8.Timestamp.TotalDistanceTraveled finalData.PKL9.Timestamp.TotalDistanceTraveled finalData.PKL10.Timestamp.TotalDistanceTraveled finalData.PKL11.Timestamp.TotalDistanceTraveled finalData.PKL12.Timestamp.TotalDistanceTraveled finalData.PKL13.Timestamp.TotalDistanceTraveled finalData.PKL14.Timestamp.TotalDistanceTraveled finalData.PKL16.Timestamp.TotalDistanceTraveled finalData.PKL17.Timestamp.TotalDistanceTraveled finalData.PKL18.Timestamp.TotalDistanceTraveled finalData.PKL19.Timestamp.TotalDistanceTraveled finalData.PKL21.Timestamp.TotalDistanceTraveled finalData.PKL25.Timestamp.TotalDistanceTraveled finalData.PKL30.Timestamp.TotalDistanceTraveled finalData.PKL35.Timestamp.TotalDistanceTraveled finalData.PKL38.Timestamp.TotalDistanceTraveled  finalData.PKL43.Timestamp.TotalDistanceTraveled];
DistanceHO = [finalData.PKL4.Timestamp.TotalDistanceTraveled finalData.PKL15.Timestamp.TotalDistanceTraveled finalData.PKL22.Timestamp.TotalDistanceTraveled finalData.PKL27.Timestamp.TotalDistanceTraveled finalData.PKL28.Timestamp.TotalDistanceTraveled finalData.PKL29.Timestamp.TotalDistanceTraveled finalData.PKL31.Timestamp.TotalDistanceTraveled finalData.PKL32.Timestamp.TotalDistanceTraveled finalData.PKL34.Timestamp.TotalDistanceTraveled finalData.PKL36.Timestamp.TotalDistanceTraveled finalData.PKL40.Timestamp.TotalDistanceTraveled finalData.PKL44.Timestamp.TotalDistanceTraveled];

DistanceWT = DistanceWT.*0.11 %cm/pixel conversion
DistanceHT = DistanceHT.*0.11 %cm/pixel conversion
DistanceHO = DistanceHO.*0.11 %cm/pixel conversion

avgDistanceWT = sum(DistanceWT)/length(DistanceWT)
avgDistanceHT = sum(DistanceHT)/length(DistanceHT) 
avgDistanceHO = sum(DistanceHO)/length(DistanceHO)

% avgDistanceWTcm = (sum(DistanceWT)/length(DistanceWT))/0.13 % Need to
% translate pixels into cm... fastest 

[hHT,pHT,ciHT,statsHT] = ttest2(DistanceWT,DistanceHT)
[hHO,pHO,ciHO,statsHO] = ttest2(DistanceWT,DistanceHO)

stderrWT=std(DistanceWT)./sqrt(length(DistanceWT));
stderrHT=std(DistanceHT)./sqrt(length(DistanceHT));
stderrHO=std(DistanceHO)./sqrt(length(DistanceHO));


figure 
Distance =[avgDistanceWT, avgDistanceHT, avgDistanceHO];
stderr =[stderrWT stderrHT stderrHO]
errorbar(Distance, stderr)
yl = get(gca,'ylim');
hold on
bar(Distance)
hold off
title 'Distance Traveled by Genotype'
xlabel Genotype
ylabel ('Distance (cm)')
set(gca,'XTickLabel',{'','Wild Type','', 'Hetero','', 'Homo'})
ylim(yl)
% saveas(gcf,'Distance.pdf')
%% Create distance bar graph with only Wild type and Homozygotes.
% figure 
% Distance =[avgDistanceWT, avgDistanceHO];
% stderr =[stderrWT stderrHO]
% errorbar(Distance, stderr)
% yl = get(gca,'ylim');
% hold on
% bar(Distance)
% hold off
% title 'Distance Traveled by Genotype'
% xlabel Genotype
% ylabel ('Distance (cm)')
% set(gca,'XTickLabel',{'','Wild Type','', 'Homo'})
% ylim(yl)
% saveas(gcf,'DistanceWTHO.pdf')

