function p = bootStrapP(met1, met2, nBoots)


    nPts1 = size(met1,1);
    nPts2 = size(met2,1);
    nTot = nPts1 + nPts2;
    metTotal = cat(1,met1,met2);
    
    testStat = meanDiff(met1,met2);
    bootStat = zeros(nBoots,1);
    for bootN = 1:nBoots
        rs1 = metTotal(randi(nTot,nPts1,1),:);
        rs2 = metTotal(randi(nTot,nPts2,1),:);
        bootStat(bootN) = meanDiff(rs1,rs2);
    end
    
    p = (sum(bootStat >= testStat) + 1)/(nBoots + 1);
    
   
%     figure;
%     hist(bootStat,100); title(num2str(p)); hold on;
%     plot([1 1]*testStat,ylim(),'r');
%     pause()
%     close gcf;

    
    
    
function out = meanDiff(met1, met2)

    out = sum(abs(nanmean(met2,1) - nanmean(met1,1)));


    
               
