% tAIX = [powerN, LR, flyN, trialN]
function [meanMetrics, meanIX] = stratifyByFly(metrics, metricIX, tAIX)

    flyList = unique(tAIX(:,3));
    powerList = unique(tAIX(:,1));
    
    meanMetrics = [];
    meanIX = [];
    for flyNn = 1:length(flyList)
        flyN = flyList(flyNn);
        for powerNn = 1:length(powerList)
            powerN = powerList(powerNn);
        
            ix = find((tAIX(metricIX,3) == flyN) & (tAIX(metricIX,1) == powerN));
            meanVal = nanmean(metrics(ix,:),1);
            
            meanMetrics = cat(1, meanMetrics,meanVal);
            meanIX      = cat(1, meanIX, [powerN flyN]);
            
        end
    end
    
    
