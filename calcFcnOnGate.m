function [perGateMeas, measIX] = calcFcnOnGate(tA, tAIX, gIX, stSamp, enSamp, calcFcn)
    
    measQty = calcFcn( tA, tAIX);
    
    % For each gating event, take the avg. of measQty over the window
    for gateN = 1:length(gIX)
        IX = gIX(gateN);
        
        perGateMeas(gateN,1) = nanmean(measQty(IX, [stSamp(gateN):enSamp(gateN)]));
    end
    
    measIX = gIX;
    
    
    
    
