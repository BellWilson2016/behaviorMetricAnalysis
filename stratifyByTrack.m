function [meanMetrics, meanIX] = stratifyByTrack(perGateMeas, measIX, tAIX)

   for trackIX = 1:size(tAIX,1)
		ix = find(measIX == trackIX);
		meanMetrics(trackIX) = nanmean(perGateMeas(ix));
	end

	meanIX = tAIX;
    
