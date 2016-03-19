function timeAnalysis()

	load('~/RTFW/Code/kennedyAnalysis/PTkennedyAnalysis.mat');

	allGenos = unique(PTgenoList);
	filtWidth = 8;
	filtRange = [(1+filtWidth):(128-filtWidth)];

	for metricN = 1:14
		figure;
		for genoNn = 1:length(allGenos)
			genoN = allGenos(genoNn);

			metVal = zeros(128,1);
			for trialN = (1+filtWidth):(128-filtWidth)
				trialBot = trialN - filtWidth;
				trialTop = trialN + filtWidth;
				ix = find((PTgenoList == genoN) & (PTpowerNList >= 7) &...
						  (PTtrialNList >= trialBot) & (PTtrialNList <= trialTop));
				metVal(trialN) = nanmean(scoresByTrack(ix,metricN));
			end
			
			ffsubplot(2,1,1);
			plot(filtRange, metVal(filtRange),'.-','Color',pretty(genoNn)); hold on;
			title('High Power');
			xlabel('Trial #');

			metVal = zeros(128,1);
			for trialN = (1+filtWidth):(128-filtWidth)
				trialBot = trialN - filtWidth;
				trialTop = trialN + filtWidth;
				ix = find((PTgenoList == genoN) & (PTpowerNList <= 2) &...
						  (PTtrialNList >= trialBot) & (PTtrialNList <= trialTop));
				metVal(trialN) = nanmean(scoresByTrack(ix,metricN));
			end

			FF = ffsubplot(2,1,2);
			plot(filtRange, metVal(filtRange),'.-','Color',pretty(genoNn)); hold on;
			title('Low power');
			xlabel('Trial #');
		end

		FF.setTitle(metricLabels{metricN});
		FF.PDF(['KenTimeSer-',num2str(metricN),'.pdf']);
	end

			

			

