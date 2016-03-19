function rawMetricAnalysis()

	load('~/RTFW/Code/kennedyAnalysis/PTkennedyAnalysis.mat');
	allGenos = unique(genoList);


	decPI = scoresByFly(:,2);
	% Find an ordering for the genoList
	for genoNn = 1:length(allGenos)
		genoN = allGenos(genoNn);

		ix = find((genoList == genoN) & (powerNList == 8));
		sortVals(genoNn) = nanmean(decPI(ix));
	end
	[B, genoSortOrder] = sort(sortVals,'descend');
	colorSpace = flipud(jet);
	colorIXList = round(linspace(1,64,length(allGenos)));

	for metricN = [3:14]

		ffsubplot(4,3,metricN - 2); 

		for genoNn = 1:length(allGenos)
			genoN = allGenos(genoSortOrder(genoNn));
			plotColor = colorSpace(colorIXList(genoNn),:);

			ix1 = find((powerNList == 1) & (genoList == genoN));
			ix8 = find((powerNList == 8) & (genoList == genoN));

			vals1 = scoresByFly(ix1, metricN);
			vals8 = scoresByFly(ix8, metricN);

			plot([1, 2],[nanmean(vals1(:)),nanmean(vals8(:))],'o-','Color',plotColor);
			xlim([.75 2.25]);
			ylabel(metricLabels{metricN});
			set(gca,'XTick',[1 2],'XTickLabel',[]);

			hold on;
		end

	end

