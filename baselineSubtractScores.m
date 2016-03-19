function blScores = baselineSubtractScores(scoresByFly, powerNList, genoList)

	allGenos = unique(genoList);

	for genoNn = 1:length(allGenos)
		genoN = allGenos(genoNn);

		ix1 = find((powerNList == 1) & (genoList == genoN));
		ix2 = find(genoList == genoN);
		meanVals = nanmean(scoresByFly(ix1,:),1);

		blScores(ix2,:) = scoresByFly(ix2,:) - ones(length(ix2(:)),1)*meanVals;
	end




