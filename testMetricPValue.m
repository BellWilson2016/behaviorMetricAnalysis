function testMetricPValue()

	load('~/RTFW/Code/kennedyAnalysis/PTkennedyAnalysis83bTarget.mat');

	testMetricList = [11,12,13,14];

	disp('Uncorrected P');
	for metricNn = 1:length(testMetricList)
		metricN = testMetricList(metricNn);

		ix = find((genoList == 18) & (powerNList == 1));
		ctrlVals = scoresByFly(ix,metricN);
		ix = find((genoList == 18) & (powerNList == 8));
		maxVals  = scoresByFly(ix,metricN);
	
		pVals(metricNn) = bootStrapP(maxVals, ctrlVals, 10000);
		disp([metricLabels{metricN},'  p = ',num2str(pVals(metricNn))]);
	end	

	pAdj = bonferroniHolm(pVals);
	disp('Corrected P');
	for metricNn = 1:length(testMetricList)
		metricN = testMetricList(metricNn);
		disp([metricLabels{metricN},'  p = ',num2str(pAdj(metricNn))]);
	end

					% Light ,  Ctrl
function p = bootStrapP(met1, met2, nBoots)

	totalMet = cat(1,met1(:),met2(:));
	bootStats = zeros(nBoots,1);
	testStat = testStatistic(met1, met2);
	for bootN = 1:nBoots
		bMet1 = totalMet(randi(length(totalMet),length(met1),1));
		bMet2 = totalMet(randi(length(totalMet),length(met2),1));
		bootStats(bootN) = testStatistic(bMet1,bMet2);
	end

	p = (sum(bootStats >= testStat) + 1)/(nBoots + 1);

					  % Light,  Ctrl
function out = testStatistic(met1, met2)

	out = abs(nanmean(met1 - nanmean(met2)));



