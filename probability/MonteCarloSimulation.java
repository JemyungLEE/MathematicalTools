package probability;

import java.util.Random;
import probability.distribution.*;

public class MonteCarloSimulation {

	public MonteCarloSimulation(){}
	

	public void mergeDistributions(Frequency fqd, Normal nd1, Normal nd2,
															int iteration){
		int i,j;
		int intervalSize;
		double min, max; 
		double span = 5.0;
		double nd1Mean, nd2Mean;
		double nd1Std, nd2Std;
		double nd1Value, nd2Value, sum;
		Random rnd1, rnd2;
		
		//set properties		
		intervalSize = fqd.intervalSize;		
		nd1Mean = nd1.mean;
		nd2Mean = nd2.mean;
		nd1Std = nd1.standardDeviation;
		nd2Std = nd2.standardDeviation;
		
		//set boundaries
		min = nd1.mean - (nd1.standardDeviation * span) +
				nd2.mean - (nd2.standardDeviation * span);		
		max = nd1.mean + (nd1.standardDeviation * span) +
				nd2.mean + (nd2.standardDeviation * span);
	
		//set intervals
		fqd.setInterval(intervalSize, min, max);
		
		//set seeds
		rnd1 = new Random(System.currentTimeMillis());
		rnd2 = new Random((long)System.currentTimeMillis()/2);
		
		//process		
		for(i=0 ; i<iteration ; i++){
			nd1Value = rnd1.nextGaussian() * nd1Std + nd1Mean;
			nd2Value = rnd2.nextGaussian() * nd2Std + nd2Mean;
			sum = nd1Value + nd2Value;
			
			if(sum<min) fqd.frequency[0]++;
			else if(sum>=max) fqd.frequency[fqd.intervalSize+1]++; 
			else
				for(j=0 ; j<intervalSize ; j++)
					if(sum>=fqd.interval[j] && sum<fqd.interval[j+1])
						fqd.frequency[j+1]++;
		}
		
		//make distribution
		fqd.population = iteration;
		fqd.normalize();
	}
	
	public Frequency mergeDistributions(int intervalSize, 
									Normal nd1, Normal nd2,	int iteration){
		
		Frequency fqd = new Frequency(intervalSize);
		mergeDistributions(fqd, nd1, nd2, iteration);

		return fqd;
	}

	
	public static void main(String[] args) {
		int iteration = 1000000;
		int size = 1000;
		double mean1 = 10.0;
		double mean2 = 4.0;
		double std1 = 2.0;
		double std2 = 2.0;
		Frequency fqd = new Frequency();
		Normal nd1 = new Normal(mean1, std1);
		Normal nd2 = new Normal(mean2, std2);
		
		MonteCarloSimulation mcs = new MonteCarloSimulation();
		fqd = mcs.mergeDistributions(size, nd1, nd2, iteration);
		
		System.out.println("MCS iteration: "+iteration);
		System.out.println("interval\tfrequency\tdistribution");
		System.out.println("~\t"+fqd.frequency[0]+"\t"+fqd.distribution[0]);
		for(int i=1 ; i<size+2 ; i++)
			System.out.println(fqd.interval[i-1]+"\t"+fqd.frequency[i]+"\t"
								+fqd.distribution[i]);
		System.out.println("~\t"+fqd.frequency[size+1]+"\t"
								+fqd.distribution[size+1]);
		
		System.out.println("pdf: "+fqd.getPdf(14));
		System.out.println("cdf: "+fqd.getCdf(14));
	}

}
