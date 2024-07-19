/**
 *  Subject: Convolition
 *  Developer: Jemyung Lee
 *  Developed Data: 2010. 8. 30
 *  Last Modified Data: 2012. 8. 2 
 *  Department: Seoul Nat. Univ. depart. of Rural Systems Engineering
 *  Description: 
 */

package probability;

import probability.distribution.*;


public class Convolution {

	public Convolution(){
		
	}

	public Convolution( Normal nd1, Normal nd2){
	}
		
	public Convolution(PowerLaw pl1, PowerLaw pl2){
	}
	
	public double getPdf(double x, Normal nd1, Normal nd2){
		double i;
		double delta = 0.0001, interval;
		double pdf = 0.0;
		double lowBoundary, upperBoundary, tmpBoundary, span;
		
		//set boundaries
		span = 5.0;
		lowBoundary = nd1.mean - (nd1.standardDeviation * span);
		tmpBoundary = nd2.mean - (nd2.standardDeviation * span);
		if(lowBoundary > tmpBoundary) lowBoundary = tmpBoundary;
		upperBoundary = nd1.mean + (nd1.standardDeviation * span);
		tmpBoundary = nd2.mean + (nd2.standardDeviation * span);
		if(upperBoundary < tmpBoundary) upperBoundary = tmpBoundary;
		
		//calculate pdf
		if(x < lowBoundary || x > upperBoundary) return pdf;
		else{
			interval = delta * (upperBoundary - lowBoundary); 
			for(i=lowBoundary ; i<upperBoundary ; i+=interval){
				pdf += nd1.getPDF(i) * nd2.getPDF(x-i) * interval;
			}
			return pdf;
		}		
		
	}
	
	public double getCdf(double x, Normal nd1, Normal nd2){
		double i;
		double delta = 0.0001;
		double cdf = 0.0;
		double lowBoundary, upperBoundary, tmpBoundary, span;
		
		//set boundaries
		span = 5.0;
		lowBoundary = nd1.mean - (nd1.standardDeviation * span);
		tmpBoundary = nd2.mean - (nd2.standardDeviation * span);
		if(lowBoundary > tmpBoundary) lowBoundary = tmpBoundary;
		upperBoundary = nd1.mean + (nd1.standardDeviation * span);
		tmpBoundary = nd2.mean + (nd2.standardDeviation * span);
		if(upperBoundary < tmpBoundary) upperBoundary = tmpBoundary;

		//calculate cdf
		if(x < lowBoundary) return 0.0;
		else if( x > upperBoundary) return 1.0;
		else{
			delta *= (upperBoundary - lowBoundary); 
			for(i=lowBoundary ; i<upperBoundary ; i+=delta){
				cdf += getPdf(x,nd1,nd2) * delta;
			}
			return cdf;
		}
	}
	
	public double getPdf(double x, PowerLaw pl1, PowerLaw pl2){
		double i;
		double max;
		double interval;
		double pdf = 0.0;
				
		//set boundary
		if(pl1.getMax() > pl2.getMax()) max = pl1.getMax();
		else max = pl2.getMax();
				
		if(pl1.getInterval() < pl2.getInterval()) interval = pl1.getInterval();
		else interval = pl2.getInterval();
		
		//calculate pdf		
		for(i=-max ; i<max ; i+=interval)
			pdf += pl1.getPDF(i) * pl2.getPDF(x-i) * interval;
		
		return pdf;
	}
	
	public double get2ndPdf(double x, PowerLaw pl1, PowerLaw pl2){
		double i;
		double max;
		double interval;
		double pdf = 0.0;
				
		//set boundary
		if(pl1.getMax() > pl2.getMax()) max = pl1.getMax();
		else max = pl2.getMax();
				
		if(pl1.getInterval() < pl2.getInterval()) interval = pl1.getInterval();
		else interval = pl2.getInterval();
		
		//calculate pdf		
		for(i=-max ; i<max ; i+=interval)
			pdf += this.getPdf(i, pl1, pl2) * pl2.getPDF(x-i) * interval;
		
		return pdf;
	}
	
	public double get3rdPdf(double x, PowerLaw pl1, PowerLaw pl2){
		double i;
		double max;
		double interval;
		double pdf = 0.0;
				
		//set boundary
		if(pl1.getMax() > pl2.getMax()) max = pl1.getMax();
		else max = pl2.getMax();
				
		if(pl1.getInterval() < pl2.getInterval()) interval = pl1.getInterval();
		else interval = pl2.getInterval();
		
		//calculate pdf		
		for(i=-max ; i<max ; i+=interval)
			pdf += this.get2ndPdf(i, pl1, pl2) * pl2.getPDF(x-i) * interval;
		
		return pdf;
	}
	
	public double getPdf(double x, double gap, PowerLaw pl1, PowerLaw pl2){
		double i;
		double max;
		double interval;
		double pdf = 0.0;
				
		//set boundary
		if(pl1.getMax() > pl2.getMax()) max = pl1.getMax();
		else max = pl2.getMax();
				
		if(pl1.getInterval() < pl2.getInterval()) interval = pl1.getInterval();
		else interval = pl2.getInterval();
		
		//calculate pdf		
		for(i=-max ; i<max ; i+=interval)
			pdf += pl1.getPDF(i) * pl2.getPDF(x-i-gap) * interval;
		
		return pdf;
	}
	
	public double getCdf(double x, PowerLaw pl1, PowerLaw pl2){
		double i;
		double max;
		double interval;
		double cdf = 0.0;
		
		//set boundary
		if(pl1.getMax() > pl2.getMax()) max = pl1.getMax();
		else max = pl2.getMax();		
		
		if(pl1.getInterval() < pl2.getInterval()) interval = pl1.getInterval();
		else interval = pl2.getInterval();

		//calculate cdf
		for(i=-max ; i<=x ; i+=interval)
			cdf += this.getPdf(i, pl1, pl2) * interval;
		
		return cdf;

	}
	
	public static void main(String[] args) {
		
		double i,j;
		/*
		Normal nd1 = new Normal(3,1);
		Normal nd2 = new Normal(3,1);
		Convolution cv = new Convolution();
		System.out.println(cv.getCdf(6, nd1, nd2));
		*/
		//System.out.println(cv.getPdf(3, nd1, nd2));
		
		
		PowerLaw pl1 = new PowerLaw(1.1);
		PowerLaw pl2 = new PowerLaw(1.1);
		Convolution cv = new Convolution();
		
		System.out.println("i\tpdf\t2nd_pdf\t3rd_pdf");		
		for(i = 0 ; i< 101 ; i++){
			System.out.print(i+"\t");			
			System.out.print(cv.getPdf(i, pl1, pl2)+"\t");
			System.out.print(cv.get2ndPdf(i, pl1, pl2)+"\t");
			System.out.print(cv.get3rdPdf(i, pl1, pl2)+"\t");
			System.out.println();
		}
				
		
		/*
		System.out.println("i\tpdf\tcdf");		
		for(i = 0 ; i< 200 ; i++){
			System.out.print(i+"\t");	
			for(j=0 ; j<100 ; j+=20){			
				System.out.print(cv.getPdf(i, j, pl1, pl2)+"\t");
				//System.out.println(cv.getCdf(i, pl1, pl2));
			
			}	
			System.out.println();
		}
				
		for(i = 0 ; i< 100 ; i++){
			System.out.print(i+"\t");
			System.out.print(cv.getPdf(i, pl1, pl2)+"\t");
		}
		*/		
	}
	
}