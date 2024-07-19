/**
 *  Subject: Power Law Distribution
 *  Developer: Jemyung Lee
 *  Developed Data: 2012.8.1
 *  Last Modified Data: 2012.8.2 
 *  Department: Seoul Nat. Univ. depart. of Rural Systems Engineering
 *  Description: 
 */

package probability.distribution;

public class PowerLaw {

	double scale;
	double normConst; 
	double max;
	double interval;
	
	public PowerLaw(){}
	
	public PowerLaw(double scale){
		initiate(scale);
		setNormConst();
	}
	
	public void initiate(double scale){
		this.scale = scale;		
		this.max = 100;
		this.interval = 0.1;
	}
	
	public double getScale(){
		return this.scale;
	}
	
	public double getMax(){
		return this.max;
	}
	
	public double getInterval(){
		return this.interval;
	}
	
	public double getNormConst(){
		return this.normConst;
	}
	
	public void setNormConst(){
		this.normConst = 0;
		for(double i=1 ; i<max+this.interval ; i+=this.interval)
			this.normConst += Math.pow(i, -1.0*this.scale)*this.interval;
	}
	
	public void setScale(double scale){
		this.scale = scale;
	}
	
	public double getPDF(double x){
		//get probability density function
		
		if(x<1) return 0.0;
		else return Math.pow(x,-1.0*this.scale)/this.normConst;
	}
	
	public double getCDF(double x){
		//get cumulative density function
		
		double i;
		double cdf = 0.0;
		
		for(i=1 ; i<x+this.interval ; i+=this.interval) 
			cdf += this.getPDF(i)*this.interval;
				
		return cdf;		
	}
	
	public static void main(String[] args) {
		
		double scale=1.1;
		double x=5;
		PowerLaw pl = new PowerLaw(scale);
		System.out.println("<power law distribution>");
		System.out.println("scale: "+pl.scale);
		System.out.println("normal constant: "+pl.normConst);
		System.out.printf("pdf at %4.3f: %6.5f",x,pl.getPDF(x));
		System.out.println();
		System.out.printf("cdf at %4.3f: %6.5f",x,pl.getCDF(x));
	}

}
