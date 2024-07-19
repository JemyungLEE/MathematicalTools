package probability.distribution;

public class Normal {

	public double mean;
	public double variance;
	public double standardDeviation; //standard deviation
	
	public Normal(){}
	
	public Normal(double mean, double std){
		initiate(mean, std);
	}
	
	public void initiate(double mean, double std){
		this.mean = mean;		
		this.standardDeviation = std;
		this.variance = Math.pow(std,2);
	}
	
	public double getMean(){
		return this.mean;
	}
	
	public double getSTD(){
		return this.standardDeviation;
	}
	
	public double getVariance(){
		return this.variance;
	}
	
	public void setMean(double mean){
		this.mean = mean;
	}
	
	public void setSTD(double std){
		this.standardDeviation = std;
	}
	
	public void setVariance(double var){
		this.variance = var;
	}
	
	public double getPDF(double x){
		//get probability density function
		
		return Math.exp(-(Math.pow((x-mean),2))/(2.0*variance))
										/Math.sqrt(2*Math.PI*variance);
	}
	
	public double getCDF(double x){
		//get cumulative density function
		
		int i;
		double cdf = 0.0;
		
	    double [] a = { 1.24818987e-4, -1.075204047e-3, 5.198775019e-3,
			     -1.9198292004e-2, 5.9054035642e-2, -1.51968751364e-1,
			     3.19152932694e-1, -5.319230073e-1, 7.97884560593e-1  };
		double [] b = { -4.5255659e-5, 1.5252929e-4, -1.9538132e-5, 
			      -6.76904986e-4, 1.390604284e-3, -7.9462082e-4,
			      -2.034254874e-3, 6.549791214e-3, -1.0557625006e-2,
			      1.1630447319e-2, -9.279453341e-3, 5.353579108e-3,
			      -2.141268741e-3, 5.35310849e-4, 9.99936657524e-1 };
		
		double z, w, y, pz;

		z = (x - mean)/standardDeviation ;

		if(z != 0.0) {
			y = Math.abs(z)*0.5;
			if(y < 3.0) {
				if(y < 1.0) {
					w = y*y;
					pz = a[0];
					for(i=1; i<9 ; i++) pz = pz*w+a[i];
					pz *= y*2.0;
			      	}
				else {
					y -= 2.0;
					pz = b[0];
					for(i=1 ; i<15 ; i++) pz = pz * y + b[i];
				}
			}
			else pz = 1.0;
		}
		else pz = 0.0;
			  
		if(z<0.0) cdf = (1.0-pz)*0.5;
		else cdf = (1.0 + pz)*0.5;
		
		return cdf;		
	}
	
	public static void main(String[] args) {

		double mean=3;
		double std=1;
		double x=8;
		Normal nm = new Normal(mean,std);
		System.out.println("<normal distribution>");
		System.out.println("mean:"+nm.mean+"\tstd:"+nm.standardDeviation);
		System.out.printf("pdf at %4.3f: %6.5f",x,nm.getPDF(x));
		System.out.println();
		System.out.printf("cdf at %4.3f: %6.5f",x,nm.getCDF(x));
	}

}
