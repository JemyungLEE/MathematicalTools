package probability.distribution;

public class Frequency {

	public int population;
	public int intervalSize;
	public double intervalSpan;
	public double min, max;
	public double[] interval;	//[0]:min, [size]:max
	public int[] frequency;	//[0]: low than min.,[size+1]: upper than max.
	public double[] distribution; //[0]: low than min.,[size+1]: upper than max.
	
	public Frequency(){}
	
	public Frequency(int size){
		initiate(size);
	}
	
	public Frequency(int size, double min, double max){
		initiate(size, min, max);
	}
	
	public void initiate(int size){		
		int i;
		intervalSize = size;
		interval = new double[size+1];
		frequency = new int[size+2];
		distribution = new double[size+2];
		for(i=0 ; i<(size+1) ; i++)	
			interval[i] = 0.0;
		for(i=0 ; i<(size+2) ; i++){
			frequency[i] = 0;
			distribution[i] = 0.0;
		}
		
	}
	
	public void initiate(int size, double min, double max){
		initiate(size);
		setInterval(size, min, max);
	}
	
	public void setInterval(int size){
		intervalSpan = (max-min)/size;
		for(int i=0 ; i<=size ; i++) interval[i] = min+intervalSpan*i;		
	}
	
	public void setInterval(int size, double min, double max){
		this.min = min;
		this.max = max;
		setInterval(size);	
	}
	
	public void normalize(){
		for(int i=0 ; i<intervalSize+2 ; i++) 
			distribution[i] = (double) frequency[i]/population/intervalSpan;
	}
	
	public double getPdf(double x){
		int i;
		
		if(x<min) return distribution[0];
		else if(x>=max) return distribution[intervalSize+1];
		else
			for(i=0 ; i<intervalSize ; i++)
				if(x>=interval[i] && x<interval[i+1]) return distribution[i+1];
		
		return -1; //return error.
	}
	
	public double getCdf(double x){
		int i;
		double cdf = 0.0;
		
		cdf += distribution[0]*intervalSpan;
		if(x<min) return cdf;
		else if(x>=max) return 1.0;
		else{
			for(i=0 ; i<intervalSize ; i++){
				cdf += distribution[i+1]*intervalSpan;
				if(x>=interval[i] && x<interval[i+1]) return cdf;
			}
		}
		return -1; //return error.
	}
}
