package optimization;

public class GoldenSectionSearch {

	/**
	 *  Subject: Numerical Methods for Engineers
	 *  Name: Jemyung Lee
	 *  ID: 2008-30334
	 *  
	 *  Homework Number: 1
	 */
	
	public GoldenSectionSearch(int iter, double xl, double xu, double es){
		algorithm(iter, xl, xu, es);
	}
	
	public double function(double x){
		double fx;
		
		fx = 2.0*Math.sin(x)-Math.pow(x, 2)/10.0;
		
		return fx;
	}
	
	public void algorithm(int iter, double xl, double xu, double es){
		double x1, x2;
		double fx1, fx2;
		double fxl, fxu;
		double d;
		double goldRatio = (Math.sqrt(5)-1)/2;
		//time variable
		long startTime, endTime;
		double operationTime;
		
		//printing variables
		System.out.printf(" i %8s %8s %8s %8s %8s %8s %8s %8s %8s"
				,"xl","f(xl)","x2","f(x2)","x1","f(x1)","xu","f(xu)","d");
		System.out.println();
		
		//checking start time
		startTime = System.currentTimeMillis();
		
		//initiate
		d = goldRatio * (xu-xl);		
		x1 = xl + d;
		x2 = xu - d;

		fxl = function(xl);
		fxu = function(xu);
		fx1 = function(x1);
		fx2 = function(x2);
		
		//checking operation time
		endTime = System.currentTimeMillis();
		operationTime = (double)(endTime - startTime)/1000.0;
		
		System.out.printf(" 1 %8.6f %8.6f %8.6f %8.6f %8.6f"+
				" %8.6f %8.6f %8.6f %8.6f %8.3f"
				,xl,fxl,x2,fx2,x1,fx1,xu,fxu,d, operationTime);
		System.out.println();
		
		for(int i=1 ; i<iter ; i++){			
			if(fx1 > fx2){
				xl = x2;
				d = goldRatio * (xu-xl);
				x2 = x1;
				x1 = xl + d;
			}
			else if(fx2 > fx1){	
				xu = x1;
				d = goldRatio * (xu-xl);
				x1 = x2;
				x2 = xu - d;				
			}
			fxl = function(xl);
			fxu = function(xu);
			fx1 = function(x1);
			fx2 = function(x2);
			
			//checking operation time
			endTime = System.currentTimeMillis();
			operationTime = (double)(endTime - startTime)/1000.0;
			
			System.out.printf("%2d %8.6f %8.6f %8.6f %8.6f %8.6f"+
					" %8.6f %8.6f %8.6f %8.6f %8.3f",
					(i+1),xl,fxl,x2,fx2,x1,fx1,xu,fxu,d,operationTime);
			System.out.println();
		}
		
	}
	
	public static void main(String[] args) {
		
		new GoldenSectionSearch(40, 0, 4, 3);

	}

}
