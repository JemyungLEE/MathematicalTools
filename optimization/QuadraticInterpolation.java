package optimization;

public class QuadraticInterpolation {

	/**
	 *  Subject: Numerical Methods in Environmental Engineering
	 *  Name: Jemyung Lee
	 *  ID: 2008-30334
	 *  
	 *  Homework Number: 2
	 *  Description: Quadratic Interpolation Method
	 *  			(One-dimensional unconstrained optimization)
	 */
	
	public QuadraticInterpolation(int iter, double x0, double x1, double x2){
		algorithm(iter, x0, x1, x2);
	}
	
	public double function(double x){
		//definition of a function
		double fx;
		
		fx = 2.0*Math.sin(x)-Math.pow(x, 2)/10.0;
		
		return fx;
	}

	public void algorithm(int iter, double x0, double x1, double x2){
		//variables
		double x3;
		double fx0, fx1, fx2, fx3;
		//time variable
		long startTime, endTime;
		double operationTime;
		
		//printing variables
		System.out.printf(" i %8s %8s %8s %8s %8s %8s %8s %8s   time"
				,"x0","f(x0)","x1","f(x1)","x2","f(x2)","x3","f(x3)");
		System.out.println();

		//checking start time
		startTime = System.currentTimeMillis();
		
		//operation
		for(int i=0 ; i<iter ; i++){
			//values of guesses
			fx0 = function(x0);
			fx1 = function(x1);
			fx2 = function(x2);
			
			//maximum point x
			x3 = (fx0*(x1*x1-x2*x2)+fx1*(x2*x2-x0*x0)+fx2*(x0*x0-x1*x1))/
					(2*(fx0*(x1-x2)+fx1*(x2-x0)+fx2*(x0-x1)));
			//maximum value f(x)
			fx3 = function(x3);
			
			//checking operation time
			endTime = System.currentTimeMillis();
			operationTime = (double)(endTime - startTime)/1000.0;
			
			//print values
			System.out.printf("%2d %8.6f %8.6f %8.6f %8.6f %8.6f %8.6f"+
					" %8.6f %8.6f %6.3f"
					,i+1,x0,fx0,x1,fx1,x2,fx2,x3,fx3,operationTime);
			System.out.println();
			
			//next guesses
			if(x3>x1){
				x0 = x1;
				x1 = x3;
			}
			else if(x3<x1){
				x2 = x1;
				x1 = x3;
			}
		}
	}
	
	public static void main(String[] args) {
		
		new QuadraticInterpolation(20, 0, 1, 4);

	}
}
