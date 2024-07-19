package optimization;

public class NewtonMethod {

	/**
	 *  Subject: Numerical Methods in Environmental Engineering
	 *  Name: Jemyung Lee
	 *  ID: 2008-30334
	 *  
	 *  Homework Number: 2
	 *  Description: Newton's Method
	 *  			(One-dimensional unconstrained optimization)
	 */
	
	public NewtonMethod(int iter, double x){
		algorithm(iter, x);
	}
	
	public double function(double x){
		//definition of a function
		double fx;
		
		fx = 2.0*Math.sin(x)-Math.pow(x, 2)/10.0;
		
		return fx;
	}

	public double fdFunction(double x){
		//1st derivative of the function
		
		return 2.0*Math.cos(x)-x/5.0;
	}
	
	public double sdFunction(double x){
		//2nd derivatives of the function
		
		return -2.0*Math.sin(x) - 1.0/5.0;
	}
	
	public void algorithm(int iter, double x){
		//variables
		double fx, fdfx, sdfx;
		//time variable
		long startTime, endTime;
		double operationTime;
		
		//printing variables
		System.out.printf(" i %8s %8s %8s %8s","x","f(x)","f'(x)","f''(x)");
		System.out.println();

		//checking start time
		startTime = System.currentTimeMillis();
		
		//operation
		for(int i=0 ; i<iter ; i++){
			fx = function(x);
			fdfx = fdFunction(x);
			sdfx = sdFunction(x);
			
			//checking operation time
			endTime = System.currentTimeMillis();
			operationTime = (double)(endTime - startTime)/1000.0;
			
			//print values
			System.out.printf("%2d %8.6f %8.6f %8.6f %8.6f %8.3f"
					,i+1,x,fx,fdfx,sdfx, operationTime);
			System.out.println();
			
			//next guess
			x = x - fdfx/sdfx;			
		}
	}
		
	public static void main(String[] args) {
			
		new NewtonMethod(20, 2.5);

	}
}
