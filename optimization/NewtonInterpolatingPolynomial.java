package optimization;

import java.io.File;
import java.io.IOException;
import java.util.Scanner;

public class NewtonInterpolatingPolynomial {

	/**
	 *  Subject: Numerical Methods for Engineers
	 *  Name: Jemyung Lee
	 *  ID: 2008-30334
	 *  
	 *  Homework Number: 4
	 *  Description: Newton's Interpolating Polynomials
	 */
	
	public NewtonInterpolatingPolynomial(String source, int order){
		
		algorithm(source, order);
	}
	
	public void algorithm(String source, int order){
		int i,j;
		int dataN=0;
		double xCal=0;
		double[] x = new double[dataN];
		double[] fx = new double[dataN];
		double[][] fdd;		
		double xterm;
		double[] y;
		double[] ea;
		
		/*** read data from file ***/
		try{
			File file = new File(source);
			Scanner scan = new Scanner(file);
			
			
			//read calculate order
			scan.next();
			dataN = scan.nextInt();
			
			//read calculate point x
			scan.next();
			xCal = scan.nextDouble();
			
			//read data			
			for(i=0 ; i<3 ; i++) scan.next();
			x = new double[dataN];
			fx = new double[dataN];
			for(i=0 ; i<dataN ; i++){
				x[i] = scan.nextDouble();
				fx[i] = scan.nextDouble();
			}
			
			scan.close();			
		} catch(IOException e) {}
		
		
		/*** calculation part ***/

		//set variables
		fdd = new double[order][order];
		y = new double[order];
		ea = new double[order-1];
		
		//interpolating polynomials
		for(i=0 ; i<order ; i++){
			fdd[i][0] = fx[i];
		}
		for(j=1 ; j<order ; j++){
			for(i=0 ; i<order-j ; i++){
				fdd[i][j] = (fdd[i+1][j-1] - fdd[i][j-1])/(x[i+j]-x[i]); 
			}
		}
	
		//calculate f_n(x) and error
		xterm = 1;
		y[0] = fdd[0][0];		
		for(i=1 ; i<order ; i++){
			xterm *= xCal - x[i-1];
			y[i] = y[i-1] + fdd[0][i] * xterm;
			ea[i-1] = y[i] - y[i-1];
		}
		
		
		/*** printing part ***/
		
		//print the inputs
		System.out.println("data sets: "+dataN);
		System.out.println("interpolating x: "+xCal);
		System.out.println("Data");
		System.out.println("x\tf(x)");
		for(i=0 ; i<dataN ; i++){
			System.out.println(x[i]+"\t"+fx[i]);
		}
		
		//print the results
		System.out.println("order\tf("+xCal+")\t\terror");
		for(i=0 ; i<order ; i++){
			System.out.printf("%d\t%.6f",i+1,y[i]);
			if(i<order-1){
				System.out.printf("\t%.6f",ea[i]);
				System.out.println();
			}
		}
	}
	
	
	public static void main(String[] args) {
		String source;
		source = "c:/optimization/interpolatingInput2.txt";
		
		new NewtonInterpolatingPolynomial(source, 6);
	}

}
