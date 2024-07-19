package optimization;

import java.io.File;
import java.io.IOException;
import java.util.Scanner;

public class LagrangeInterpolatingPolynomial {

	/**
	 *  Subject: Numerical Methods for Engineers
	 *  Name: Jemyung Lee
	 *  ID: 2008-30334
	 *  
	 *  Homework Number: 4
	 *  Description: Lagrange Interpolating Polynomials
	 */
	
	public LagrangeInterpolatingPolynomial(String source, int order){
		
		algorithm(source, order);
	}
	
	public void algorithm(String source, int order){
		int i,j;
		int dataN=0;
		double xCal=0;
		double[] x = new double[order];
		double[] fx = new double[order];
		double sum;
		double product;
		double lagrange;
		
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
		sum = 0;
		
		for(i=0 ; i<order+1 ; i++){
			product = fx[i];
			for(j=0 ; j<order+1 ; j++){
				if(i!=j){
					product *= (xCal - x[j])/(x[i]-x[j]);
				}
			}
			sum += product;
		}
		lagrange = sum;
		
		
		/*** printing part ***/
		
		//print the inputs
		System.out.println("data sets: "+dataN);
		System.out.println("interpolating x: "+xCal);
		System.out.println("Data");
		System.out.println("x\tf(x)");
		for(i=0 ; i<dataN ; i++){
			System.out.println(x[i]+"\t"+fx[i]);
		}
		System.out.println();
		
		//print the results
		System.out.printf("order: %d",order);
		System.out.println();
		System.out.printf("lagrange: %.6f",lagrange);
		System.out.println();
	}
	
	public static void main(String[] args) {
		String source;
		source = "c:/optimization/interpolatingInput2.txt";
		
		new LagrangeInterpolatingPolynomial(source, 1);
	}

}
