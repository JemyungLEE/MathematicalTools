package optimization.function;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Scanner;

import Jama.Matrix;

public class NewmarkMethod {
	/**
	 *  Subject: Structural Optimization
	 *  Developer: Jemyung Lee (ID: 2008-30334)
	 *    
	 *  Description: Newmark's method to solve second differential equation
	 */
	
	public Matrix readMatrix(String filename){
		int i,j;
		int size = 9;
		double[][] m = new double[size][size];
		
		try{
			File file = new File(filename);
			Scanner scan = new Scanner(file);
						
			for(i=0 ; i<size ; i++)
				for(j=0 ; j<size ; j++) m[i][j] = scan.nextDouble();	
			
			scan.close();			
		} catch(IOException e) {}
		
		return new Matrix(m);
	}
	
	public Matrix readEarthquake(String filename){
		int size = 2502;
		double[][] eq = new double[size][2];
		
		try{
			File file = new File(filename);
			Scanner scan = new Scanner(file);
						
			for(int i=0 ; i<size ; i++){
				eq[i][0] = scan.nextDouble();
				eq[i][1] = scan.nextDouble();
			}
			
			scan.close();			
		} catch(IOException e) {}
		
		return new Matrix(eq);
	}
	
	public Matrix makeKmatrix(double[] k){
		int i,j;
		int size = 9;
		double[][] kmatrix = new double[size][size];
		
		for(i=0 ; i<size ; i++) for(j=0 ; j<size ; j++) kmatrix[i][j] = 0.0;
			
		kmatrix[0][0] = k[0]+k[1];
		kmatrix[0][1] = -k[1];
		for(i=1 ; i<size-1 ; i++){
			kmatrix[i][i-1] = -k[i];
			kmatrix[i][i] = k[i]+k[i+1];
			kmatrix[i][i+1] = -k[i+1];
		}
		kmatrix[size-1][size-2] = -k[size-1];
		kmatrix[size-1][size-1] = k[size-1];
		
		return new Matrix(kmatrix);
	}
	
	public double process(double alpha, 
						Matrix M, Matrix C, Matrix eqe, double[] k){
		int i,j;
		int ndof = 9;
		double gam = 0.5;
		double beta = 0.25;
		double delt = 0.02;
		double[][] ones = new double[ndof][1];
		double[][] zeros = new double[ndof][1];
		
		double K0 = 17.7*Math.pow(10, 8);
		double D0 = 0.4188;
		double sumK=0, sumD=0, maxD;
		double[][] d = new double[ndof][eqe.getRowDimension()-1];
		
		for(i=0 ; i<ndof ; i++){
			ones[i][0] = 1.0;
			zeros[i][0] = 0.0;
		}
		
		Matrix K = makeKmatrix(k);
		Matrix b = M.times(-1.0).times(new Matrix(ones));
		Matrix P0 = b.times(eqe.get(0, 0));
		Matrix x = new Matrix(zeros);
		Matrix xd = new Matrix(zeros);
		Matrix tmpMatrix = P0.minus(C.times(xd)).minus(K.times(x));
		Matrix xdd = M.inverse().times(tmpMatrix);
		Matrix Kh = K.plus(C.times(gam/(beta*delt)))
						.plus(M.times(1.0/(beta*delt*delt)));
		Matrix iKh = Kh.inverse();
		Matrix aa = M.times(1.0/(beta*delt)).plus(C.times(gam/beta));
		Matrix bb = M.times(1.0/(2.0*beta))
							.plus(C.times(delt*(gam/(2.0*beta)-1.0)));
		Matrix deqe = new Matrix(eqe.getRowDimension()-1,1);
		for(i=0 ; i<deqe.getRowDimension() ; i++)
			deqe.set(i, 0, eqe.get(i+1, 0)-eqe.get(i, 0));
		
		Matrix dPi;
		Matrix dPih;
		Matrix dxi;
		Matrix dxdi;
		Matrix dxddi;
		for(int ii=0 ; ii<deqe.getRowDimension() ; ii++){
			dPi = b.times(deqe.get(ii, 0));
			dPih = dPi.plus(aa.times(xd)).plus(bb.times(xdd));
			dxi = iKh.times(dPih);
			dxdi = dxi.times(gam/(beta*delt)).minus(xd.times(gam/beta))
								.plus(xdd.times(delt*(1.0-gam/(2.0*beta))));
			dxddi = dxi.times(1.0/(beta*delt*delt))
											.minus(xd.times(1.0/(beta*delt)))
											.minus(xdd.times(1.0/(2.0*beta)));
			x = x.plus(dxi);
			xd = xd.plus(dxdi);
			xdd = xdd.plus(dxddi);
			
			d[0][ii] = x.get(0, 0);
			for(i=1 ; i<ndof ; i++) d[i][ii]=Math.abs(x.get(i,0)-x.get(i-1,0)); 
		}		
		
		for(i=0 ; i<ndof ; i++){
			sumK += k[i];
			maxD = d[i][0];
			for(j=1 ; j<deqe.getRowDimension() ; j++)
				if(maxD<d[i][j]) maxD = d[i][j];
			sumD += maxD;
		}
		
		return alpha*sumK/K0+(1.0-alpha)*sumD/D0;
	}
	
	public static void main(String[] args) {
		NewmarkMethod nm = new NewmarkMethod();
		Matrix M = nm.readMatrix("c:/optimization/hw8/mass.txt");
		Matrix C = nm.readMatrix("c:/optimization/hw8/damping.txt");
		Matrix eqe = nm.readEarthquake("c:/optimization/hw8/elcentro.txt");
		double[] k = {2.1*Math.pow(10, 8),
				2.1*Math.pow(10, 8),
				2.1*Math.pow(10, 8),
				2.1*Math.pow(10, 8),
				2.1*Math.pow(10, 8),
				2.1*Math.pow(10, 8),
				1.7*Math.pow(10, 8),
				1.7*Math.pow(10, 8),
				1.7*Math.pow(10, 8)};
		System.out.println(nm.process(0.1, M, C, eqe, k));		
	}
	
	
}
