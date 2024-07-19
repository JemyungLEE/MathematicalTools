package optimization;

import java.io.File;
import java.io.IOException;
import java.util.Scanner;

import Jama.*;

public class InteriorPointMethod {

	Matrix A;
	Matrix b;
	Matrix c;
	Matrix x0;	//initial value
	Matrix x1;
	
	Matrix D;
	Matrix B;
	Matrix p;
	Matrix y0;
	Matrix y1;
	
	double Z0;
	double Z1;
	
	Matrix P;
	
	Matrix d0;	//direction vector
	double ramda;
	double ratio=0.9;	//ramda:percentage of max.
	
	Matrix w0;	//value of dual problem
	
	double dualityGap;
	double error=0.01;
	
	public void readInputData(String inputFile){
		int i,j;
		int row, column;
		double[][] temp;		
		
		try{
			File file = new File(inputFile);
			Scanner scan = new Scanner(file);
			
			//read matrix A
			scan.next();
			row = scan.nextInt();
			column = scan.nextInt();
			temp = new double[row][column];
			for(i=0 ; i<row ; i++)
				for(j=0 ; j<column ; j++)
					temp[i][j] = scan.nextDouble();
			A = new Matrix(temp);
			
			//read vector b
			scan.next();
			row = scan.nextInt();
			temp = new double[row][1];
			for(i=0 ; i<row ; i++)
				temp[i][0] = scan.nextDouble();
			b = new Matrix(temp);
			
			//read vector c
			scan.next();
			row = scan.nextInt();
			temp = new double[row][1];
			for(i=0 ; i<row ; i++)
				temp[i][0] = scan.nextDouble();
			c = new Matrix(temp);
			
			//read initial point x0
			scan.next();
			row = scan.nextInt();
			temp = new double[row][1];
			for(i=0 ; i<row ; i++)
				temp[i][0] = scan.nextDouble();
			x0 = new Matrix(temp);			
			
			scan.close();			
		} catch(IOException e){}
		

	}

	public void setInputData(){
		
		double[][] inputA = {{2,1,1,0},{0,1,0,1}};
		double[][] inputB = {{24},{8}};
		double[][] inputC = {{-1},{-1},{0},{0}};
		double[][] inputX0 = {{2},{4},{16},{4}};
		A = new Matrix(inputA);
		b = new Matrix(inputB);
		c = new Matrix(inputC);
		x0 = new Matrix(inputX0);
	}	
	
	public void initiate(){
		int i;
		int size;
		
		//define D
		size=x0.getRowDimension();
		D = new Matrix(size,size);
		for(i=0 ; i<size ; i++) D.set(i,i,x0.get(i,0));

		//define y0
		y0 = D.inverse().times(x0);		
	}
	
	public void calculate(){
		int i,j;
		int size=x0.getRowDimension();;
		double maxRamda;
		int chk;
		double temp;
				
		//calculate B
		B = A.times(D);
		
		//calculate p
		p = D.times(c);
		
		//calculate P
		Matrix I = new Matrix(size,size);
		Matrix BT = B.transpose();
		for(i=0 ; i<size ; i++) I.set(i,i,1.0);
		P = I.minus(BT.times((B.times(BT)).inverse()).times(B));
			
		//calculate d0
		d0 = (P.times(p)).times(-1.0);
				
		//calculate ramda
		chk=0;
		maxRamda=0;
		for(i=0 ; i<size ; i++){
			temp = d0.get(i,0);
			if(temp<0){				
				temp = 1.0/(-1.0*temp);
				if(chk==0){
					maxRamda = temp;
					chk=1;					
				}
				else if(maxRamda>temp){
					maxRamda = temp;
				}
			}			
		}
		if(chk==0) System.err.println("LP problrm is unbounded");
		else ramda = ratio*maxRamda;
		
		//calculate y1
		y1 = y0.plus(d0.times(ramda));
		
		//calculate x1
		x1 = D.times(y1);
		
		//calculate Z0, Z1
		Z0 = c.transpose().times(x0).get(0,0);
		Z1 = c.transpose().times(x1).get(0,0);
	}
	
	public void checkDualityGap(){
		//calculate w0
		Matrix AT = A.transpose();
		Matrix D2 = D.times(D);
		w0 = (((A.times(D2.times(AT)).inverse()).times(A)).times(D2)).times(c);
				
		//calculate duality gap
		Matrix cT = c.transpose();
		Matrix bT = b.transpose();
		Matrix gap = cT.times(x0).minus(bT.times(w0));
		
		dualityGap = Math.abs(gap.get(0,0));
	}

	public void printResult(){
		int i,j;
	
		/***
		System.out.println("A");
		for(i=0 ; i<A.getRowDimension() ; i++){
			for(j=0 ; j<A.getColumnDimension() ; j++){
				System.out.print(""+A.get(i,j)+"\t");
			}
			System.out.println();
		}
		
		System.out.println("b");
		for(i=0 ; i<b.getRowDimension() ; i++){
			for(j=0 ; j<b.getColumnDimension() ; j++){
				System.out.print(""+b.get(i,j)+"\t");
			}
			System.out.println();
		}
		
		System.out.println("c");
		for(i=0 ; i<c.getRowDimension() ; i++){
			for(j=0 ; j<c.getColumnDimension() ; j++){
				System.out.print(""+c.get(i,j)+"\t");
			}
			System.out.println();
		}
		
		System.out.println("x0");
		for(i=0 ; i<x0.getRowDimension() ; i++){
			for(j=0 ; j<x0.getColumnDimension() ; j++){
				System.out.print(""+x0.get(i,j)+"\t");
			}
			System.out.println();
		}
		
		System.out.println("D");
		for(i=0 ; i<D.getRowDimension() ; i++){
			for(j=0 ; j<D.getColumnDimension() ; j++){
				System.out.print(""+D.get(i,j)+"\t");
			}
			System.out.println();
		}
		
		System.out.println("B");
		for(i=0 ; i<B.getRowDimension() ; i++){
			for(j=0 ; j<B.getColumnDimension() ; j++){
				System.out.print(""+B.get(i,j)+"\t");
			}
			System.out.println();
		}
		
		System.out.println("p");
		for(i=0 ; i<p.getRowDimension() ; i++){
			for(j=0 ; j<p.getColumnDimension() ; j++){
				System.out.print(""+p.get(i,j)+"\t");
			}
			System.out.println();
		}
		
		System.out.println("y0");
		for(i=0 ; i<y0.getRowDimension() ; i++){
			for(j=0 ; j<y0.getColumnDimension() ; j++){
				System.out.print(""+y0.get(i,j)+"\t");
			}
			System.out.println();
		}
		
		System.out.println("P");
		for(i=0 ; i<P.getRowDimension() ; i++){
			for(j=0 ; j<P.getColumnDimension() ; j++){
				System.out.print(""+P.get(i,j)+"\t");
			}
			System.out.println();
		}
		
		System.out.println("d0");
		for(i=0 ; i<d0.getRowDimension() ; i++){
			for(j=0 ; j<d0.getColumnDimension() ; j++){
				System.out.print(""+d0.get(i,j)+"\t");
			}
			System.out.println();
		}
		
		System.out.println("y1");
		for(i=0 ; i<y1.getRowDimension() ; i++){
			for(j=0 ; j<y1.getColumnDimension() ; j++){
				System.out.print(""+y1.get(i,j)+"\t");
			}
			System.out.println();
		}
		
		***/
		
		System.out.println("x1");
		for(i=0 ; i<x1.getRowDimension() ; i++){
			for(j=0 ; j<x1.getColumnDimension() ; j++){
				System.out.print(""+x1.get(i,j)+"\t");
			}
			System.out.println();
		}
		System.out.println();
		System.out.println("Z0: "+Z0);
		System.out.println("Z1: "+Z1);
		System.out.println("duality gap: "+dualityGap);
	}
	
	public void process(String inputFile){
		int iteration=1;
		long startTime, endTime, processTime;		
		startTime = System.currentTimeMillis();	//check start time
		
		readInputData(inputFile);
		do{
			initiate();
			calculate();
			checkDualityGap();
			System.out.println("iteration: "+iteration);
			printResult();
			x0=x1;
			iteration++;
		}while(error<=dualityGap);
		
		/*** check processing time ***/
		endTime = System.currentTimeMillis();	//check end time
		processTime = (endTime - startTime)/1000;
		System.out.println();
		System.out.print("process time: "+processTime+" sec");
		System.out.println(" ("+(endTime - startTime)+" ms)");
		System.out.println();
	}

	public void process(){
		int iteration=1;
		long startTime, endTime, processTime;		
		startTime = System.currentTimeMillis();	//check start time
		
		setInputData();
		do{
			initiate();
			calculate();
			checkDualityGap();
			System.out.println("iteration: "+iteration);
			printResult();
			x0=x1;			
			iteration++;
		}while(error<=dualityGap);

		
		/*** check processing time ***/
		endTime = System.currentTimeMillis();	//check end time
		processTime = (endTime - startTime)/1000;
		System.out.println();
		System.out.print("process time: "+processTime+" sec");
		System.out.println(" ("+(endTime - startTime)+" ms)");
		System.out.println();
	}
	
	public static void main(String[] args) {
		
		InteriorPointMethod ipm = new InteriorPointMethod();
		ipm.process("c:/optimization/input.txt");

	}

}
