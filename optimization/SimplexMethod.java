package optimization;

import java.io.File;
import java.io.IOException;
import java.util.Scanner;

public class SimplexMethod {

	/**
	 *  Subject: Numerical Methods for Engineers
	 *  Name: Jemyung Lee
	 *  ID: 2008-30334
	 *  
	 *  Homework Number: 3
	 *  Description: Simplex Method (Constrained optimization)
	 */
	
	public SimplexMethod(String source){
		
		algorithm(source);
	}
	
	public void algorithm(String source){
		
		int i,j,k;
		int cnt;	
		int chk;
		int leaving=0;
		int entering=0;
		double co;		//coefficient of entering variable
		double tmp=0.0;
		int surplus=0;
		int solution=0;
		double[] intercept;

		int[] nb=null;
		int[] order=null;
		int eq=0; 
		int nv=0; 		//number of variables
		int nsv=0; 		//number of slack variables
		double[][] table=null;
		
		/*** read data from file ***/
		try{
			File file = new File(source);
			Scanner scan = new Scanner(file);
			
			//read equation number
			eq = scan.nextInt();
			//read variable number
			nv = scan.nextInt();
			//read slack variable number
			nsv = scan.nextInt();
			//set surplus number
			surplus = nv+nsv-eq+1;
			//set solution
			solution = nv+nsv+2;
			//read start point
			nb = new int[surplus];
			for(i=0 ; i<surplus ; i++) nb[i] = scan.nextInt();
			
			//read coefficients
			table = new double[eq][solution];
			for(i=0;i<eq;i++){
				//first row is objective function 
				for(j=0 ; j<solution ; j++){					
					table[i][j] = scan.nextDouble();
				}
			}
			
			scan.close();			
		} catch(IOException e) {}

		
		/*** calculation part ***/
				
		order = new int[eq-1];
		intercept = new double[eq];
		
		//set initial entering variable
		cnt=0;
		for(i=0; i<solution ; i++){
			if(table[0][i]<0){
				entering=i;
				cnt=1;
				break;
			}
		}
		
		//set initial variable order
		chk=0;
		for(i=0,k=0 ; i<nv+nsv ; i++){
			for(j=0 ; j<surplus ; j++) if((nb[j]-1)==i) chk = 1;
			if(chk == 0){
				order[k] = i;
				k++;
			}
			chk=0;
		}
		
		//iteration
		while(cnt>0){
			
			//calculation of intercept
			for(i=1 ; i<eq ; i++){
				if(table[i][entering]!=0) 
					intercept[i] = table[i][nv+nsv+1]/table[i][entering];
				else if(table[i][entering]==0)
					intercept[i] = 10000;
			}
						
			//define leaving variable
			chk=0;
			for(i=1 ; i<eq ; i++){
				if(chk==0 && intercept[i]>=0){
					tmp = intercept[i];
					leaving = i;
					chk = 1;
				}
				else if(intercept[i]<tmp && intercept[i]>=0){
					tmp = intercept[i];
					leaving = i;
				}
			}
			
			//set variable order
			order[leaving-1]=entering-2;
			
			//Gauss-Jordan Method
			co = table[leaving][entering];
			for(i=0 ; i<nv+nsv+2 ; i++)	
				table[leaving][i] = table[leaving][i]/co;
			for(i=0 ; i<eq ; i++){
				if(i != leaving){				
					tmp = table[i][entering];
					for(j=0 ; j<nv+nsv+2; j++) 
						table[i][j] -= table[leaving][j]*tmp;
				}
			}
				
			//set next entering variable
			cnt=0;
			for(i=0; i<solution ; i++){
				if(table[0][i]<0){
					entering=i;
					cnt=1;
					break;
				}
			}
		
			/*** print the tableau ***/
			
			System.out.print("Basic\tZ");
			for(i=0 ; i<nv ; i++) System.out.print("\tX"+(i+1));
			for(i=0 ; i<nsv ; i++) System.out.print("\tS"+(i+1));
			System.out.println("\tSolution");
			for(i=0 ; i<eq ; i++){
				if(i==0) System.out.print("Z");
				else if(order[i-1]<nv) System.out.print("X"+(order[i-1]+2));
				else System.out.print("S"+(order[i-1]-nv+1));
				for(j=0 ; j<nv+nsv+2 ; j++) 
					System.out.printf("\t%6.4f",table[i][j]);
				System.out.println();
			}	
			System.out.println();
		}

		
		/*** print the result ***/	
		System.out.println("RESULT");
		for(i=0 ; i<nv ; i++){
			for(j=1 ; j<eq ; j++){
				if(order[j-1]+1==i) 
					System.out.println("X"+(i+1)+" = "+table[j][nv+nsv+1]);
			}
		}
		System.out.println("Z  = "+table[0][nv+nsv+1]);
		
	}
	
	public static void main(String[] args) {
		String source; 
		
		source = new String("c:/input.txt");
		new SimplexMethod(source);
	}

}
