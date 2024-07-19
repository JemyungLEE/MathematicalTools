package optimization;

import Jama.Matrix;

public class QuadraticProgramming {

	public Matrix quadraticProgramming(){
		double[][] xx={{1},{1}};
		Matrix x = new Matrix(xx);
		
		double[][] a={{2,1},{1,-4}};
		Matrix A = new Matrix(a);
		double[][] b={{6},{0}};
		Matrix B = new Matrix(b);
		double[][] c={{-4},{0}};
		Matrix C = new Matrix(c);
		double[][] d={{2,-2},{-2,4}};
		Matrix D = new Matrix(d);
		
		
		//Matrix A = derivativeInequality(x);
		//Matrix B = inequalityConstraint(x);
		//Matrix C = derivativeObjective(x);
		//Matrix D = hessian;		
		int xRow = x.getRowDimension();
		int bRow = B.getRowDimension();
		int simplexRow = bRow+ xRow + 1;
		int simplexColumn = 2 * xRow + 2 * bRow + 1;  
		double[][] simplex = new double[simplexRow][simplexColumn];
			//row: constrains, Q
			//column: s_i, lamda_i, Y_i, z_i, b_i
			
		int[] basis = new int[simplexRow-1];
		int entering, leaving;	//entering next basis
		int chk; //0:contains no artificial variables
		int kt;	//for check K-T condition
		double tmp, BperA, tmpBperA, tmpValue;
		Matrix result = new Matrix(bRow+ xRow,1);

		
		int i,j;
		
		//for(i=0 ; i<B.getRowDimension() ; i++)
		//	if( B.get(i, 0) > 0 ) B.set(i, 0, 0.9 * B.get(i,0));
		
		//initiate simplex tableau
		for(i=0 ; i<simplexRow ; i++)
			for(j=0 ; j<simplexColumn ; j++) simplex[i][j] = 0.0;
		for(i=0 ; i<bRow ; i++){
			for(j=0 ; j<xRow ; j++)	simplex[i][j] = A.get(i, j);
			simplex[i][xRow+bRow+i] = 1.0;
			simplex[i][simplexColumn-1] = B.get(i, 0);			
		}
		for(i=0 ; i<xRow ; i++){
			for(j=0 ; j<xRow ; j++)	simplex[i+bRow][j] = D.get(i, j);
			for(j=0 ; j<bRow ; j++) simplex[i+bRow][j+xRow] = A.get(j,i);
			simplex[i+bRow][xRow+2*bRow+i] = 1.0;
			simplex[i+bRow][simplexColumn-1] = -1.0*C.get(i, 0);
		}
				
		for(i=0 ; i<xRow ; i++){
			for(j=0 ; j<xRow+bRow ; j++) 
				simplex[simplexRow-1][j] -= simplex[i+bRow][j];
			simplex[simplexRow-1][simplexColumn-1] 
			                      -= simplex[i+bRow][simplexColumn-1];
		}
		
		//set initial basic variable list
		for(i=0 ; i<bRow ; i++) basis[i] = xRow+bRow+i;
		for(i=bRow ; i<bRow+xRow ; i++) basis[i] = xRow+bRow+i;
		
		//quadratic programming process
		do{			
			//select entering basis
			tmp = 0.0;
			entering = -1;
			for(i=0 ; i<simplexColumn-1 ; i++){
				if(i>=xRow && i<xRow+bRow){
					for(kt=0,j=0 ; j<xRow+bRow ; j++){
						if(basis[j]>=xRow+bRow && basis[j]<xRow+2*bRow 
								&& basis[j]==i+bRow) kt++;
						else if(basis[j]>=xRow && basis[j]<xRow+bRow 
								&& basis[j]==i-bRow) kt++;
					}
					if(kt==0){
						if(simplex[simplexRow-1][i]<tmp){
							tmp = simplex[simplexRow-1][i];
							entering = i;
						}	
					}					
				}
				else{
					if(simplex[simplexRow-1][i]<tmp){
						tmp = simplex[simplexRow-1][i];
						entering = i;
					}
				}								
			}
			
			System.out.println("entering:"+entering);
			
			//select leaving basis
			leaving = -1;
			BperA = 10000;	//set some large value			
			for(i=0 ; i<simplexRow-1 ; i++){
				if(simplex[i][entering]>0){
					tmpBperA = simplex[i][simplexColumn-1]/simplex[i][entering];
					if(tmpBperA<BperA){
						BperA = tmpBperA;
						leaving = i;
					}
				}
			}
				
			System.out.println("leaving:"+leaving);
			
			//update simplex tableau
			basis[leaving] = entering;	//arrange basis list
			tmpValue = simplex[leaving][entering];
			for(i=0 ; i<simplexColumn ; i++) simplex[leaving][i]/=tmpValue;			
			for(i=0 ; i<simplexRow ; i++){
				if(i!=leaving){
					tmpValue = simplex[i][entering];
					for(j=0 ; j<simplexColumn ; j++)
						simplex[i][j] -= tmpValue*simplex[leaving][j];					
				}
			}
			
			
			System.out.print("basis;");
			for(i=0 ; i<bRow+xRow ; i++) System.out.print("\t"+basis[i]);
			System.out.println();
			for(i=0 ; i<simplexRow ; i++){
				for(j=0 ; j<simplexColumn ; j++){
					System.out.print(simplex[i][j]+"\t");
				}
				System.out.println();
			}			
			
			
			//check whether basis contains artificial variables z_i 
			chk=0;
			for(i=bRow ; i<bRow+xRow ; i++)
				if(basis[i]>=xRow+bRow+bRow) chk++;
			
		}while(chk>0);
		
		//set the result values
		//s_i, lamda_i
		for(i=0 ; i<bRow+xRow ; i++){
			for(j=0 ; j<xRow ; j++)
				if(basis[i]==j) result.set(j, 0, simplex[i][simplexColumn-1]);
			for(j=0 ; j<bRow ; j++){
				if(basis[i]==xRow+j) 
					result.set(xRow+j,0,simplex[i][simplexColumn-1]);
				if(basis[i]==xRow+bRow+j) 
					result.set(xRow+j,0,0.0);
			}
		}
		/*
		for(i=0 ; i<xRow ; i++) 
			System.out.println("x"+(i+1)+":"+result.get(i,0));
		for(i=0 ; i<bRow ; i++)
			System.out.println("lamda"+(i+1)+":"+result.get(i+xRow,0));
		*/
		return result;
	}

	public static void main(String[] args) {
		QuadraticProgramming qp = new QuadraticProgramming();
		qp.quadraticProgramming();
		System.out.println("process complete");
	}
}
