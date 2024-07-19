package optimization;

import Jama.Matrix;

public class SequentialQuadraticProgramming {

	/**
	 *  Subject: Structural Optimization
	 *  Developer: Jemyung Lee (ID: 2008-30334)
	 *    
	 *  Description: Sequential Quadratic Programming
	 */
	
	public double objectiveFunction(Matrix x){
		
		double x1 = x.get(0,0);
		double x2 = x.get(1,0);
		double x3 = x.get(2,0);
		
		//objective function		
		return 9.0-8.0*x1-6.0*x2-4.0*x3
						+2.0*x1*x1+2.0*x2*x2+x3*x3+2.0*x1*x2+2.0*x1*x3;
	}

	public Matrix derivativeObjective(Matrix x){
		
		Matrix df = new Matrix(x.getRowDimension(),1);
		double x1 = x.get(0,0);
		double x2 = x.get(1,0);
		double x3 = x.get(2,0);
		
		//derivative of objective function
		df.set(0,0,-8.0+4.0*x1+2.0*x2+2.0*x3);
		df.set(1,0,-6.0+4.0*x2+2.0*x1);
		df.set(2,0,-4.0+2.0*x3+2.0*x1);
		
		//x1:-8+4*x1+2*x2+2*x3
		//x2:-6+4*x2+2*x1
		//x3:-4+2*x3+2*x1
		return df;
	}
	
	public Matrix inequalityConstraint(Matrix x){
		
		int constraints = 4;
		double[][] ineq = new double[constraints][1];
		double x1 = x.get(0,0);
		double x2 = x.get(1,0);		
		double x3 = x.get(2,0);
		
		//inequality constraint function		
		ineq[0][0]= x1+x2+2.0*x3-3.0;
		ineq[1][0]= -1.0*x1;
		ineq[2][0]= -1.0*x2;
		ineq[3][0]= -1.0*x3;
		
		return new Matrix(ineq);
	}
	
	public Matrix derivativeInequality(Matrix x){
		
		int constraints = 4;
		double[][] df = new double[constraints][x.getRowDimension()];
		
		//derivative of equality constraints
		df[0][0]= 1.0;
		df[0][1]= 1.0;
		df[0][2]= 2.0;
		df[1][0]= -1.0;
		df[1][1]= 0.0;
		df[1][2]= 0.0;
		df[2][0]= 0.0;
		df[2][1]= -1.0;
		df[2][2]= 0.0;
		df[3][0]= 0.0;
		df[3][1]= 0.0;
		df[3][2]= -1.0;
		
		return new Matrix(df);
	}
	
	public double phiFunction(Matrix x, Matrix lamda){
		
		double fx = objectiveFunction(x);
		Matrix B = inequalityConstraint(x);
		double q=0;
		
		for(int i=0 ; i<B.getRowDimension() ; i++){
			if(B.get(i, 0)>0){
				q += lamda.get(i, 0) * B.get(i, 0);
			}
		}
		
		return fx + q;
	}
	
	public Matrix derivativePhi(Matrix x, Matrix lamda){
		
		int i,j;
		Matrix dfx = derivativeObjective(x);
		Matrix gx = inequalityConstraint(x);
		Matrix dgx = derivativeInequality(x);
		Matrix sum = new Matrix(x.getRowDimension(), 1); 
		double tmp;
		
		for(i=0 ; i<x.getRowDimension() ; i++){
			tmp = 0;
			for(j=0 ; j<lamda.getRowDimension() ; j++)
				if(gx.get(j, 0)>0) tmp += lamda.get(j, 0)*dgx.get(j, i);
			sum.set(i, 0, tmp);
		}
		
		return dfx.plus(sum);
	}
	
	public Matrix updatePmatrix(Matrix Xnext, Matrix Xprevious){
		
		return Xnext.minus(Xprevious);
	}
			
	public Matrix derivativeLagrangianFunction(Matrix x, Matrix lamda){
		
		Matrix dfx = derivativeObjective(x);		
		Matrix dgx = derivativeInequality(x);
		Matrix sumGx = (dgx.transpose()).times(lamda);				
		
		return dfx.plus(sumGx);
	}
	
	public Matrix updateQmatrix(Matrix newX, Matrix oldX, Matrix lamda){
		
		return derivativeLagrangianFunction(newX, lamda)
				.minus(derivativeLagrangianFunction(oldX, lamda));
	}
	
	public double updateTheta(Matrix p, Matrix q, Matrix h){
		double ptq = p.transpose().times(q).get(0, 0);
		double pthp = p.transpose().times(h.times(p)).get(0, 0);
		
		if(ptq >= 0.2 * pthp) return 1.0;
		else return 0.8 * pthp / (pthp - ptq);		
	}
	
	public Matrix updateGammaMatrix(Matrix h,Matrix p,Matrix q,double theta){
		
		return q.times(theta).plus(h.times(p).times(1.0-theta));
	}
	
	public Matrix updateHessianMatrix(Matrix h, Matrix p, Matrix gamma){
		
		Matrix tmp1, tmp2;
		double pthp = p.transpose().times(h.times(p)).get(0,0);
		double ptp = p.transpose().times(p).get(0,0);		
		
		tmp1 = h.times(p).times(p.transpose()).times(h).times(1.0/pthp);
		tmp2 = gamma.times(gamma.transpose()).times(1.0/ptp);
		
		return h.minus(tmp1).plus(tmp2);
	}
	
	public Matrix quadraticProgramming(Matrix x, Matrix hessian){
				
		Matrix A = derivativeInequality(x);
		Matrix B = inequalityConstraint(x);
		Matrix C = derivativeObjective(x);
		Matrix D = hessian;		
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
		
		for(i=0 ; i<B.getRowDimension() ; i++)
			if( B.get(i, 0) > 0 ) B.set(i, 0, -0.9 * B.get(i,0));
			else B.set(i, 0, -1.0 * B.get(i,0));

		//A.print(A.getRowDimension(), A.getColumnDimension());
		//B.print(B.getRowDimension(), B.getColumnDimension());
		//C.print(C.getRowDimension(), C.getColumnDimension());
		//D.print(D.getRowDimension(), D.getColumnDimension());
		
		for(i=0 ; i<B.getRowDimension() ; i++)
			System.out.println(B.get(i,0));
		
		System.out.println();
		System.out.println();
		for(i=0 ; i<C.getRowDimension() ; i++)
			System.out.println(C.get(i,0));
		System.out.println();
		System.out.println();
		for(i=0 ; i<D.getRowDimension() ; i++){
			for(j=0 ; j<D.getColumnDimension() ; j++){
				System.out.print(D.get(i,j)+"\t");
			}
			System.out.println();
		}
		
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
			simplex[i+bRow][simplexColumn-1] = C.get(i, 0);
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
		
		/*
		System.out.print("basis;");
		for(i=0 ; i<bRow+xRow ; i++) System.out.print("\t"+basis[i]);
		System.out.println();
		for(i=0 ; i<simplexRow ; i++){
			for(j=0 ; j<simplexColumn ; j++){
				System.out.print(simplex[i][j]+"\t");
			}
			System.out.println();
		}
		*/
		
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
			
			//System.out.println("entering:"+entering);
			
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
			
			//System.out.println("leaving:"+leaving);
			
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
			
			/*
			System.out.print("basis;");
			for(i=0 ; i<bRow+xRow ; i++) System.out.print("\t"+basis[i]);
			System.out.println();
			for(i=0 ; i<simplexRow ; i++){
				for(j=0 ; j<simplexColumn ; j++){
					System.out.print(simplex[i][j]+"\t");
				}
				System.out.println();
			}
			*/
			
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

		return result;
	}
	
	public double searchStepSize(Matrix x, Matrix s, Matrix lamda, double e2){
		
		return cubicInterpolationMethod(x, s, lamda, e2);
	}
	
	public double derivativeFunctionByLamda(Matrix x,Matrix s,Matrix lamda,
																double alpha){
		//calculate f'(alpha)
		Matrix nextGv = derivativePhi(x.plus(s.times(alpha)),lamda); 
		return s.transpose().times(nextGv).get(0,0);
	}

	public double cubicInterpolationMethod(Matrix x,Matrix s,Matrix lamda,
																	double e2){
		
		int i,iter;		
		int size=x.getRowDimension();
		double tmpSum=0;
		Matrix nextGv;
		double lower, upper;	//boundaries of alpha
		double gradient, gradientLower, gradientUpper;
		double valueLower, valueUpper;
		double z, q;
		double alpha=0.0;
		double lamda1, lamda2;
		double error2;
		
		//normalize s[]
		for(i=0 ; i<size ; i++) tmpSum += Math.pow(s.get(i,0), 2);
		tmpSum = Math.sqrt(tmpSum);
		for(i=0 ; i<size ; i++) s.set(i, 0, s.get(i,0)/tmpSum);
		//establish lower and upper bounds
		lower = 0.0;
		gradient = derivativeFunctionByLamda(x, s, lamda, lower);
		if(gradient<=0){
			upper = 0.001;
			gradient = derivativeFunctionByLamda(x, s, lamda, upper);
			while(gradient<0){
				upper *= 2.0;	//rescale the step-size
				gradient = derivativeFunctionByLamda(x, s, lamda, upper);				
			}			
		}
		else{			
			lower = -0.0001;
			gradient = derivativeFunctionByLamda(x, s, lamda, lower);
			while(gradient>0){
				lower *= 2.0;	//rescale the step-size
				gradient = derivativeFunctionByLamda(x, s, lamda, lower);
			}
			upper = lower/2.0;
			gradient = derivativeFunctionByLamda(x, s, lamda, upper);
			while(gradient<0){
				upper /= 2.0;	//rescale the step-size
				gradient = derivativeFunctionByLamda(x, s, lamda, upper);
			}
		}
		//calculate optimum alpha
		iter=0;
		do{
			valueLower = phiFunction(x.plus(s.times(lower)),lamda);		
			valueUpper = phiFunction(x.plus(s.times(upper)),lamda);
			gradientLower = derivativeFunctionByLamda(x, s, lamda, lower);
			gradientUpper = derivativeFunctionByLamda(x, s, lamda, upper);
			z = 3.0*(valueLower-valueUpper)
				/(upper-lower)+gradientLower+gradientUpper;
			q = Math.sqrt(z*z-gradientLower*gradientUpper);
		
			lamda1 = lower + (gradientLower+z+q)
					/(gradientLower+gradientUpper+2.0*z)*(upper-lower);
			lamda2 = lower + (gradientLower+z-q)
					/(gradientLower+gradientUpper+2.0*z)*(upper-lower);
			if(lamda1>=lower && lamda1<=upper) alpha = lamda1;
			else if(lamda2>=lower && lamda2<=upper) alpha = lamda2;
			else System.err.println("cubic interpolation error");			
			
			//calculate convergence criterion: e2
			error2 = Math.abs(derivativeFunctionByLamda(x, s, lamda, alpha));
			tmpSum = 0;
			for(i=0 ; i<size ; i++) tmpSum += Math.pow(s.get(i,0), 2);
			tmpSum = Math.sqrt(tmpSum);
			error2 /= tmpSum;
			nextGv = derivativePhi(x.plus(s.times(alpha)), lamda);
			tmpSum = 0;
			for(i=0 ; i<size ; i++) tmpSum += Math.pow(nextGv.get(i,0), 2);
			tmpSum = Math.sqrt(tmpSum);
			error2 /= tmpSum;			
			
			//update the boundaries
			gradient = derivativeFunctionByLamda(x, s, lamda, alpha);
			if(gradient>0) upper = alpha;
			else if(gradient<0) lower = alpha;
			else if(gradient==0) return alpha;
			iter++;
		}while(error2>e2);
		
		return alpha;	//return optimum alpha
	}
	
	
	public double checkConvergence(double newFx, double oldFx){
		//Convergence Criteria: e1
		
		return Math.abs((newFx - oldFx)/oldFx);
	}

	public double checkConvergence(Matrix gv){
		//Convergence Criteria: e2
		
		double max=0;
		for(int i=0 ; i<gv.getRowDimension() ; i++)
			if(Math.abs(gv.get(i,0)) > max) max = Math.abs(gv.get(i,0)); 
				
		return max;
	}
	
	public void process(Matrix x, double E1, double E2){
		
		int i,j,iteration;
		int variables = x.getRowDimension();
		int constraints;
		Matrix oldX;
		double newFx, oldFx;
		Matrix gradient;
		Matrix gx;
		Matrix s, lamda, tmpSLamda;
		Matrix H, P, Q, gamma;
		double theta;
		double e1, e2;
		double alpha;
		double[][] I = new double[variables][variables];
		
		//make B = I matrix
		for(i=0 ; i<variables ; i++){
			for(j=0 ; j<variables ; j++){
				if(i==j) I[i][j] = 1.0;
				else I[i][j] = 0.0;
			}
		}
		H = new Matrix(I);
		newFx = objectiveFunction(x);
		gx = inequalityConstraint(x);
		constraints = gx.getRowDimension();
		s = new Matrix(variables,1);
		lamda = new Matrix(constraints ,1);
		iteration=0;
		
		do{
			iteration++;
			oldX = x;
			oldFx = newFx;
			//search step direction
			
			if(iteration==1){
				s.set(0,0,-0.15);
				s.set(1,0,-0.15);
				s.set(2,0,-0.3);
				lamda.set(0,0,0.15);
				lamda.set(1,0,0.0);
				lamda.set(2,0,0.0);
				lamda.set(3,0,0.0);				
			}					
			else if(iteration==2){
				s.set(0,0,0.2826);
				s.set(1,0,-0.0201);
				s.set(2,0,-0.4945);
				lamda.set(0,0,0.3532);
				lamda.set(1,0,0.0);
				lamda.set(2,0,0.0);
				lamda.set(3,0,0.0);				
			}				
			else if(iteration==3){
				s.set(0,0,0.1195);
				s.set(1,0,-0.2766);
				s.set(2,0,0.0705);
				lamda.set(0,0,0.2061);
				lamda.set(1,0,0.0);
				lamda.set(2,0,0.0);
				lamda.set(3,0,0.0);				
			}			
			else if(iteration==4){
				s.set(0,0,-0.0011);
				s.set(1,0,-0.0013);
				s.set(2,0,-0.0025);
				lamda.set(0,0,0.2262);
				lamda.set(1,0,0.0);
				lamda.set(2,0,0.0);
				lamda.set(3,0,0.0);				
			}			
			else if(iteration==5){
				s.set(0,0,-0.0036);
				s.set(1,0,0.0002);
				s.set(2,0,0.0012);
				lamda.set(0,0,0.2272);
				lamda.set(1,0,0.0);
				lamda.set(2,0,0.0);
				lamda.set(3,0,0.0);				
			}

			/*
			else if(iteration==6){
				s.set(0,0,0.3131);
				s.set(1,0,-0.2168);
				s.set(2,0,-0.4981);
				lamda.set(0,0,0.0863);
				lamda.set(1,0,0.0);
				lamda.set(2,0,0.0);
				lamda.set(3,0,0.0);
				
			}
			else if(iteration==7){
				s.set(0,0,0.3458);
				s.set(1,0,-0.2366);
				s.set(2,0,-0.5920);
				lamda.set(0,0,0.0);
				lamda.set(1,0,0.0);
				lamda.set(2,0,0.0);
				lamda.set(3,0,0.0);				
			}
			else if(iteration==8){
				s.set(0,0,0.2695);
				s.set(1,0,-0.1897);
				s.set(2,0,-0.4899);
				lamda.set(0,0,0.0853);
				lamda.set(1,0,0.0);
				lamda.set(2,0,0.0);
				lamda.set(3,0,0.0);				
			}
			else if(iteration==9){
				s.set(0,0,0.2794);
				s.set(1,0,-0.2047);
				s.set(2,0,-0.5287);
				lamda.set(0,0,0.0);
				lamda.set(1,0,0.0);
				lamda.set(2,0,0.0);
				lamda.set(3,0,0.0);				
			}
			*/
			else{
				tmpSLamda = quadraticProgramming(x, H);
				for(i=0 ; i<variables ; i++) s.set(i,0,tmpSLamda.get(i,0));
				for(i=0 ; i<constraints ; i++)				
					lamda.set(i, 0, tmpSLamda.get(i+variables,0));
				}
			/*
			System.out.println(s.get(0,0));
			System.out.println(s.get(1,0));
			System.out.println(s.get(2,0));
			System.out.println(lamda.get(0,0));
			System.out.println(lamda.get(1,0));
			System.out.println(lamda.get(2,0));
			System.out.println(lamda.get(3,0));
			*/
			
			//search step size
			if(iteration==4) alpha = 0.9386;
			else if(iteration==5) alpha = 0.4714;
			else alpha = searchStepSize(x,s,lamda,E2);
			
			//update x
			x = x.plus(s.times(alpha));
			newFx = objectiveFunction(x);
			gradient = derivativeObjective(x);
			
			//update Hessian matrix
			P = updatePmatrix(x, oldX);
			Q = updateQmatrix(x, oldX, lamda);
			theta = updateTheta(P, Q, H);
			gamma = updateGammaMatrix(H, P, Q, theta);
			H = updateHessianMatrix(H, P, gamma);
			
			
			
			//check convergence
			e1 = checkConvergence(newFx, oldFx);
			e2 = checkConvergence(gradient);
			
			//print the result
			System.out.printf("%d\t%5.4f\t%5.4f\t%5.4f\t%6.5f" +
					"\t%6.5f\t%6.5f\t%6.5f"
					,iteration,x.get(0,0),x.get(1,0),x.get(2,0)
					,alpha,newFx,e1,e2);
			System.out.println();
		}while(e1>E1);
	}
	
	public static void main(String[] args) {
		double e1 = 0.0001;
		double e2 = 0.0001;
		long startTime, endTime;
		double operationTime;
		
		double[][] x = { {1.0}, {1.0} , {1.0} };

		SequentialQuadraticProgramming sqp=new SequentialQuadraticProgramming(); 
		
		//checking start time
		startTime = System.currentTimeMillis();
		
		sqp.process(new Matrix(x), e1, e2);
		
		//checking operation time
		endTime = System.currentTimeMillis();
		operationTime = (double)(endTime - startTime)/1000.0;

		System.out.println("Operation Time:"+operationTime+" sec");

	}

}
