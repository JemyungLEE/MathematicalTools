package optimization;

import Jama.Matrix;

public class BFGStest {

	/**
	 *  Subject: Structural Optimization
	 *  Developer: Jemyung Lee (ID: 2008-30334)
	 *    
	 *  Description: Broyden-Fletcher-Goldfarb-Shanno Method
	 */

	public double objectiveFunction(Matrix x){
		
		double r=1.0;
		
		//objective function
		
		return Math.pow(x.get(0,0),3)-6.0*Math.pow(x.get(0,0),2)+11.0*x.get(0,0)+x.get(2,0)
		+r*(1.0/(Math.pow(x.get(0,0),2)+Math.pow(x.get(1,0),2)-Math.pow(x.get(2,0),2))
			+1.0/(4.0-Math.pow(x.get(0,0),2)-Math.pow(x.get(1,0),2)-Math.pow(x.get(2,0),2))
			+1.0/(x.get(2,0)-5.0));
		
		/*
		return 1.0/3.0*Math.pow(x.get(0,0)+1.0,3)+x.get(1,0)
						-r*(1/(1.0-x.get(0,0))-1.0/x.get(1,0));
		 */
	}

	public Matrix gradientVector(Matrix x){
		
		double r=1.0;
		
		Matrix gv = new Matrix(x.getRowDimension(),1);	//gv[i-1]=df/dx[i]
		
		//partial derivative of objective function

		gv.set(0,0,3*Math.pow(x.get(0,0),2)-12.0*x.get(0,0)+11.0
				-r*(-2.0*x.get(0,0)/(Math.pow(Math.pow(x.get(0,0),2)+Math.pow(x.get(1,0),2)-Math.pow(x.get(2,0),2), 2))
						+2.0*x.get(0,0)/Math.pow((4.0-Math.pow(x.get(0,0),2)-Math.pow(x.get(1,0),2)-Math.pow(x.get(2,0),2)),2)
						));
		gv.set(1,0,-1.0*r*(-2.0*x.get(1,0)/(Math.pow(Math.pow(x.get(0,0),2)+Math.pow(x.get(1,0),2)-Math.pow(x.get(2,0),2), 2))
				+2.0*x.get(1,0)/Math.pow((4.0-Math.pow(x.get(0,0),2)-Math.pow(x.get(1,0),2)-Math.pow(x.get(2,0),2)),2)
				));
		
		gv.set(2,0,1.0-r*(-2.0*x.get(2,0)/(Math.pow(Math.pow(x.get(0,0),2)+Math.pow(x.get(1,0),2)-Math.pow(x.get(2,0),2), 2))
				+2.0*x.get(2,0)/Math.pow((4.0-Math.pow(x.get(0,0),2)-Math.pow(x.get(1,0),2)-Math.pow(x.get(2,0),2)),2)
				-1.0/Math.pow(x.get(2,0)-5.0, 2)));
		//3*x1^2-12*x1+11-r*(-2/(x1^2+x2^2-x3^2)^2*x1+2/(4-x1^2-x2^2-x3^2)^2*x1)
		
		//-r*(-2/(x1^2+x2^2-x3^2)^2*x2+2/(4-x1^2-x2^2-x3^2)^2*x2)
		
		//1-r*(2/(x1^2+x2^2-x3^2)^2*x3+2/(4-x1^2-x2^2-x3^2)^2*x3-1/(x3-5)^2)
		
		/*
		gv.set(0,0,Math.pow((x.get(0,0)+1.0),2)-r/Math.pow((1.0-x.get(0,0)),2));

		gv.set(1,0,1.0-r/Math.pow(x.get(1,0), 2));
		 */
		return gv;
	}
	
	public Matrix searchDirection(Matrix gv, Matrix B){
				
		return B.times(gv).times(-1);
	}
	
	public Matrix makeGmatrix(Matrix oldGv, Matrix newGv){
		
		return newGv.minus(oldGv);
	}
	
	public Matrix makeDmatrix(double lamda, Matrix S){
		
		return S.times(lamda);
	}
	
	public Matrix updateBmatrix(Matrix B, Matrix D, Matrix g){
		Matrix tD = D.transpose();
		Matrix tg = g.transpose();
		Matrix tmp1, tmp2, tmp3, tmp4;
		double temp;
		
		temp = tD.times(g).get(0,0);
		tmp1 = D.times(tD).times(1.0/temp);
		tmp2 = tg.times(B).times(g).times(1.0/temp);
		tmp3 = B.times(g.times(tD)).times(1.0/temp);
		tmp4 = D.times(tg.times(B)).times(1.0/temp);
		
		return B.plus(tmp1.times(1+tmp2.get(0,0))).minus(tmp3).minus(tmp4);
	}
	
	public double searchStepSize(Matrix x, Matrix s, double e2){
		
		return cubicInterpolationMethod(x, s, e2);
	}
	
	public double derivativeFunctionByLamda(Matrix x,Matrix s, double lamda){

		//calculate f'(lamda)
		Matrix nextGv = gradientVector(x.plus(s.times(lamda))); 
		return s.transpose().times(nextGv).get(0,0);
	}
	
	public double cubicInterpolationMethod(Matrix x,Matrix s,double e2){
	
		int i,iter;		
		int size=x.getRowDimension();
		double tmpSum=0;
		Matrix nextGv;
		double lower, upper;	//boundaries of lamda
		double gradient, gradientLower, gradientUpper;
		double valueLower, valueUpper;
		double z, q;
		double lamda=0.0;
		double lamda1, lamda2;
		double error2;
		
		//normalize s[]
		for(i=0 ; i<size ; i++) tmpSum += Math.pow(s.get(i,0), 2);
		tmpSum = Math.sqrt(tmpSum);
		for(i=0 ; i<size ; i++) s.set(i, 0, s.get(i,0)/tmpSum);
		//establish lower and upper bounds
		lower = 0.0;
		gradient = derivativeFunctionByLamda(x, s, lower);
		if(gradient<=0){
			upper = 0.1;
			gradient = derivativeFunctionByLamda(x, s, upper);
			while(gradient<0){
				upper *= 2.0;	//rescale the step-size
				gradient = derivativeFunctionByLamda(x, s, upper);				
			}			
		}
		else{			
			lower = -0.0001;
			gradient = derivativeFunctionByLamda(x, s, lower);
			while(gradient>0){
				lower *= 2.0;	//rescale the step-size
				gradient = derivativeFunctionByLamda(x, s, lower);
			}
			upper = lower/2.0;
			gradient = derivativeFunctionByLamda(x, s, upper);
			while(gradient<0){
				upper /= 2.0;	//rescale the step-size
				gradient = derivativeFunctionByLamda(x, s, upper);
			}
		}
		System.out.println("upper"+upper+"\tlower:"+lower);
		//calculate optimum lamda
		iter=0;
		do{
			valueLower = objectiveFunction(x.plus(s.times(lower)));		
			valueUpper = objectiveFunction(x.plus(s.times(upper)));
			gradientLower = derivativeFunctionByLamda(x, s, lower);
			gradientUpper = derivativeFunctionByLamda(x, s, upper);
			z = 3.0*(valueLower-valueUpper)
				/(upper-lower)+gradientLower+gradientUpper;
			q = Math.sqrt(z*z-gradientLower*gradientUpper);
		
			lamda1 = lower + (gradientLower+z+q)
					/(gradientLower+gradientUpper+2.0*z)*(upper-lower);
			lamda2 = lower + (gradientLower+z-q)
					/(gradientLower+gradientUpper+2.0*z)*(upper-lower);
			if(lamda1>=lower && lamda1<=upper) lamda = lamda1;
			else if(lamda2>=lower && lamda2<=upper) lamda = lamda2;
		//	else System.err.println("cubic interpolation error");			
			
			//calculate convergence criterion: e2
			error2 = Math.abs(derivativeFunctionByLamda(x, s, lamda));
			tmpSum = 0;
			for(i=0 ; i<size ; i++) tmpSum += Math.pow(s.get(i,0), 2);
			tmpSum = Math.sqrt(tmpSum);
			error2 /= tmpSum;
			nextGv = gradientVector(x.plus(s.times(lamda)));
			tmpSum = 0;
			for(i=0 ; i<size ; i++) tmpSum += Math.pow(nextGv.get(i,0), 2);
			tmpSum = Math.sqrt(tmpSum);
			error2 /= tmpSum;			
			
			//update the boundaries
			gradient = derivativeFunctionByLamda(x, s, lamda);
			if(gradient>0) upper = lamda;
			else if(gradient<0) lower = lamda;
			else if(gradient==0) return lamda;
			iter++;
		}while(error2>e2);
		//System.out.println("cubic iter:"+iter);
		return lamda;	//return optimum lamda
	}
	
	public double checkConvergence(double newFx, double oldFx){
		//Convergence Criteria: e1
		
		return Math.abs((newFx - oldFx)/oldFx);
	}

	public double checkConvergence(Matrix gv){
		//Convergence Criteria: e2
		
		double max=0;
		for(int i=0 ; i<4 ; i++)
			if(Math.abs(gv.get(i,0)) > max) max = Math.abs(gv.get(i,0)); 
				
		return max;
	}
	
	public Matrix process(Matrix x, double E1, double E2){
		
		int i,j,iter;
		int size = x.getRowDimension();
		double newFx, oldFx;
		Matrix newGv, oldGv;
		Matrix S;
		Matrix D;
		double[][] b = new double[size][size];
		Matrix B;
		Matrix g;
		double e1, e2;
		double lamda;
		
		//make b = I matrix
		for(i=0 ; i<size ; i++){
			for(j=0 ; j<size ; j++){
				if(i==j) b[i][j] = 1.0;
				else b[i][j] = 0.0;
			}
		}
		B = new Matrix(b);
		oldFx = objectiveFunction(x);
		newGv = gradientVector(x);
		S = searchDirection(newGv, B);
		lamda = searchStepSize(x, S, E2);
		x = x.plus(S.times(lamda));
		newFx = objectiveFunction(x);
		e1 = checkConvergence(newFx, oldFx);
		e2 = checkConvergence(newGv);
		iter=1;
		
		System.out.println("iter\tx[0]\tx[1]\tx[2]\tx[3]"+
							"\tlamda\tnewFx\te1\te2");
		System.out.printf("%d\t%5.4f\t%5.4f\t%5.4f\t%5.4f\t%6.5f" +
				"\t%6.5f\t%6.5f\t%6.5f"
							,iter,x.get(0,0),x.get(1,0),x.get(2,0),x.get(3,0)
							,lamda,newFx,e1,e2);
		System.out.println();
		
		for(iter=2; e2>E2; iter++){
			oldFx = newFx;
			oldGv = newGv;
			newGv = gradientVector(x);
			g = makeGmatrix(oldGv, newGv);
			D = makeDmatrix(lamda, S);
			B = updateBmatrix(B, D, g);
			S = searchDirection(newGv, B);
			lamda = searchStepSize(x, S, E2);
			x = x.plus(S.times(lamda));
			newFx = objectiveFunction(x);
			e1 = checkConvergence(newFx, oldFx);
			e2 = checkConvergence(newGv);
			System.out.printf("%d\t%5.4f\t%5.4f\t%5.4f\t%5.4f\t%6.5f" +
					"\t%6.5f\t%6.5f\t%6.5f"
					,iter,x.get(0,0),x.get(1,0),x.get(2,0),x.get(3,0)
					,lamda,newFx,e1,e2);
			System.out.println();	
		}			
		
		return x;
	}
	
	public static void main(String[] args) {

		double[] x = new double[2];
		double e1 = 0.01;
		double e2;
		long startTime, endTime;
		double operationTime;
		
		double[][] xx = { {0.5}, {0.5} , {3.0} };

		BFGStest bfgs = new BFGStest(); 
		
		//checking start time
		startTime = System.currentTimeMillis();

		System.out.println("a=0.01");
		e2 = 0.001*e1;
		bfgs.process(new Matrix(xx), e1, e2);
		
		//checking operation time
		endTime = System.currentTimeMillis();
		operationTime = (double)(endTime - startTime)/1000.0;

		System.out.println("Operation Time:"+operationTime+" sec");
	}
	
}
