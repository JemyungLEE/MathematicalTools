package optimization;

import Jama.Matrix;

public class FRmethod2 {

	/**
	 *  Subject: Structural Optimization
	 *  Developer: Jemyung Lee (ID: 2008-30334)
	 *    
	 *  Description: Fletcher-Reeves (conjugate gradient) Method
	 */

	public double objectiveFunction(Matrix x){
		
		//objective function
		return 100.0*Math.pow(x.get(1,0)-x.get(0,0)*x.get(0,0),2)
				+Math.pow(1.0-x.get(0,0)*x.get(0,0),2)+
				90.0*Math.pow(x.get(3,0)-x.get(2,0)*x.get(2,0),2)
				+Math.pow(1-x.get(2,0),2)+10.1*(Math.pow(x.get(1,0)-1.0, 2)
				+Math.pow(x.get(3,0)-1.0,2))
				+19.8*(x.get(1,0)-1.0)*(x.get(3,0)-1.0);
	}

	public Matrix gradientVector(Matrix x){
		
		Matrix gv = new Matrix(x.getRowDimension(),1);	//gv[i-1]=df/dx[i]
		
		//partial derivative of objective function
		gv.set(0,0,
				-400.0*(x.get(1,0)
				-x.get(0,0)*x.get(0,0))*x.get(0,0)-2.0+2.0*x.get(0,0));
		gv.set(1,0,
				1101.0/5.0*x.get(1,0)
				-200.0*x.get(0,0)*x.get(0,0)-40.0+99.0/5.0*x.get(3,0));
		gv.set(2,0,
				-360.0*(x.get(3,0)
				-x.get(2,0)*x.get(2,0))*x.get(2,0)-2.0+2.0*x.get(2,0));
		gv.set(3,0,
				1001.0/5.0*x.get(3,0)
				-180.0*x.get(2,0)*x.get(2,0)-40.0+99.0/5.0*x.get(1,0));
			
		return gv;
	}

	public Matrix searchDirection(Matrix gv){

		int size = gv.getRowDimension();
		Matrix s = new Matrix(size,1);		
		for(int i=0 ; i<size ; i++) s.set(i,0, -gv.get(i,0));

		return s;
	}
	
	public Matrix searchDirection(Matrix oldS,Matrix oldGv,Matrix newGv){
		
		int size = oldS.getRowDimension();
		Matrix s = new Matrix(size,1);
		for(int i=0 ; i<size ; i++)
			s.set(i,0,-1.0*newGv.get(i,0)
					+Math.pow(newGv.get(i,0)/oldGv.get(i,0), 2)*oldS.get(i,0));

		return s;
	}
	
	public double searchStepSize(Matrix x, Matrix s, double e2){
		
		return cubicInterpolationMethod(x, s, e2);
		//return goldenSectionSearch(x, s, e2);
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
			else System.err.println("cubic interpolation error");			
			
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
	
	public double goldenSectionSearch(Matrix x,Matrix s,double e2){
		
		int i,iter;
		double lamda1, lamda2;
		double lamdaL, lamdaU;
		Matrix x1;
		Matrix x2;
		Matrix xl;
		Matrix xu;

		double fx1, fx2;
		double fxL, fxU;
		double d;
		double goldRatio = (Math.sqrt(5.0)-1.0)/2.0;
		double error1, error2;		
	
		int size=x.getRowDimension();
		double tmpSum=0;
		Matrix nextGv;
		double lower, upper;	//boundaries of lamda
		double gradient, gradientLower, gradientUpper;
		double valueLower, valueUpper;

		double lamda=0.0;
		
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
		lamdaL = lower;
		lamdaU = upper;

		//initiate
		d = goldRatio * (lamdaU-lamdaL);
		lamda1 = lamdaL + d;
		lamda2 = lamdaU - d;
		x1 = x.plus(s.times(lamda1));
		x2 = x.plus(s.times(lamda2));
		xl = x.plus(s.times(lamdaL));
		xu = x.plus(s.times(lamdaU));			
		fx1 = objectiveFunction(x1);
		fx2 = objectiveFunction(x2);
		
		//iteration
		do{			
			if(fx1 < fx2){
				lamdaL = lamda2;
				d = goldRatio * (lamdaU-lamdaL);
				lamda2 = lamda1;
				lamda1 = lamdaL + d;
				x1 = x.plus(s.times(lamda1));
				x2 = x.plus(s.times(lamda2));
				fx1 = objectiveFunction(x1);
				fx2 = objectiveFunction(x2);
				error1 = checkConvergence(fx1,fx2);
				error2 = checkConvergence(gradientVector(x1));
			}
			else {
				lamdaU = lamda1;
				d = goldRatio * (lamdaU-lamdaL);
				lamda1 = lamda2;
				lamda2 = lamdaU - d;
				x1 = x.plus(s.times(lamda1));
				x2 = x.plus(s.times(lamda2));
				fx1 = objectiveFunction(x1);
				fx2 = objectiveFunction(x2);
				error1 = checkConvergence(fx2,fx1);
				error2 = checkConvergence(gradientVector(x2));
			}

		}while(error2>e2);	
		
		//return optimum lamda
		if(fx1<fx2) return lamda1;
		else return lamda2;
	}
	
	public double checkConvergence(double newFx, double oldFx){
		//Convergence Criteria: e1
		
		return Math.abs((newFx - oldFx)/oldFx);
	}

	public double checkConvergence(double[] gv){
		//Convergence Criteria: e2
		
		double max=0;
		for(int i=0 ; i<4 ; i++)
			if(Math.abs(gv[i]) > max) max = Math.abs(gv[i]); 
				
		return max;
	}
	
	public double checkConvergence(Matrix gv){
		//Convergence Criteria: e2
		
		double e=0;
		for(int i=0 ; i<4 ; i++) e += Math.pow(gv.get(i,0),2);
		e = Math.sqrt(e);
				
		return e;
	}
	
	
	public Matrix process(Matrix x, double E1, double E2){
		
		int i,iter;
		double newFx;
		double oldFx;
		Matrix newGv;
		Matrix oldGv;
		Matrix newS;
		Matrix oldS;
		double e1, e2;
		double lamda;
		
		oldFx = objectiveFunction(x);
		newGv = gradientVector(x);		
		newS = searchDirection(newGv);
		lamda = searchStepSize(x, newS, E2);
		x = x.plus(newS.times(lamda));
		newFx = objectiveFunction(x);
		e1 = checkConvergence(newFx, oldFx);
		e2 = checkConvergence(newGv);
		oldS = newS;
		iter=1;		
		
		//print the results
		//System.out.println(oldFx);
		System.out.println("iter\tx[0]\tx[1]\tx[2]\tx[3]"+
				"\tlamda\tnewFx\te1\te2");
		System.out.printf("%d\t%5.4f\t%5.4f\t%5.4f\t%5.4f\t%6.5f" +
				"\t%6.5f\t%6.5f\t%6.5f"
				,iter,x.get(0,0),x.get(1,0),x.get(2,0),x.get(3,0)
				,lamda,newFx,e1,e2);
		System.out.println();
		
		for(iter=2; e2>E2 ; iter++){
			oldFx = newFx;
			oldGv = newGv;
			oldS = newS;
			newGv = gradientVector(x);
			newS = searchDirection(oldS, oldGv, newGv);
			lamda = searchStepSize(x, newS, E2);
			x = x.plus(newS.times(lamda));
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

		double[][] x = { {-3.0}, {-1.0} , {-3.0}, {-1.0} };
		double e1 = 0.01;
		double e2;
		long startTime, endTime;
		double operationTime;
		
		FRmethod2 fr = new FRmethod2(); 
		
		//checking start time
		startTime = System.currentTimeMillis();
		

		
		System.out.println("a=0.01");
		e2 = 10.0*e1;
		fr.process(new Matrix(x), e1, e2);
		System.out.println();
		
		/*
		System.out.println("a=0.1");
		e2 = 0.1*e1;
		fr.process(x, e1, e2);
		System.out.println();
		
		System.out.println("a=1");
		e2 = e1;
		fr.process(x, e1, e2);
		System.out.println();
		
		System.out.println("a=10");
		e2 = 10.*e1;
		fr.process(x, e1, e2);
		*/
		
		//checking operation time
		endTime = System.currentTimeMillis();
		operationTime = (double)(endTime - startTime)/1000.0;

		System.out.println("Operation Time:"+operationTime+" sec");
	}
	
}
