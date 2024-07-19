package optimization;

import Jama.Matrix;

public class BFGSmethod2 {

	/**
	 *  Subject: Structural Optimization
	 *  Developer: Jemyung Lee (ID: 2008-30334)
	 *    
	 *  Description: Broyden-Fletcher-Goldfarb-Shanno Method
	 */

	
	public double objectiveFunction(double[] x){
		
		//objective function
		return 100.0*Math.pow(x[1]-x[0]*x[0],2)+Math.pow(1.0-x[0]*x[0],2)+
				90.0*Math.pow(x[3]-x[2]*x[2],2)+Math.pow(1-x[2],2)+
				10.1*(Math.pow(x[1]-1.0, 2)+Math.pow(x[3]-1.0,2))+
				19.8*(x[1]-1.0)*(x[3]-1.0);
	}
	
	public Matrix gradientVector(double[] x){
		
		double[] gv = new double[4];	//gv[i-1]=df/dx[i]
		
		//partial derivative of objective function
		gv[0] = -400.0*(x[1]-x[0]*x[0])*x[0]-2.0+2.0*x[0];
		gv[1] = 1101.0/5.0*x[1]-200.0*x[0]*x[0]-40.0+99.0/5.0*x[3];
		gv[2] = -360.0*(x[3]-x[2]*x[2])*x[2]-2.0+2.0*x[2];
		gv[3] = 1001.0/5.0*x[3]-180.0*x[2]*x[2]-40.0+99.0/5.0*x[1];
			
		return new Matrix(gv,1).transpose();
	}
	
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
	
	public double searchStepSize(double[] x, Matrix S, double e1, double e2){
		
		double[] s = new double[4];
		for(int i=0 ; i<4 ; i++) s[i] = S.get(i,0);
		
		return goldenSectionSearch(x, s, e1, e2);
	}
	
	public double goldenSectionSearch(double[] x,double[] s,double e1,double e2){
		
		int i;
		double lamda1, lamda2;
		double lamdaL, lamdaU;
		double[] x1 = new double[4];
		double[] x2 = new double[4];
		double[] xl = new double[4];
		double[] xu = new double[4];

		double fx1, fx2;
		double fxL, fxU;
		double d;
		double goldRatio = (Math.sqrt(5.0)-1.0)/2.0;
		double error1, error2;
		
		//bracket
		double lamdaTmp;
		double dist = Math.sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]+x[3]*x[3]);
		//dist = 0.00001;
		lamdaL = 0;
		lamdaU = 0;
		
		do{
			lamdaTmp = lamdaL;
			lamdaL = lamdaU;			
			lamdaU += dist;
			for(i=0 ; i<4 ; i++){
				xl[i] = x[i] + lamdaL*s[i];
				xu[i] = x[i] + lamdaU*s[i];
			}		
			fxL = objectiveFunction(xl);
			fxU = objectiveFunction(xu);
		}while(fxU < fxL);
		lamdaL = lamdaTmp;
		
		System.out.println("lamdaLU:"+lamdaL+"\t"+lamdaU);
		System.out.println("fxLU:"+fxL+"\t"+fxU);
		
		//initiate
		d = goldRatio * (lamdaU-lamdaL);
		lamda1 = lamdaL + d;
		lamda2 = lamdaU - d;
		for(i=0 ; i<4 ; i++){
			x1[i] = x[i] + lamda1*s[i];
			x2[i] = x[i] + lamda2*s[i];
			xl[i] = x[i] + lamdaL*s[i];
			xu[i] = x[i] + lamdaU*s[i];
		}		
		fx1 = objectiveFunction(x1);
		fx2 = objectiveFunction(x2);
		//System.out.println("x:"+x1[0]+"\t"+x1[1]+"\t"+x1[2]+"\t"+x1[3]);
		//System.out.println("fx:"+fx1+"\t"+fx2);
		
		//iteration
		do{			
			if(fx1 < fx2){
				lamdaL = lamda2;
				d = goldRatio * (lamdaU-lamdaL);
				lamda2 = lamda1;
				lamda1 = lamdaL + d;
				for(i=0 ; i<4 ; i++){
					x1[i] = x[i] + lamda1*s[i];
					x2[i] = x[i] + lamda2*s[i];
				}
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
				for(i=0 ; i<4 ; i++){
					x1[i] = x[i] + lamda1*s[i];
					x2[i] = x[i] + lamda2*s[i];
				}
				fx1 = objectiveFunction(x1);
				fx2 = objectiveFunction(x2);
				error1 = checkConvergence(fx2,fx1);
				error2 = checkConvergence(gradientVector(x2));
			}

		}while(error1>e1 || error2>e2);	

		//System.out.println("x:"+x1[0]+"\t"+x1[1]+"\t"+x1[2]+"\t"+x1[3]);
		//System.out.println("lamda:"+lamda1+"\t"+lamda2);
		//System.out.println("fx:"+fx1+"\t"+fx2);
		//System.out.println("error;"+error1+"\t"+error2);
		
		//return optimum lamda
		if(fx1<fx2) return lamda1;
		else return lamda2;
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
	
	public double[] process(double[] x, double E1, double E2){
		
		int i,j,iter;
		double newFx, oldFx;
		Matrix newGv, oldGv;
		Matrix S;
		Matrix D;
		double[][] b = new double[4][4];
		Matrix B;
		Matrix g;
		double e1, e2;
		double lamda;
		
		//make b = I matrix
		for(i=0 ; i<4 ; i++){
			for(j=0 ; j<4 ; j++){
				if(i==j) b[i][j] = 1.0;
				else b[i][j] = 0.0;
			}
		}
		B = new Matrix(b);
		oldFx = objectiveFunction(x);
		newGv = gradientVector(x);
		S = searchDirection(newGv, B);
		lamda = searchStepSize(x, S, E1, E2);		
		x[0] += lamda*S.get(0,0);
		x[1] += lamda*S.get(1,0);
		newFx = objectiveFunction(x);
		e1 = checkConvergence(newFx, oldFx);
		e2 = checkConvergence(newGv);
		iter=1;
		
		System.out.printf("%d\t%5.4f\t%5.4f\t%5.4f" +
				"\t%5.4f\t%5.4f",iter,x[0],x[1],lamda,newFx,e1);
		System.out.println();
		
		for(iter=2; e1>E1 || e2>E2 ; iter++){
			oldFx = newFx;
			oldGv = newGv;
			newGv = gradientVector(x);
			g = makeGmatrix(oldGv, newGv);
			D = makeDmatrix(lamda, S);
			B = updateBmatrix(B, D, g);
			S = searchDirection(newGv, B);
			lamda = searchStepSize(x, S, E1, E2);
			x[0] += lamda*S.get(0,0);
			x[1] += lamda*S.get(1,0);
			newFx = objectiveFunction(x);
			e1 = checkConvergence(newFx, oldFx);
			e2 = checkConvergence(newGv);
			System.out.printf("%d\t%5.4f\t%5.4f\t%5.4f" +
					"\t%5.4f\t%5.4f",iter,x[0],x[1],lamda,newFx,e1);
			System.out.println();	
		}			
		
		return x;
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
		
		double[][] xx = { {-3.0}, {-1.0} , {-3.0}, {-1.0} };

		BFGSmethod2 bfgs = new BFGSmethod2(); 
		
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
