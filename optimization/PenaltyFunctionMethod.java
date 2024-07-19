package optimization;

import Jama.Matrix;

public class PenaltyFunctionMethod {
	/**
	 *  Subject: Structural Optimization
	 *  Developer: Jemyung Lee (ID: 2008-30334)
	 *    
	 *  Description: Interior and Exterior Penalty Function method
	 *  				with BFGS method and
	 *  				with Cubic Interpolation method
	 */

	public double Q(Matrix x){
				
		double x1 = x.get(0,0);
		double x2 = x.get(1,0);
		
		//objective function
		
		return Math.pow(x1,4)-2.0*x1*x1*x2+x1*x1+x1*x2*x2-2.0*x1+4.0;
	}

	public Matrix derivativeObjective(Matrix x){
		
		Matrix df = new Matrix(x.getRowDimension(),1);
		double x1 = x.get(0,0);
		double x2 = x.get(1,0);
		
		//derivative of objective function
		//odd: x1, even: x2
		df.set(0,0,4.0*(x1*x1*x1-x1*x2)+2.0*x1+x2*x2-2.0);
		df.set(1,0,-2.0*x1*x1+2.0*x1*x2);
		
		//x1:4*x1^3-4*x1*x2+2*x1+x2^2-2
		//x2:-2*x1^2+2*x1*x2
		return df;
	}
	
	public double equalityConstraint(Matrix x, double r){
		
		double x1 = x.get(0,0);
		double x2 = x.get(1,0);		
		
		//equality constraint function		
		return Math.pow(x1*x1+x2*x2-2.0,2);
	}
	
	public Matrix derivativeEquality(Matrix x){
		
		Matrix df = new Matrix(x.getRowDimension(),1);
		double x1 = x.get(0,0);
		double x2 = x.get(1,0);
		
		//derivative of equality constraints
		//odd: x1, even: x2
		double tmp = 4.0*(x1*x1+x2*x2-2.0);
		df.set(0,0,tmp*x1);
		df.set(1,0,tmp*x2);
		
		return df;
	}
	
	public Matrix inequalityConstraint(Matrix x){

		int constraints = 5;
		double[][] ineq = new double[constraints][1];
		double x1 = x.get(0,0);
		double x2 = x.get(1,0);	
		double tmp=0;
		
		//inequality constraint function		
		ineq[0][0]= 0.25*Math.pow(x1,2)+0.75*Math.pow(x2,2)-1.0;
		ineq[1][0]= x1-5.0;
		ineq[2][0]= -1.0*x1;
		ineq[3][0]= x2-5.0;
		ineq[4][0]= -1.0*x2;
		
		return new Matrix(ineq);
	}
	
	public Matrix derivativeInteriorInequality(Matrix x, double r){
		
		Matrix df = new Matrix(x.getRowDimension(),1);
		double x1 = x.get(0,0);
		double x2 = x.get(1,0);
		double[]tmpdf = new double[x.getRowDimension()*5];
		
		//derivative of inequality constraints
		//even: x1, odd: x2
		double tmp = -1.0/2.0/Math.pow((0.25*x1*x1+0.75*x2*x2-1.0),2);
		tmpdf[0]= tmp*x1;
		tmpdf[1]= 3.0*tmp*x2;
		tmpdf[2]= -1.0/Math.pow(x1-5.0,2);
		tmpdf[3]= 0.0;
		tmpdf[4]= 1.0/x1/x1;
		tmpdf[5]= 0.0;
		tmpdf[6]= 0.0;
		tmpdf[7]= -1.0/Math.pow(x2-5.0,2);
		tmpdf[8]= 0.0;
		tmpdf[9]= 1.0/x2/x2;
		
		df.set(0,0,tmpdf[0]+tmpdf[2]+tmpdf[4]+tmpdf[6]+tmpdf[8]);
		df.set(1,0,tmpdf[1]+tmpdf[3]+tmpdf[5]+tmpdf[7]+tmpdf[9]);
		
		return df;
	}
	
	public Matrix derivativeExteriorInequality(Matrix x, double r){
		
		Matrix ineq = inequalityConstraint(x);
		Matrix df = new Matrix(x.getRowDimension(),1);
		double x1 = x.get(0,0);
		double x2 = x.get(1,0);
		double[]tmpdf = new double[x.getRowDimension()*5];
		double tmpx1=0, tmpx2=0;
		
		//derivative of inequality constraints
		//even: x1, odd: x2
		double tmp = 0.25*x1*x1+0.75*x2*x2-1.0;
		tmpdf[0]= tmp*x1;
		tmpdf[1]= 3.0*tmp*x2;
		tmpdf[2]= 2.0*x1-10.0;
		tmpdf[3]= 0.0;
		tmpdf[4]= 2.0*x1;
		tmpdf[5]= 0.0;
		tmpdf[6]= 0.0;
		tmpdf[7]= 2.0*x2-10.0;
		tmpdf[8]= 0.0;
		tmpdf[9]= 2.0*x2;	
		
		for(int i=0 ; i<ineq.getRowDimension() ; i++){
			if(ineq.get(i, 0) > 0){
				tmpx1 += tmpdf[i*2];
				tmpx2 += tmpdf[i*2+1];
			}
		}
		
		df.set(0,0,tmpx1);
		df.set(1,0,tmpx2);
		
		return df;
	}
	
	public double interiorPenaltyFunction(Matrix x, double r){
		
		Matrix ineq = inequalityConstraint(x);
		double inequalityFunction = 0;
		double equalityFunction = 1.0/Math.sqrt(r)*equalityConstraint(x,r);
		
		for(int i=0 ; i<ineq.getRowDimension() ;i++) 
			inequalityFunction += 1.0/ineq.get(i,0);
		
		inequalityFunction *= -r;		
		
		return Q(x)+ inequalityFunction + equalityFunction;
	}
	
	public Matrix derivativeInteriorPenalty(Matrix x, double r){
		
		return derivativeObjective(x)
				.plus(derivativeInteriorInequality(x,r).times(-1.0*r))
				.plus(derivativeEquality(x).times(1.0/Math.sqrt(r)));
	}	
	
	public double exteriorPenaltyFunction(Matrix x, double r){
				
		Matrix ineq = inequalityConstraint(x);
		double inequalityFunction = 0;
		double equalityFunction = r*equalityConstraint(x,r);
		
		for(int i=0 ; i<ineq.getRowDimension() ;i++){
			if(ineq.get(i, 0)> 0) 
				inequalityFunction += ineq.get(i,0) * ineq.get(i,0); 
		}
		inequalityFunction *= r;	
		
		return Q(x)+ inequalityFunction + equalityFunction;
	}
	
	public Matrix derivativeExteriorPenalty(Matrix x, double r){
		
		return derivativeObjective(x)
				.plus(derivativeExteriorInequality(x,r).times(r))
				.plus(derivativeEquality(x).times(r));
	}
	
	public Matrix gradientVector(int type, Matrix x, double r){
		
		if(type==0) return derivativeInteriorPenalty(x,r);
		else return derivativeExteriorPenalty(x,r);
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
	
	public double searchStepSize(int type, Matrix x, Matrix s, 
											double r, double e2){
		
		return cubicInterpolationMethod(type, x, s, r, e2);
	}
	
	public double derivativeFunctionByLamda(int type, Matrix x, Matrix s,
											double r, double lamda){

		//calculate f'(lamda)
		Matrix nextGv = gradientVector(type, x.plus(s.times(lamda)),r); 
		return s.transpose().times(nextGv).get(0,0);
	}
	
	public double cubicInterpolationMethod(int type, Matrix x,Matrix s,
											double r, double e2){
	
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
		gradient = derivativeFunctionByLamda(type, x, s, r, lower);
		if(gradient<=0){
			upper = 0.1;
			gradient = derivativeFunctionByLamda(type, x, s, r, upper);
			while(gradient<0){
				upper *= 2.0;	//rescale the step-size
				gradient = derivativeFunctionByLamda(type, x, s, r, upper);				
			}			
		}
		else{			
			lower = -0.0001;
			gradient = derivativeFunctionByLamda(type, x, s, r, lower);
			while(gradient>0){
				lower *= 2.0;	//rescale the step-size
				gradient = derivativeFunctionByLamda(type, x, s, r, lower);
			}
			upper = lower/2.0;
			gradient = derivativeFunctionByLamda(type, x, s, r, upper);
			while(gradient<0){
				upper /= 2.0;	//rescale the step-size
				gradient = derivativeFunctionByLamda(type, x, s, r, upper);
			}
		}
		
		//calculate optimum lamda
		iter=0;
		do{
			valueLower = Q(x.plus(s.times(lower)));		
			valueUpper = Q(x.plus(s.times(upper)));
			gradientLower = derivativeFunctionByLamda(type, x, s, r, lower);
			gradientUpper = derivativeFunctionByLamda(type, x, s, r, upper);
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
			error2 = Math.abs(derivativeFunctionByLamda(type, x, s, r, lamda));
			tmpSum = 0;
			for(i=0 ; i<size ; i++) tmpSum += Math.pow(s.get(i,0), 2);
			tmpSum = Math.sqrt(tmpSum);
			error2 /= tmpSum;
			nextGv = gradientVector(type, x.plus(s.times(lamda)), r);
			tmpSum = 0;
			for(i=0 ; i<size ; i++) tmpSum += Math.pow(nextGv.get(i,0), 2);
			tmpSum = Math.sqrt(tmpSum);
			error2 /= tmpSum;			
			
			//update the boundaries
			gradient = derivativeFunctionByLamda(type, x, s, r, lamda);
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
		for(int i=0 ; i<gv.getRowDimension() ; i++)
						max += Math.pow(gv.get(i,0),2); 
				
		max = Math.sqrt(max);
		
		return max;
	}
	
	public Matrix BFGSmethod(int type,Matrix x, double r, double E1, double E2){
		//type=0 : interior penalty function method
		//type=1 : exterior penalty function method
		
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
		
		//make B = I matrix
		for(i=0 ; i<size ; i++){
			for(j=0 ; j<size ; j++){
				if(i==j) b[i][j] = 1.0;
				else b[i][j] = 0.0;
			}
		}		
		B = new Matrix(b);
		
		if(type==0) oldFx = interiorPenaltyFunction(x,r);
		else oldFx = exteriorPenaltyFunction(x,r);		
		newGv = gradientVector(type, x,r);
		S = searchDirection(newGv, B);
		lamda = searchStepSize(type, x, S, r, E2);
		x = x.plus(S.times(lamda));
		if(type==0) newFx = interiorPenaltyFunction(x,r);
		else newFx = exteriorPenaltyFunction(x,r);
		e1 = checkConvergence(newFx, oldFx);
		e2 = checkConvergence(newGv);
		iter=1;

		for(iter=2; e1>E1; iter++){
			oldFx = newFx;
			oldGv = newGv;
			newGv = gradientVector(type, x,r);
			g = makeGmatrix(oldGv, newGv);
			D = makeDmatrix(lamda, S);
			B = updateBmatrix(B, D, g);
			S = searchDirection(newGv, B);
			lamda = searchStepSize(type, x, S, r, E2);
			x = x.plus(S.times(lamda));
			if(type==0) newFx = interiorPenaltyFunction(x,r);
			else newFx = exteriorPenaltyFunction(x,r);
			e1 = checkConvergence(newFx, oldFx);
			e2 = checkConvergence(newGv);			
		}
		
		System.out.print(iter+"\t");
		
		return x;
	}
	
	public Matrix process(int type, Matrix x, double r, double c, 
											double e1, double e2){		
		//type=0 : interior penalty function method
		//type=1 : exterior penalty function method
		
		double previousF, nextF;
		double penaltyF=0;
		Matrix gradient;
		double convergence1, convergence2;
		
		nextF = Q(x);
		
		System.out.println("iter\tr\tx1\tx2\tpenalty\tf\te1\te2");
		do{
			previousF = nextF; 
			x = BFGSmethod(type,x,r,e1,e2);
			nextF = Q(x);
			if(type==0) penaltyF = interiorPenaltyFunction(x,r);
			else if(type==1) penaltyF = exteriorPenaltyFunction(x,r);
			gradient = derivativeObjective(x);
			
			convergence1 = checkConvergence(nextF, previousF);
			convergence2 = checkConvergence(gradient);			
			
			//print the result			
			System.out.printf("%f\t%6.5f\t%6.5f\t%6.5f\t%6.5f\t%6.5f\t%6.5f"
								,r,x.get(0,0),x.get(1,0),penaltyF,nextF
								,convergence1,convergence2);
			System.out.println();			
			
			//update r(t)
			r *= c;
			
		}while(convergence1>e1);
		
		return x;
	}
	
	public static void main(String[] args) {
		//type=0 : interior penalty function method
		//type=1 : exterior penalty function method		
		
		//set r(t) and c(t)
		double r = 10.0;
		double c = 0.1;
		
		//set convergence criteria
		double e1 = 0.0001;
		double e2 = 0.01;
		long startTime, endTime;
		double operationTime;
		
		///process
		PenaltyFunctionMethod pfm = new PenaltyFunctionMethod(); 

		//operate 'interior' penalty method		
		startTime = System.currentTimeMillis();	//check start time
		r = 10.0;	//set r(t)
		c = 0.1;	//set c(t)
		double[][] in = { {3.0}, {3.0} };	//set initial value
		System.out.println("interior penalty method");
		pfm.process(0, new Matrix(in), r, c, e1, e2);

		//check operation time
		endTime = System.currentTimeMillis();
		operationTime = (double)(endTime - startTime)/1000.0;
		System.out.println("Operation Time:"+operationTime+" sec");
		
		//operate 'exterior' penalty method
		startTime = System.currentTimeMillis();	//check start time
		r = 1.0;	//set r(t)
		c = 10;	//set c(t)
		double[][] ex = { {-0.5}, {-0.5} };	//set initial value
		System.out.println("exterior penalty method");
		pfm.process(1, new Matrix(ex), r, c, e1, e2);
		
		//check operation time
		endTime = System.currentTimeMillis();
		operationTime = (double)(endTime - startTime)/1000.0;
		System.out.println("Operation Time:"+operationTime+" sec");
	}
}
