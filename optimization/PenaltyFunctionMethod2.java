package optimization;

import Jama.Matrix;

public class PenaltyFunctionMethod2 {
	/**
	 *  Subject: Structural Optimization
	 *  Developer: Jemyung Lee (ID: 2008-30334)
	 *    
	 *  Description: Interior and Exterior Penalty Function method
	 *  				with BFGS method and
	 *  				with Cubic Interpolation method
	 */

	public double objectiveFunction(Matrix x){
				
		double x1 = x.get(0,0);
		double x2 = x.get(1,0);
		double x3 = x.get(2,0);
		double x4 = x.get(3,0);
		
		//objective function		
		return 2.0*x2*x4+(x1-2.0*x4)*x3;
	}

	public Matrix derivativeObjective(Matrix x){
		
		Matrix df = new Matrix(x.getRowDimension(),1);
		double x1 = x.get(0,0);
		double x2 = x.get(1,0);
		double x3 = x.get(2,0);
		double x4 = x.get(3,0);
		
		//derivative of objective function
		df.set(0,0,x3);
		df.set(1,0,2.0*x4);
		df.set(2,0,x1-2.0*x4);
		df.set(3,0,2.0*x2-2.0*x3);
		
		return df;
	}
	
	public double momentOfInertia(Matrix x){
		double x1 = x.get(0,0);
		double x2 = x.get(1,0);
		double x3 = x.get(2,0);
		double x4 = x.get(3,0);		
		
		return (x2*Math.pow(x1,3)-(x2-x4)*Math.pow(x1-2.0*x3,3))/12.0;
	}
	
	public double firstMoment(Matrix x){
		double x1 = x.get(0,0);
		double x2 = x.get(1,0);
		double x3 = x.get(2,0);
		double x4 = x.get(3,0);		
		
		return 0.5*(x2*x4+x3*(x1-x4))*(x1-x4);
	}
	
	public Matrix inequalityConstraint(Matrix x){

		int constraints = 15;
		double[][] ineq = new double[constraints][1];
		double x1 = x.get(0,0);
		double x2 = x.get(1,0);
		double x3 = x.get(2,0);
		double x4 = x.get(3,0);
		double Ix = momentOfInertia(x);		
		double Qx = firstMoment(x);
		
		//inequality constraint function		
		ineq[0][0]= 0.0000205078125/Ix-6.0;
		ineq[1][0]= 3.9375*x1/Ix-250000;
		ineq[2][0]= 6.3*Qx/Ix/x3-145000;
		ineq[3][0]= x1-3.0*x2;
		ineq[4][0]= 2.0*x2-x1;
		ineq[5][0]= x3-1.5*x4;
		ineq[6][0]= 0.5*x4-x3;
		ineq[7][0]= 0.075-x1;
		ineq[8][0]= x1-0.5;
		ineq[9][0]= 0.050-x2;
		ineq[10][0]= x2-0.380;
		ineq[11][0]= 0.003-x3;
		ineq[12][0]= x3-0.019;
		ineq[13][0]= 0.006-x4;
		ineq[14][0]= x4-0.032;
		
		return new Matrix(ineq);
	}
	
	public Matrix derivativeInteriorInequality(Matrix x, double r){
		
		int i,j;
		int constraint = 15;
		Matrix ineq = inequalityConstraint(x);
		Matrix df = new Matrix(x.getRowDimension(),1);
		double x1 = x.get(0,0);
		double x2 = x.get(1,0);
		double x3 = x.get(2,0);
		double x4 = x.get(3,0);
		double[][]tmpdf = new double[constraint][x.getRowDimension()];
		double tmp;
		
		//derivative of inequality constraints
		tmpdf[0][0]= 1.0/Math.pow((0.0000205078125/(1.0/12.0*x2*x1*x1*x1-1.0/12.0*(x2-x4)*Math.pow((x1-2.0*x3),3))-6.0),2)*0.0000205078125/Math.pow((1.0/12.0*x2*x1*x1*x1-1.0/12.0*(x2-x4)*Math.pow((x1-2.0*x3),3)),2)*(1.0/4.0*x2*x1*x1-1.0/4.0*(x2-x4)*(x1-2*x3)*(x1-2*x3));
		tmpdf[0][1]= 1.0/Math.pow((0.0000205078125/(1.0/12.0*x2*x1*x1*x1-1.0/12.0*(x2-x4)*Math.pow((x1-2.0*x3),3))-6.0),2)*0.0000205078125/Math.pow((1.0/12.0*x2*x1*x1*x1-1.0/12.0*(x2-x4)*Math.pow((x1-2.0*x3),3)),2)*(1.0/12.0*x1*x1*x1-1.0/12.0*Math.pow((x1-2.0*x3),3));
		tmpdf[0][2]= 1.0/2.0/Math.pow((0.0000205078125/(1.0/12.0*x2*x1*x1*x1-1.0/12.0*(x2-x4)*Math.pow((x1-2.0*x3),3))-6.0),2)*0.0000205078125/Math.pow((1.0/12.0*x2*x1*x1*x1-1.0/12.0*(x2-x4)*Math.pow((x1-2.0*x3),3)),2)*(x2-x4)*Math.pow((x1-2*x3),2);
		tmpdf[0][3]= 1.0/12.0/Math.pow((0.0000205078125/(1.0/12.0*x2*x1*x1*x1-1.0/12.0*(x2-x4)*Math.pow((x1-2.0*x3),3))-6.0),2)*0.0000205078125/Math.pow((1.0/12.0*x2*x1*x1*x1-1.0/12.0*(x2-x4)*Math.pow((x1-2.0*x3),3)),2)*Math.pow((x1-2.0*x3),3);
		tmpdf[1][0]= -1.0/Math.pow((0.0000205078125*x1/(1.0/12.0*x2*x1*x1*x1-1.0/12.0*(x2-x4)*Math.pow((x1-2.0*x3),3))-250000.0),2)*(0.0000205078125/(1.0/12.0*x2*x1*x1*x1-1.0/12.0*(x2-x4)*Math.pow((x1-2.0*x3),3))-0.0000205078125*x1/Math.pow((1.0/12.0*x2*x1*x1*x1-1.0/12.0*(x2-x4)*Math.pow((x1-2*x3),3)),2)*(1.0/4.0*x2*x1*x1-1.0/4.0*(x2-x4)*Math.pow((x1-2.0*x3),2)));
		tmpdf[1][1]= 1.0/Math.pow((3.9375*x1/(1.0/12.0*x2*x1*x1*x1-1.0/12.0*(x2-x4)*Math.pow((x1-2.0*x3),3))-250000.0),2)*3.9375*x1/Math.pow((1.0/12.0*x2*x1*x1*x1-1.0/12.0*(x2-x4)*Math.pow((x1-2.0*x3),3)),2)*(1.0/12.0*x1*x1*x1-1.0/12.0*Math.pow((x1-2*x3),3));
		tmpdf[1][2]= 1.0/2.0/Math.pow((3.9375*x1/(1.0/12.0*x2*x1*x1*x1-1.0/12.0*(x2-x4)*Math.pow((x1-2.0*x3),3))-250000.0),2)*3.9375*x1/Math.pow((1.0/12.0*x2*x1*x1*x1-1.0/12.0*(x2-x4)*Math.pow((x1-2.0*x3),3)),2)*(x2-x4)*Math.pow((x1-2*x3),2);
		tmpdf[1][3]= 1.0/12.0/Math.pow((3.9375*x1/(1.0/12.0*x2*x1*x1*x1-1.0/12.0*(x2-x4)*Math.pow((x1-2.0*x3),3))-250000.0),2)*3.9375*x1/Math.pow((1.0/12.0*x2*x1*x1*x1-1.0/12.0*(x2-x4)*Math.pow((x1-2.0*x3),3)),2)*Math.pow((x1-2.0*x3),3);
		tmpdf[2][0]= -1.0/Math.pow((6.3*(1.0/2.0*x2*x4*(x1-x4)+1.0/2.0*x3*Math.pow((x1-x4),2))/x3/(1.0/12.0*x2*x1*x1*x1-1.0/12.0*(x2-x4)*Math.pow((x1-2.0*x3),3))-145000.0),2)*(6.3*(1.0/2.0*x2*x4+x3*(x1-x4))/x3/(1.0/12.0*x2*x1*x1*x1-1.0/12.0*(x2-x4)*Math.pow((x1-2.0*x3),3))-6.3*(1.0/2.0*x2*x4*(x1-x4)+1.0/2.0*x3*Math.pow((x1-x4),2))/x3/Math.pow((1.0/12.0*x2*x1*x1*x1-1.0/12.0*(x2-x4)*Math.pow((x1-2.0*x3),3)),2)*(1.0/4.0*x2*x1*x1-1.0/4.0*(x2-x4)*Math.pow((x1-2.0*x3),2)));
		tmpdf[2][1]= -1.0/Math.pow((6.3*(1.0/2.0*x2*x4*(x1-x4)+1.0/2.0*x3*Math.pow((x1-x4),2))/x3/(1.0/12.0*x2*x1*x1*x1-1.0/12.0*(x2-x4)*Math.pow((x1-2.0*x3),3))-145000.0),2)*(1/2*6.3*x4*(x1-x4)/x3/(1.0/12.0*x2*x1*x1*x1-1.0/12.0*(x2-x4)*Math.pow((x1-2.0*x3),3))-6.3*(1.0/2.0*x2*x4*(x1-x4)+1.0/2.0*x3*Math.pow((x1-x4),2))/x3/Math.pow((1.0/12.0*x2*x1*x1*x1-1.0/12.0*(x2-x4)*Math.pow((x1-2.0*x3),3)),2)*(1.0/12.0*x1*x1*x1-1.0/12.0*Math.pow((x1-2.0*x3),3)));
		tmpdf[2][2]= -1.0/Math.pow((6.3*(1.0/2.0*x2*x4*(x1-x4)+1.0/2.0*x3*Math.pow((x1-x4),2))/x3/(1.0/12.0*x2*x1*x1*x1-1.0/12.0*(x2-x4)*Math.pow((x1-2.0*x3),3))-145000.0),2)*(1/2*6.3*Math.pow((x1-x4),2)/x3/(1.0/12.0*x2*x1*x1*x1-1.0/12.0*(x2-x4)*Math.pow((x1-2.0*x3),3))-6.3*(1.0/2.0*x2*x4*(x1-x4)+1.0/2.0*x3*Math.pow((x1-x4),2))/x3/x3/(1.0/12.0*x2*x1*x1*x1-1.0/12.0*(x2-x4)*Math.pow((x1-2.0*x3),3))-1.0/2.0*6.3*(1.0/2.0*x2*x4*(x1-x4)+1.0/2.0*x3*Math.pow((x1-x4),2))/x3/Math.pow((1.0/12.0*x2*x1*x1*x1-1.0/12.0*(x2-x4)*Math.pow((x1-2.0*x3),3)),2)*(x2-x4)*Math.pow((x1-2.0*x3),2));
		tmpdf[2][3]= -1.0/Math.pow((6.3*(1.0/2.0*x2*x4*(x1-x4)+1.0/2.0*x3*Math.pow((x1-x4),2))/x3/(1.0/12.0*x2*x1*x1*x1-1.0/12.0*(x2-x4)*Math.pow((x1-2.0*x3),3))-145000.0),2)*(6.3*(1/2*x2*(x1-x4)-1/2*x2*x4-x3*(x1-x4))/x3/(1.0/12.0*x2*x1*x1*x1-1.0/12.0*(x2-x4)*Math.pow((x1-2.0*x3),3))-1/12*6.3*(1.0/2.0*x2*x4*(x1-x4)+1.0/2.0*x3*Math.pow((x1-x4),2))/x3/Math.pow((1.0/12.0*x2*x1*x1*x1-1.0/12.0*(x2-x4)*Math.pow((x1-2.0*x3),3)),2)*Math.pow((x1-2*x3),3));
		tmpdf[3][0]= -1.0/Math.pow((x1-3.0*x2),2);
		tmpdf[3][1]= 3.0/Math.pow((x1-3.0*x2),2);
		tmpdf[3][2]= 0.0;
		tmpdf[3][3]= 0.0;
		tmpdf[4][0]= 1.0/Math.pow((2.0*x2-x1),2);
		tmpdf[4][1]= -2.0/Math.pow((2.0*x2-x1),2);
		tmpdf[4][2]= 0.0;
		tmpdf[4][3]= 0.0;
		tmpdf[5][0]= 0.0;
		tmpdf[5][1]= 0.0;
		tmpdf[5][2]= -1.0/Math.pow((x3-3.0/2.0*x4),2);
		tmpdf[5][3]= 3.0/2.0/Math.pow((x3-3.0/2.0*x4),2);
		tmpdf[6][0]= 0.0;
		tmpdf[6][1]= 0.0;
		tmpdf[6][2]= -1.0/Math.pow((1.0/2.0*x4+x3),2);
		tmpdf[6][3]= -1.0/2.0/Math.pow((1.0/2.0*x4+x3),2);
		tmpdf[7][0]= 1.0/Math.pow((0.075-x1),2);
		tmpdf[7][1]= 0.0;
		tmpdf[7][2]= 0.0;
		tmpdf[7][3]= 0.0;
		tmpdf[8][0]= -1.0/Math.pow((x1-0.500),2);
		tmpdf[8][1]= 0.0;
		tmpdf[8][2]= 0.0;
		tmpdf[8][3]= 0.0;
		tmpdf[9][0]= 0.0;
		tmpdf[9][1]= 1.0/Math.pow((0.050-x2),2);
		tmpdf[9][2]= 0.0;
		tmpdf[9][3]= 0.0;
		tmpdf[10][0]= 0.0;
		tmpdf[10][1]= -1.0/Math.pow((x2-0.380),2);
		tmpdf[10][2]= 0.0;
		tmpdf[10][3]= 0.0;
		tmpdf[11][0]= 0.0;
		tmpdf[11][1]= 0.0;
		tmpdf[11][2]= 1.0/Math.pow((0.003-x3),2);
		tmpdf[11][3]= 0.0;
		tmpdf[12][0]= 0.0;
		tmpdf[12][1]= 0.0;
		tmpdf[12][2]= -1.0/Math.pow((x3-0.019),2);
		tmpdf[12][3]= 0.0;
		tmpdf[13][0]= 0.0;
		tmpdf[13][1]= 0.0;
		tmpdf[13][2]= 0.0;
		tmpdf[13][3]= 1.0/Math.pow((0.006-x4),2);
		tmpdf[14][0]= 0.0;
		tmpdf[14][1]= 0.0;
		tmpdf[14][2]= 0.0;
		tmpdf[14][3]= -1.0/Math.pow((x4-0.032),2);
		
		for(i=0 ; i<x.getRowDimension() ; i++){
			tmp = 0;
			for(j=0 ; j<constraint ; j++) tmp += tmpdf[j][i];
			df.set(i, 0, tmp);
		}
		
		return df;
	}
	
	public Matrix derivativeExteriorInequality(Matrix x, double r){
		
		int constraint = 15;
		Matrix ineq = inequalityConstraint(x);
		Matrix df = new Matrix(x.getRowDimension(),1);
		double x1 = x.get(0,0);
		double x2 = x.get(1,0);
		double x3 = x.get(2,0);
		double x4 = x.get(3,0);
		double[][]tmpdf = new double[constraint][x.getRowDimension()];
		double tmpx1=0, tmpx2=0;
		
		//derivative of inequality constraints
		//even: x1, odd: x2
		/*
		double tmp = 0.25*x1*x1+0.75*x2*x2-1.0;
		tmpdf[0][0]= ;
		tmpdf[0][1]= ;
		tmpdf[0][2]= ;
		tmpdf[0][3]= ;
		tmpdf[1][0]= ;
		tmpdf[1][1]= ;
		tmpdf[1][2]= ;
		tmpdf[1][3]= ;
		tmpdf[2][0]= ;
		tmpdf[2][1]= ;
		tmpdf[2][2]= ;
		tmpdf[2][3]= ;
		tmpdf[3][0]= ;
		tmpdf[3][1]= ;
		tmpdf[3][2]= ;
		tmpdf[3][3]= ;
		tmpdf[4][0]= ;
		tmpdf[4][1]= ;
		tmpdf[4][2]= ;
		tmpdf[4][3]= ;
		tmpdf[5][0]= ;
		tmpdf[5][1]= ;
		tmpdf[5][2]= ;
		tmpdf[5][3]= ;
		tmpdf[6][0]= ;
		tmpdf[6][1]= ;
		tmpdf[6][2]= ;
		tmpdf[6][3]= ;
		tmpdf[7][0]= ;
		tmpdf[7][1]= ;
		tmpdf[7][2]= ;
		tmpdf[7][3]= ;
		tmpdf[8][0]= ;
		tmpdf[8][1]= ;
		tmpdf[8][2]= ;
		tmpdf[8][3]= ;
		tmpdf[9][0]= ;
		tmpdf[9][1]= ;
		tmpdf[9][2]= ;
		tmpdf[9][3]= ;
		tmpdf[10][0]= ;
		tmpdf[10][1]= ;
		tmpdf[10][2]= ;
		tmpdf[10][3]= ;
		tmpdf[11][0]= ;
		tmpdf[11][1]= ;
		tmpdf[11][2]= ;
		tmpdf[11][3]= ;
		tmpdf[12][0]= ;
		tmpdf[12][1]= ;
		tmpdf[12][2]= ;
		tmpdf[12][3]= ;
		tmpdf[13][0]= ;
		tmpdf[13][1]= ;
		tmpdf[13][2]= ;
		tmpdf[13][3]= ;
		tmpdf[14][0]= ;
		tmpdf[14][1]= ;
		tmpdf[14][2]= ;
		tmpdf[14][3]= ;
		
		for(int i=0 ; i<ineq.getRowDimension() ; i++){
			if(ineq.get(i, 0) > 0){
				tmpx1 += tmpdf[i*2];
				tmpx2 += tmpdf[i*2+1];
			}
		}
		*/
		df.set(0,0,0.0);
		df.set(1,0,0.0);
		
		return df;
	}
	
	public double interiorPenaltyFunction(Matrix x, double r){
		
		Matrix ineq = inequalityConstraint(x);
		double inequalityFunction = 0;
		
		for(int i=0 ; i<ineq.getRowDimension() ;i++) 
			inequalityFunction += 1.0/ineq.get(i,0);
		
		inequalityFunction *= -r;		
		
		return objectiveFunction(x)+ inequalityFunction;
	}
	
	public Matrix derivativeInteriorPenalty(Matrix x, double r){
		
		return derivativeObjective(x)
				.plus(derivativeInteriorInequality(x,r).times(-1.0*r));
	}	
	
	public double exteriorPenaltyFunction(Matrix x, double r){
				
		Matrix ineq = inequalityConstraint(x);
		double inequalityFunction = 0;
		double equalityFunction = r;
		
		for(int i=0 ; i<ineq.getRowDimension() ;i++){
			if(ineq.get(i, 0)> 0) 
				inequalityFunction += ineq.get(i,0) * ineq.get(i,0); 
		}
		inequalityFunction *= r;	
		
		return objectiveFunction(x)+ inequalityFunction + equalityFunction;
	}
	
	public Matrix derivativeExteriorPenalty(Matrix x, double r){
		
		return derivativeObjective(x)
				.plus(derivativeExteriorInequality(x,r).times(r));
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
			valueLower = objectiveFunction(x.plus(s.times(lower)));		
			valueUpper = objectiveFunction(x.plus(s.times(upper)));
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
		
		nextF = objectiveFunction(x);
		
		System.out.println("iter\tr\tx1\tx2\tpenalty\tf\te1\te2");
		do{
			previousF = nextF; 
			x = BFGSmethod(type,x,r,e1,e2);
			nextF = objectiveFunction(x);
			if(type==0) penaltyF = interiorPenaltyFunction(x,r);
			else if(type==1) penaltyF = exteriorPenaltyFunction(x,r);
			gradient = derivativeObjective(x);
			
			convergence1 = checkConvergence(nextF, previousF);
			convergence2 = checkConvergence(gradient);			
			
			//print the result			
			System.out.printf("%f\t%6.5f\t%6.5f\t%6.5f\t%6.5f\t%6.5f\t%6.5f\t%6.5f\t%6.5f"
								,r,x.get(0,0),x.get(1,0),x.get(2,0),x.get(3,0),penaltyF,nextF
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
		double e2 = 0.001;
		long startTime, endTime;
		double operationTime;
		
		///process
		PenaltyFunctionMethod2 pfm = new PenaltyFunctionMethod2(); 

		//operate 'interior' penalty method		
		startTime = System.currentTimeMillis();	//check start time
		r = 10.0;	//set r(t)
		c = 0.1;	//set c(t)
		double[][] in = { {0.1}, {0.1} , {0.01} , {0.01} };	//set initial value
		System.out.println("interior penalty method");
		pfm.process(0, new Matrix(in), r, c, e1, e2);

		//check operation time
		endTime = System.currentTimeMillis();
		operationTime = (double)(endTime - startTime)/1000.0;
		System.out.println("Operation Time:"+operationTime+" sec");
		/*
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
		*/
	}
}
