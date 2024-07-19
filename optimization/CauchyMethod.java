package optimization;

public class CauchyMethod {

	/**
	 *  Subject: Structural Optimization
	 *  Developer: Jemyung Lee (ID: 2008-30334)
	 *    
	 *  Description: Cauchy (Steepest descent) Method
	 */

	
	public double objectiveFunction(double[] x){
		
		double x1 = x[0];
		double x2 = x[1];
		
		double dx1 = Math.pow(x1, 2);	//x1^2
		double dx2 = Math.pow(x2, 2);	//x2^2
		double qx1 = Math.pow(x1, 4);	//x1^4
		double qx2 = Math.pow(x2, 4);	//x2^4
		
		//objective function
		return 0.1*(12.0+dx1+(1.0+dx2)/dx1+(dx1+dx2+100.0)/(qx1*qx2));
	}
	
	public double[] gradientVector(double[] x){
		
		double x1 = x[0];
		double x2 = x[1];
		
		double[] gv = new double[2];	//gv[0]=df/dx1, gv[1]=df/dx2
		double dx1 = Math.pow(x1, 2);	//x1^2
		double tx1 = Math.pow(x1, 3);	//x1^3
		double qx1 = Math.pow(x1, 4);	//x1^4
		double fx1 = Math.pow(x1, 5);	//x1^5
		double dx2 = Math.pow(x2, 2);	//x2^2
		double tx2 = Math.pow(x2, 3);	//x2^3
		double qx2 = Math.pow(x2, 4);	//x2^4
		double fx2 = Math.pow(x2, 5);	//x2^5
		
		//partial derivative of objective function
		gv[0] = 0.2*(x1-(1.0+dx2)/tx1-1.0/tx1/qx2-2.0/fx1/dx2-200.0/fx1/qx2);
		gv[1] = 0.2*(x2/dx1-2.0/dx1/fx2-1.0/qx1/tx2-200.0/qx1/fx2);
		
		return gv;
	}
	
	public double[] searchDirection(double[] x){
		
		double[] gv = gradientVector(x);
		
		double[] s = new double[2];
		s[0] = -1.0*gv[0];
		s[1] = -1.0*gv[1];
		
		return s;
	}
	
	public double searchStepSize(double[] x, double[] s, double e){
	
		return goldenSectionSearch(x, s, e);
	}
	
	public double goldenSectionSearch(double[] x, double[] s, double e){
		
		int i;
		double lamda1, lamda2;
		double lamdaL, lamdaU;
		double[] x1 = new double[2];
		double[] x2 = new double[2];
		double[] xl = new double[2];
		double[] xu = new double[2];

		double fx1, fx2;
		double fxL, fxU;
		double d;
		double goldRatio = (Math.sqrt(5.0)-1.0)/2.0;
		double error;
		
		//bracket
		double lamdaTmp;
		double dist = Math.sqrt(x[0]*x[0]+x[1]*x[1]);
		lamdaL = 0;
		lamdaU = 0;
		do{
			lamdaTmp = lamdaL;
			lamdaL = lamdaU;			
			lamdaU += dist;
			for(i=0 ; i<2 ; i++){
				xl[i] = x[i] + lamdaL*s[i];
				xu[i] = x[i] + lamdaU*s[i];
			}		
			fxL = objectiveFunction(xl);
			fxU = objectiveFunction(xu);
		}while(fxU < fxL);
		lamdaL = lamdaTmp;
		
		//initiate
		d = goldRatio * (lamdaU-lamdaL);
		lamda1 = lamdaL + d;
		lamda2 = lamdaU - d;
		for(i=0 ; i<2 ; i++){
			x1[i] = x[i] + lamda1*s[i];
			x2[i] = x[i] + lamda2*s[i];
			xl[i] = x[i] + lamdaL*s[i];
			xu[i] = x[i] + lamdaU*s[i];
		}		
		fx1 = objectiveFunction(x1);
		fx2 = objectiveFunction(x2);
		
		//iteration
		do{			
			if(fx1 < fx2){
				lamdaL = lamda2;
				d = goldRatio * (lamdaU-lamdaL);
				lamda2 = lamda1;
				lamda1 = lamdaL + d;
				for(i=0 ; i<2 ; i++){
					x1[i] = x[i] + lamda1*s[i];
					x2[i] = x[i] + lamda2*s[i];
				}
				fx1 = objectiveFunction(x1);
				fx2 = objectiveFunction(x2);
				error = checkConvergence(fx1,fx2);
			}
			else {
				lamdaU = lamda1;
				d = goldRatio * (lamdaU-lamdaL);
				lamda1 = lamda2;
				lamda2 = lamdaU - d;
				for(i=0 ; i<2 ; i++){
					x1[i] = x[i] + lamda1*s[i];
					x2[i] = x[i] + lamda2*s[i];
				}
				fx1 = objectiveFunction(x1);
				fx2 = objectiveFunction(x2);
				error = checkConvergence(fx2,fx1);
			}
			
		}while(error>e);	
		
		//return optimum lamda
		if(fx1<fx2) return lamda1;
		else return lamda2;
	}	
	
	public double checkConvergence(double newFx, double oldFx){
		//Convergence Criteria: e1
		
		return Math.abs((newFx - oldFx)/oldFx);
	}
	
	public double[] process(double[] x, double e){
		
		int iter;
		double newFx;
		double oldFx;
		double e1;
		double[] s = new double[2];
		double lamda;
		
		oldFx = objectiveFunction(x);
		s = searchDirection(x);
		lamda = searchStepSize(x, s, e);
		x[0] += lamda*s[0];
		x[1] += lamda*s[1];
		newFx = objectiveFunction(x);
		e1 = checkConvergence(newFx, oldFx);
		iter=1;
		
		System.out.printf("%d\t%5.4f\t%5.4f\t%5.4f" +
				"\t%5.4f\t%5.4f",iter,x[0],x[1],lamda,newFx,e1);
		System.out.println();
		
		for(iter=2; e1>e ; iter++){
			oldFx = newFx;
			s = searchDirection(x);
			lamda = searchStepSize(x, s, e);
			x[0] += lamda*s[0];
			x[1] += lamda*s[1];
			newFx = objectiveFunction(x);
			e1 = checkConvergence(newFx, oldFx);		
			System.out.printf("%d\t%5.4f\t%5.4f\t%5.4f" +
					"\t%5.4f\t%5.4f",iter,x[0],x[1],lamda,newFx,e1);
			//System.out.printf("iter:\t%d\tx1:\t%5.4f\tx2:\t%5.4f\tlamda:\t%5.4f"
			//		+"\tnewFx:\t%5.4f\terror:\t%5.4f"
			//		,iter,x[0],x[1],lamda,newFx,e1);
			System.out.println();	
		}			
		
		return x;
	}
		
	public static void main(String[] args) {

		double[] x = { 1.0 , 1.0 };
		double e = 0.0001;
		long startTime, endTime;
		double operationTime;
		
		x[0] = 1.0;
		x[1] = 1.0;
		
		CauchyMethod cm = new CauchyMethod(); 
		
		//checking start time
		startTime = System.currentTimeMillis();
		
		cm.process(x, e);
		
		//checking operation time
		endTime = System.currentTimeMillis();
		operationTime = (double)(endTime - startTime)/1000.0;

		System.out.println("Operation Time:"+operationTime+" sec");
	}

}
