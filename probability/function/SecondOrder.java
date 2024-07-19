package probability.function;
public class SecondOrder {
	
	private double a0, a1, a2;	// y = a0+a1*x+a2*x^2
	private double r2;	//coefficient of determination

	public SecondOrder(){}
	
	public SecondOrder(double a0, double a1, double a2){
		setVariables(a0, a1, a2);
	}
	
	public void setVariables(double a0, double a1, double a2){
		this.a0 = a0;
		this.a1 = a1;
		this.a2 = a2;
	}
	
	public void setCoefficientOfDetermination(double r2){
		this.r2 = r2;
	}
	
	public double[] getVariables(){
		
		double[] variables = new double[3];
		
		variables[0] = a0;
		variables[1] = a1;
		variables[2] = a2;
		
		return variables;
	}
	
	public double getValue(double x){
		return a0+a1*x+a2*x*x;
	}

	public double getCoefficientOfDetermination(){
		return this.r2;
	}
	
}
