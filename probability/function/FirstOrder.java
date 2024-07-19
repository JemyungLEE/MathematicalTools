package probability.function;

public class FirstOrder {
	
	private double a0, a1;	// y = a0+a1*x
	private double r2;	//coefficient of determination
	
	public FirstOrder(){}
	
	public FirstOrder(double a0, double a1){
		setVariables(a0, a1);
	}
	
	public void setVariables(double a0, double a1){
		this.a0 = a0;
		this.a1 = a1;
	}
	
	public void setCoefficientOfDetermination(double r2){
		this.r2 = r2;
	}
	
	public double[] getVariables(){
		
		double[] variables = new double[2];
		
		variables[0] = a0;
		variables[1] = a1;
		
		return variables;
	}
	
	public double getValue(double x){
		return a0+a1*x;
	}

	public double getCoefficientOfDetermination(){
		return this.r2;
	}
	
}
