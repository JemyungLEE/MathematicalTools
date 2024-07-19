package matrix;

public class MakeInverseMatrix {
	
	public MakeInverseMatrix(){}

	public MakeInverseMatrix(int size, double[][] K, double[] P, double[] D){
		GaussJordanMethod(size, K, P, D);
	}
	
	public void GaussJordanMethod(int size, 
									double[][] K, double[] P, double[] D){
		
		int i,j,k;
		double m;
		double[][] temp = new double[size][size+1];
		
		//make augmented matrix
		for(i=0 ; i<size ; i++){
			for(j=0 ; j<size ; j++) temp[i][j] = K[i][j];
			temp[i][size] = P[i];
		}
		
		//Gauss-Jordan Elimination
		for(i=0 ; i < size - 1 ; i++){
			for(j=i+1 ; j < size ; j++){
				m = temp[j][i] / temp[i][i];
				for(k=i ; k < size+1 ; k++){
					temp[j][k] = temp[j][k]- m * temp[i][k];
				}
			}
		}
 
		for(i= size - 1 ; i > 0  ; i--){    
			for(j=i-1 ; j >= 0 ; j--){
				m = temp[j][i] / temp[i][i];
				for(k= size - i ; k < size + 1 ; k++){
					temp[j][k] = temp[j][k] - m * temp[i][k];
				}
			}    
		}	
  
		for(i=0 ; i < size ; i++){
			m = 1.0 / temp[i][i];
			for(j=0 ; j < size+1 ; j++) temp[i][j] = temp[i][j] * m;
		}
  
		//write the result
		for(i=0 ; i < size ; i++) D[i] = temp[i][size];
		 
	}
}
