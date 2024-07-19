package probability.data;

import java.io.File;
import java.util.ArrayList;
import java.util.Scanner;

public class TwoDimensionData {

	private int dataNumber;
	private ArrayList<Double> x;
	private ArrayList<Double> y;
	
	public TwoDimensionData(){
		initiate();
	}
	
	public TwoDimensionData(String source){
		initiate();
		readData(source);
	}
	
	public void initiate(){
		x = new ArrayList<Double>();
		y = new ArrayList<Double>();
	}
	
	public void addValueX(double value){
		x.add(value);
	}
	
	public void addValueY(double value){
		y.add(value);
	}
	
	public void setDataNumber(int number){
		this.dataNumber = number;
	}
	
	public double getValueX(int index){
		return x.get(index);
	}
	
	public double getValueY(int index){
		return y.get(index);
	}
	
	public int getDataNumber(){
		return this.dataNumber;
	}
	
	public void readData(String source){
		
		int cnt = 0;
		
		try{
			File file = new File(source);
			Scanner scan = new Scanner(file);
		
			if(scan.next().equals("x")){
				if(scan.next().equals("y")){
					while(scan.hasNextDouble()){
						this.addValueX(scan.nextDouble());
						this.addValueY(scan.nextDouble());
						cnt++;
					}
					this.setDataNumber(cnt);
				}
				else System.err.println("Wrong Input Data.");				
			}
			else System.err.println("Wrong Input Data.");			
			
			scan.close();	
		
		} catch (Exception e) {
			System.out.println("error");
		}
	}
}
