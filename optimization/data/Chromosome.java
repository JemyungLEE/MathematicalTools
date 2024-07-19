package optimization.data;

import java.util.ArrayList;

public class Chromosome {
	/**
	 *  Subject: Structural Optimization
	 *  Developer: Jemyung Lee (ID: 2008-30334)
	 *    
	 *  Description: Chromosome data structure for Genetic Algorithm
	 */	
	
	public int population;		//number of population
	public int chromosomeSize;	//length of chromosome
	public int[] chromosome; 	//chromosome string
	public ArrayList<int[]> chromosomeList;	//list of chromosome strings
	
	public int variables;		//number of variables
	public int[] variableLengths; //chromosome string length of each variable
	public int[] realnumber; //real number of each variable's chromosome string 
	public ArrayList<int[]> realnumberList;
	
	public Chromosome(){}
	
	public Chromosome(int population, int chromosomeSize){
		
		generate(population, chromosomeSize);
	}
	
	public Chromosome(ArrayList<int[]> list, int chromosomeSize){
		
		this.population = list.size();
		this.chromosomeSize = chromosomeSize;
		this.chromosomeList = list;
	}
	
	public void generate(int population, int chromosomeSize){
		
		int i,j;
		this.population = population;
		this.chromosomeSize = chromosomeSize;
		this.chromosomeList = new ArrayList<int[]>();
		for(i=0 ; i<population ; i++){
			chromosome = new int[chromosomeSize];
			for(j=0 ; j<chromosomeSize ; j++){
				if(Math.random()>0.5) chromosome[j] = 1;
				else chromosome[j] = 0;
			}
			chromosomeList.add(chromosome);
		}
	}
	
	public void setVariables(int size, int[] length){
		
		this.variables = size;
		this.variableLengths = length;
	}
	
	public void transform(){
		
		int i,j,k;
		int tmpPoint;
		int tmpValue;
		int[] tmpString;
		int[] tmpRealnumber;
		this.realnumberList = new ArrayList<int[]>();
		for(i=0 ; i<population ; i++){
			tmpPoint = 0;
			tmpString = chromosomeList.get(i);
			tmpRealnumber = new int[variables];
			for(j=0 ; j<variables ; j++){
				tmpValue=0;				
				for(k=0 ; k<variableLengths[j] ; k++){
					tmpValue += (int)Math.pow(2, k) * tmpString[tmpPoint]; 
					tmpPoint++;
				}
				tmpRealnumber[j] = tmpValue;
			}
			this.realnumberList.add(tmpRealnumber);
		}
	}
	
	public void transform(int size, int[] length){
		
		int i,j,k;
		int tmpPoint;
		int tmpValue;
		int[] tmpString;
		int[] tmpRealnumber;
		this.realnumberList = new ArrayList<int[]>();
		
		setVariables(size, length);
		for(i=0 ; i<population ; i++){
			tmpPoint = 0;
			tmpString = chromosomeList.get(i);
			tmpRealnumber = new int[variables];
			for(j=0 ; j<variables ; j++){
				tmpValue=0;				
				for(k=0 ; k<variableLengths[j] ; k++){
					tmpValue += (int)Math.pow(2, k) * tmpString[tmpPoint]; 
					tmpPoint++;
				}
				tmpRealnumber[j] = tmpValue;
			}
			this.realnumberList.add(tmpRealnumber);
		}
	}
}
