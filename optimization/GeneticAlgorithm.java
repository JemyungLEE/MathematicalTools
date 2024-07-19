package optimization;

import java.util.ArrayList;
import Jama.Matrix;

import optimization.data.Chromosome;
import optimization.function.NewmarkMethod;

public class GeneticAlgorithm {
	/**
	 *  Subject: Structural Optimization
	 *  Developer: Jemyung Lee (ID: 2008-30334)
	 *    
	 *  Description: Genetic Algorithm
	 */
	
	public double[] objectiveFunction(Chromosome list){
		
		int i,j;
		double[] value = new double[list.population];
		double min=2.1*Math.pow(10, 6);
		double max=2.1*Math.pow(10, 10);
		int length=list.variableLengths[0];
		double r = Math.pow(2, length);
		double resolution = (max-min)/(r-1.0);
		double[][] k = new double[list.population][list.variables];
		int[] real;
		double alpha = 0.9;
		
		NewmarkMethod nm = new NewmarkMethod();
		Matrix M = nm.readMatrix("c:/optimization/hw8/mass.txt");
		Matrix C = nm.readMatrix("c:/optimization/hw8/damping.txt");
		Matrix eqe = nm.readEarthquake("c:/optimization/hw8/elcentro.txt");
		
		list.transform();
		for(i=0 ; i<list.population ; i++){
			real = list.realnumberList.get(i);
			for(j=0 ; j<list.variables ; j++)	
				k[i][j]=min+resolution*((double)real[j]);
		}		
		for(i=0 ; i<list.population ; i++) 
			value[i] = nm.process(alpha, M, C, eqe, k[i]);
		
		return value;
	}
	
	public double[] objectiveFunction(Chromosome list, boolean onlyBest){
		
		int i,j;
		double[] value = new double[list.population];
		double min=2.1*Math.pow(10, 6);
		double max=2.1*Math.pow(10, 10);
		int length=list.variableLengths[0];
		double r = Math.pow(2, length);
		double resolution = (max-min)/(r-1.0);
		double[][] k = new double[list.population][list.variables];
		int[] real;
		double alpha = 0.9;
		
		NewmarkMethod nm = new NewmarkMethod();
		Matrix M = nm.readMatrix("c:/optimization/hw8/mass.txt");
		Matrix C = nm.readMatrix("c:/optimization/hw8/damping.txt");
		Matrix eqe = nm.readEarthquake("c:/optimization/hw8/elcentro.txt");		
		
		list.transform();
		for(i=0 ; i<list.population ; i++){
			real = list.realnumberList.get(i);
			for(j=0 ; j<list.variables ; j++)	
				k[i][j]=min+resolution*((double)real[j]);
		}		
		for(i=0 ; i<list.population ; i++) 
			value[i] = nm.process(alpha, M, C, eqe, k[i]);
		
		if(onlyBest==false){
			for(i=0 ; i<list.population ; i++){
				System.out.printf("%6.5f\t",value[i]);
				for(j=0 ; j<list.variables ; j++) 
					System.out.printf("%f\t",k[i][j]);
				System.out.println();
			}
		}
		else if(onlyBest==true){
			double best = value[0];
			for(i=1 ; i<list.population ; i++) if(best>value[i]) best=value[i];
			for(i=1 ; i<list.population ; i++){
				if(best==value[i]){
					System.out.printf("%6.5f\t",value[i]);
					for(j=0 ; j<list.variables ; j++) 
						System.out.printf("%f\t",k[i][j]);
					System.out.println();
					break;
				}				
			}
		}
		
		return value;
	}
	
	public Chromosome createInitialPopulation(int population, 
										int variables, int length){
		int i;
		int[] variableLengths = new int[population];
		for(i=0 ; i<population ; i++) variableLengths[i] = length;
		Chromosome list = new Chromosome(population, variables*length);
		list.setVariables(variables, variableLengths);
		
		return list;
	}
	
	public int[] reproduction(Chromosome list){
		//roulette wheel selection
		int i,j;
		int size = list.population;
		double k=3.0;	//pressure parameter
		int[] selection = new int[2];
		double[] fx = objectiveFunction(list);
		double[] fitness = new double[size];
		double max=fx[0], min=fx[0], sum=0, tmpSum;
		double pressureFactor;
		double selectionPoint;
		
		for(i=0 ; i<size ; i++){
			if(fx[i]>max) max = fx[i];
			if(fx[i]<min) min = fx[i];
		}
		pressureFactor = (max-min)/(k-1.0);
		for(i=0 ; i<size ; i++){
			fitness[i] = max - fx[i] + pressureFactor;
			sum += fitness[i];
		}
		for(i=0 ; i<2 ; i++){
			selectionPoint = Math.random() * sum;
			tmpSum = 0;
			for(j=0 ; j<size ; j++){
				tmpSum += fitness[i];
				if(selectionPoint<tmpSum){
					selection[i] = j;
					break;
				}
			}
		}		
		
		return selection;		
	}
	
	public int[][] crossover(int size, int[] one, int[] two, double rate){
		int crossoverPoint;
		int[][] offspring = new int[2][size];
		
		if(Math.random()<rate){			
			crossoverPoint = (int)(Math.random()*size);			
			for(int i=0 ; i<size ; i++){
				if(i<crossoverPoint){
					offspring[0][i] = one[i];
					offspring[1][i] = two[i];
				}
				else{
					offspring[0][i] = two[i];
					offspring[1][i] = one[i];					
				}
			}
		}
		
		return offspring;
	}
	
	public void mutation(int size, int[][] offspring, double rate){
		//typically mutation rate is 0.001
		int i,j;
		for(i=0 ; i<2 ; i++){
			for(j =0 ; j<size ; j++){			
				if(Math.random() < rate){
					if(offspring[i][j]==0) offspring[i][j]=1;
					else offspring[i][j]=0;
				}
			}
		}
	}
	
	public int[][] geneticOperation(Chromosome list){		
		int[][] offspring;
		int[] selection;
		double crossoverRate = 0.9; 
		double mutationRate = 0.05;
		int[] one, two;
		
		selection = reproduction(list);
		one = list.chromosomeList.get(selection[0]);
		two = list.chromosomeList.get(selection[1]);
		offspring = crossover(list.chromosomeSize, one, two, crossoverRate);
		mutation(list.chromosomeSize, offspring, mutationRate);
		
		return offspring;
	}
	
	public double checkConvergence(double newFx, double oldFx){
		//Convergence Criteria: e1		
		return Math.abs((newFx - oldFx)/oldFx);
	}
	
	public void process(Chromosome list, double E1){
		int i,j,k;
		int size = list.population;
		int generation;
		int goTimes = 5;
		int[] rank = new int[size];	//rank[0] is wosrt
		int tmprank;
		int[][][] offspring = new int[goTimes][2][list.chromosomeSize];
		double[] fx;
		double oldFx, newFx;
		double e1=E1+1.0;	//convergence criteria
				
		fx = objectiveFunction(list);
		for(newFx=fx[0], i=1 ; i<size ; i++) if(newFx>fx[i]) newFx = fx[i];
		
		System.out.println("initial:");
		objectiveFunction(list,false);
		System.out.println("iterate:");
		for(generation=0 ; generation<200 ; generation++){
			oldFx = newFx;
			for(j=0 ; j<goTimes ; j++) offspring[j] = geneticOperation(list);
			
			for(j=0 ; j<size ; j++) rank[j] = j;
			for(j=1 ; j<size ; j++){
				for(k=0 ; k<j ; k++){
					if(fx[rank[k]]<fx[rank[j]]){
						tmprank = rank[k];
						rank[k] = rank[j];
						rank[j] = tmprank;
					}
				}
			}

			for(j=0 ; j<goTimes ; j++){
				list.chromosomeList.set(rank[j*2], offspring[j][0]);
				list.chromosomeList.set(rank[j*2+1], offspring[j][1]);
			}
			
			fx = objectiveFunction(list);
			for(newFx=fx[0], i=1 ; i<size ; i++) if(newFx>fx[i]) newFx = fx[i];
			
			e1 = checkConvergence(newFx, oldFx);
			
			objectiveFunction(list,true);
		}		
		System.out.println("result:");
		objectiveFunction(list,false);
	}
	
	public static void main(String[] args) {
		
		GeneticAlgorithm ga = new GeneticAlgorithm();
		Chromosome list = ga.createInitialPopulation(20, 9, 30);
		ga.process(list, 0.0001);

	}

}
