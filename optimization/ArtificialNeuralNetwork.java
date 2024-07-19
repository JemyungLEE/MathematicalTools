package optimization;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Date;
import java.util.Random;
import java.util.Scanner;

public class ArtificialNeuralNetwork {

	/**
	 *  Subject: Artificial Neural Network simulation tool
	 *  Developer: Jemyung Lee
	 *  Developed Data: 2008.8.12
	 *  Last Modified Data: 2010.11.28 
	 *  Department: Seoul Nat. Univ. depart. of Rural Systems Engineering
	 *  Description:- input data maker, output data maker included
	 *  			- training function, simulation function included
	 *  			- input, output & hidden node number can be set by user
	 *  			- simulation function operate both simulation & training	
	 */

	
	private int dataNumber;
	private int nodeValueNumber;
	private int inputNodeSort, outputNodeSort; //outputNodeSort=outputNodeNumber
	private int inputNodeNumber, outputNodeNumber, hiddenNodeNumber;
	private int[] inputNodeNumberBySort;
	private double[][] inputNodeValue;	//node value number * input nodes
	private double[][] outputNodeValue;	//node value number * output nodes
	private double[][] observedValue;	//node value number * output nodes
	private ArrayList<Integer> dataID; 
	private ArrayList<ArrayList<Double>> inputLayerData;
	private ArrayList<ArrayList<Double>> outputLayerData;
	private double[] inputLayerMax, inputLayerMin;	//[inputNodeSort]
	private double[] outputLayerMax, outputLayerMin;	//[outputNodeSort]
	
	private double totalError;
	private double[] inputHiddenWeight;
	private double[] hiddenOutputWeight;
	
	
	public ArtificialNeuralNetwork(){
	}
		
	public ArtificialNeuralNetwork(String[] argv, 
										int[] nodeNumber, int nhid, int nout){
		process(argv, nodeNumber, nhid, nout);
	}
	
	public void process(String[] argv, int[] inputNode, 
											int hiddenNode, int outputNode) {
	
		this.setNodeValues(argv[0], inputNode, hiddenNode, outputNode);
		System.out.println("node value setting ok.");		
		this.training();
		System.out.println("training ok.");		
		this.setNodeValues(argv[1], inputNode, hiddenNode, outputNode);
		System.out.println("node value setting ok.");
		this.simulation();
		System.out.println("simulation ok.");
		this.makeOutputFile(argv[2]);
		System.out.println("output file printting ok.");
		
		System.out.println("process complete");
	}      
	
	
	private double generateRandomNumber(Random rand){	       
	    		
		return rand.nextDouble() - 0.5;
	}
	
	private double SigmoidFunction(double value){
	       
		return 1.0/(1.0+Math.exp(-1.0*value));
	}
	
	public void training(){
        // ETA: momentum rate
        // alpha: learning rate
        // AMAXE: Available maximum error
        double ETA = 0.95;
        double ALPHA = 0.7;
        double AMAXE = 0.001;
        
        this.training(ETA, ALPHA, AMAXE);
	}
	
	public void training(double ETA, double ALPHA, double AMAXE){
        // ETA: momentum rate
        // alpha: learning rate
        // AMAXE: Available maximum error
		
		int i, j, k;
	    
	    int ndata = this.nodeValueNumber;	//number of data set
	    int nin = this.inputNodeNumber;		//number of input node
	    int nhid = this.hiddenNodeNumber;	//number of hide node
	    int nout = this.outputNodeNumber;	//number of output node
	    
	    double[][] input = this.inputNodeValue.clone();
	    double[][] output = new double[ndata][nout];
	    double[][] obs = this.observedValue.clone();    
	    
	    
    	double[] W1 = new double[(nin+1)*nhid];
    	double[] W2 = new double[(nhid+1)*nout];
    	double[] olddw1 = new double[(nin+1)*nhid];
    	double[] olddw2 = new double[(nhid+1)*nout];
    	double[] dw1 = new double[(nin+1)*nhid];
    	double[] dw2 = new double[(nhid+1)*nout];
        
    	//Initialization of weight W1(input-hidden layer), W2(hidden-output layer)
        Random r = new Random();
        Date t = new Date();
        r.setSeed(t.getTime());
        for(i=0; i<=nin; i++){
        	for(j=0; j<nhid; j++){
        		W1[j+i*nhid] = generateRandomNumber(r);
        	}
        }
        for(i=0; i<=nhid; i++){
        	for(j=0; j<nout; j++){
        		W2[j+i*nout] = generateRandomNumber(r);
        	}
        }
        //Initialization of delta weight olddw1, olddw2            
        for(i=0; i<=nin; i++){
        	for(j=0; j<nhid; j++){
        		olddw1[j+i*nhid] = 0.0;
        	}
        }
        for(i=0; i<=nhid; i++){
        	for(j=0; j<nout; j++){
        		olddw2[j+i*nout] = 0.0;
        	}
        }
        //Choose an input-output pair
        int iter;
        double[] X = new double[nin+1];
        double[] Y = new double[nout];
        double[] H = new double[nhid+1];
        double[] EP = new double[nout];
        double anet, bnet;
        double[] delta1 = new double[nhid];
        double[] delta2 = new double[nout];
        double sumdw, sum_diff, total_err;
        
        iter = 0;
        do{
        	total_err = 0.0;  
        	for(i=0; i<ndata; i++){
        		X[0] = 0.0;
        		for(j=0; j<nin; j++){
        			X[j+1] = input[i][j];
        		}
        		for(j=0; j<nout; j++){
        			Y[j] = obs[i][j];
        		}

        		//Calculation of NET(P,J)
        		H[0] = 1.0;
        		for(j=0; j<nhid; j++){
        			anet = 0.0;
        			for(k=0; k<=nin ; k++){
        				anet += X[k] * W1[k+j*(nin+1)];
        			}
        			H[j+1] = SigmoidFunction(anet);
        		}
        		//Calculation of NET(P,K)                     
        		for(j=0; j<nout; j++){
        			bnet = 0.0;
        			for(k=0; k<=nhid; k++){
        				bnet += H[k] * W2[k+j*(nhid+1)];
        			}
        			output[i][j] = SigmoidFunction(bnet);
        		}
        		//Calculation of adjustment of weight (hidden-output)
        		for(j=0; j<nout; j++){
        			delta2[j] = output[i][j] * (1 - output[i][j])
        								* (Y[j] - output[i][j]);
        		}
        		for(j=0; j<nhid; j++){
        			sumdw = 0.0;
        			for(k=0; k<nout; k++){
        				sumdw += delta2[k] * W2[j+1+k*(nhid+1)];
        			}
        			delta1[j] = H[j+1] * (1 - H[j+1]) * sumdw;
        		}
        		//Weight adjustment (hidden - output)
        		for(j=0; j<nout; j++){
        			for(k=0; k<=nhid; k++){
        				dw2[k+j*(nhid+1)] = ETA * delta2[j] * H[k] + ALPHA 
        									* olddw2[k+j*(nhid+1)];
        			}
        		}
        		for(j=0; j<nhid; j++){
        			for(k=0; k<=nin; k++){
        				dw1[k+j*(nin+1)] = ETA * delta1[j] * X[k] + ALPHA
        									* olddw1[k+j*(nin+1)];
        			}
        		}
        		//Weight Adjustment (input-hidden)
        		for(j=0; j<nhid; j++){
        			for(k=0; k<=nin; k++){
        				olddw1[k+j*(nin+1)] = dw1[k+j*(nin+1)];
        				W1[k+j*(nin+1)] += dw1[k+j*(nin+1)];
        			}
        		}
        		for(j=0; j<nout; j++){
        			for(k=0; k<=nhid; k++){
        				olddw2[k+j*(nhid+1)] = dw2[k+j*(nhid+1)];
        				W2[k+j*(nhid+1)]+= dw2[k+j*(nhid+1)];
        			}
        		}
        		//Error calculation                     
        		sum_diff = 0.0;
        		for(j=0; j<nout; j++){
        			EP[j] = Math.abs(Y[j] - output[i][j]);
        			sum_diff += Math.pow(EP[j], 2);
        		}
        		total_err += 0.5 * sum_diff;
        	}
        	iter++;

        //System.out.println("iteration:"+iter+"\t error:"+total_err);
        //}while(total_err > AMAXE && iter < 10000);
        }while(iter < 100000);

        //store the results
        this.totalError = total_err;
        this.inputHiddenWeight = W1.clone();
        this.hiddenOutputWeight = W2.clone();

	}
	
	public void simulation(){
        // ETA: momentum rate
        // alpha: learning rate
        // AMAXE: Available maximum error
        double ETA = 0.95;
        double ALPHA = 0.7;
        double AMAXE = 0.001;
        
        this.simulation(ETA, ALPHA, AMAXE);
	}
	
	public void simulation(double ETA, double ALPHA, double AMAXE){
        // ETA: momentum rate
        // alpha: learning rate
        // AMAXE: Available maximum error
		
		int i, j, k;
	    
	    int ndata = this.nodeValueNumber;	//number of data set
	    int nin = this.inputNodeNumber;		//number of input node
	    int nhid = this.hiddenNodeNumber;	//number of hide node
	    int nout = this.outputNodeNumber;	//number of output node
	    
	    double[][] input = this.inputNodeValue.clone();
	    double[][] output = new double[ndata][nout];
	    double[][] obs = this.observedValue.clone();    
	    
	    
    	double[] W1 = this.inputHiddenWeight.clone();
    	double[] W2 = this.hiddenOutputWeight.clone();
    	double[] olddw1 = new double[(nin+1)*nhid];
    	double[] olddw2 = new double[(nhid+1)*nout];
    	double[] dw1 = new double[(nin+1)*nhid];
    	double[] dw2 = new double[(nhid+1)*nout];
			    
	
		// Initialization of delta weight olddw1, olddw2            
        for(i=0; i<=nin; i++){
        	for(j=0; j<nhid; j++){
        		olddw1[j+i*nhid] = 0.0;
        	}
        }
        for(i=0; i<=nhid; i++){
        	for(j=0; j<nout; j++){
        		olddw2[j+i*nout] = 0.0;
        	}
        }
		
        //Simulation
		double[] X=new double[nin+1];
		double[] Y=new double[nout];
		double[] H=new double[nhid+1];
		double anet, bnet;
      
        double[] delta1 = new double[nhid];
        double[] delta2 = new double[nout];
        double sumdw;
	
		for(i=0; i<ndata; i++){ 
			//System.out.println("data no.: "+(i+1));
			X[0] = 0.0;
			for(j=0; j<nin; j++){
				X[j+1] = input[i][j];
			}
			for(j=0; j<nout; j++){
				Y[j] = obs[i][j];
			}
           
			//Calculation of NET(P,J)
			H[0] = 1.0;
			for(j=0; j<nhid; j++){
				anet = 0.0;
				for(k=0; k<=nin; k++){
					anet += X[k] * W1[k+j*(nin+1)];
				}
				H[j+1] = SigmoidFunction(anet);
			}
          
			//Calculation of NET(P,K)
			for(j=0; j<nout; j++){
				bnet = 0.0;
				for(k=0; k<=nhid; k++){
					bnet += H[k] * W2[k+j*(nhid+1)];
				}
				output[i][j] = SigmoidFunction(bnet);
			}
			
			// Calculation of adjustment of weight (hidden-output)
    		for(j=0; j<nout; j++){
    			delta2[j] = output[i][j] * (1 - output[i][j])
    								* (Y[j] - output[i][j]);
    		}
    		for(j=0; j<nhid; j++){
    			sumdw = 0.0;
    			for(k=0; k<nout; k++){
    				sumdw += delta2[k] * W2[j+1+k*(nhid+1)];
    			}
    			delta1[j] = H[j+1] * (1 - H[j+1]) * sumdw;
    		}
    		//Weight adjustment (hidden - output)
    		for(j=0; j<nout; j++){
    			for(k=0; k<=nhid; k++){
    				dw2[k+j*(nhid+1)] = ETA * delta2[j] * H[k] + ALPHA 
    									* olddw2[k+j*(nhid+1)];
    			}
    		}
    		for(j=0; j<nhid; j++){
    			for(k=0; k<=nin; k++){
    				dw1[k+j*(nin+1)] = ETA * delta1[j] * X[k] + ALPHA
    									* olddw1[k+j*(nin+1)];
    			}
    		}
    		//Weight Adjustment (input-hidden)
    		for(j=0; j<nhid; j++){
    			for(k=0; k<=nin; k++){
    				olddw1[k+j*(nin+1)] = dw1[k+j*(nin+1)];
    				W1[k+j*(nin+1)] += dw1[k+j*(nin+1)];
    			}
    		}
    		for(j=0; j<nout; j++){
    			for(k=0; k<=nhid; k++){
    				olddw2[k+j*(nhid+1)] = dw2[k+j*(nhid+1)];
    				W2[k+j*(nhid+1)]+= dw2[k+j*(nhid+1)];
    			}
    		}
		}
		
		//set output node data, simulation result
		this.outputNodeValue = output.clone();		
	}
	

	public void readData(String source){

		//input layer's name should contain 'input'.
		//output layer's name should contain 'output'.
		//data ID should be integer.
		
		int i;
		String layerName;
		
		//initiate
		this.dataID = new ArrayList<Integer>();
		this.inputLayerData = new ArrayList<ArrayList<Double>>();
		this.outputLayerData = new ArrayList<ArrayList<Double>>();
		
		try{
			File file = new File(source);
			Scanner scan = new Scanner(file);

			//check the number of input and output layers
			scan.next();	//read name of data ID
			while(!scan.hasNextInt()){	//read name of layers				
				layerName = scan.next().toLowerCase();
				if(layerName.contains("input")){	//add input layer
					this.inputLayerData.add(new ArrayList<Double>());					
				}else if(layerName.contains("output")){	//add output layer
					this.outputLayerData.add(new ArrayList<Double>());
				}else{ 
					System.err.println("index of input file is wrong.");
					break;
				}
			}
			
			//set input,output layer sort numbers
			this.inputNodeSort = this.inputLayerData.size();
			this.outputNodeSort = this.outputLayerData.size();
						
			//read data
			while(scan.hasNext()){
				this.dataID.add(scan.nextInt());	//read data ID
				for(i=0 ; i<this.inputNodeSort ; i++) //read input layer's data
					this.inputLayerData.get(i).add(scan.nextDouble());
				for(i=0 ; i<this.outputNodeSort ; i++)//read output layer's data
					this.outputLayerData.get(i).add(scan.nextDouble());
			}
			
			//set number of data
			this.dataNumber = dataID.size();
			
			scan.close();
			
		} catch (Exception e) {
			System.err.println("error_readData");		
		}
	}

	public void setNodeNumbers(int[] inputNodeBySort, 
			int hiddenNode, int outputNode){
	
		int i;
		int inputNode = 0;
		int sort = inputNodeBySort.length;

		//initiate
		this.inputNodeNumberBySort = new int[this.inputNodeSort];
		
		//set input node number according to each node type
		if(sort==this.inputNodeSort){		
			for(i=0 ; i<sort ; i++){
				this.inputNodeNumberBySort[i] = inputNodeBySort[i];
				inputNode += inputNodeBySort[i];
			}
		}
		else if(sort==1){
			for(i=0 ; i<this.inputNodeSort ; i++){
				this.inputNodeNumberBySort[i] = inputNodeBySort[0];
				inputNode += inputNodeBySort[0];
			}
		}
		else System.err.println("input node number by sort does not match.");
		
		//set input, hidden, output node numbers
		this.inputNodeNumber = inputNode;
		this.hiddenNodeNumber = hiddenNode;
		this.outputNodeNumber = outputNode;		
		
	}
	
	public void setNodeValues(String inputFile, int[] inputNodeBySort, 
										int hiddenNode, int outputNode){
		this.readData(inputFile);
		this.setNodeNumbers(inputNodeBySort, hiddenNode, outputNode);
		this.setNodeValues();
	}
	
	public void setNodeValues(){
		//temporary variables for brief, fast calculation 
		int i,j,k;		
		int size = this.dataNumber;		
		int inputSort = this.inputNodeSort;
		int outputSort = this.outputNodeSort;
		int[] inputNodeBySort = this.inputNodeNumberBySort.clone(); 
		int maxInputNode;	//node number of largest group(sort) in input data
		int tmpSize;	//temporary size = size - maxInputNode
		int tmpIndex;
		double tmpSpan;	//span = max - min		
		double[] inputMax = new double[inputSort];
		double[] inputMin = new double[inputSort];
		double[] outputMax = new double[outputSort];
		double[] outputMin = new double[outputSort];
		double[][] inputData, outputData;	//data sort * data number		
				
		//initiate
		inputData = new double[inputSort][size];
		outputData = new double[outputSort][size];
		
		//copy data to temporary variables
		for(i=0 ; i<inputSort ; i++)
			for(j=0 ; j<size ; j++)
				inputData[i][j] = this.inputLayerData.get(i).get(j);
		for(i=0 ; i<outputSort ; i++)
			for(j=0 ; j<size ; j++)
				outputData[i][j] = this.outputLayerData.get(i).get(j);				
		
		//check maximum, minimum of data
		for(i=0 ; i<inputSort ; i++){
			inputMax[i] = inputData[i][0];
			inputMin[i] = inputData[i][0];	
			for(j=1 ; j<size ; j++){
				if(inputMax[i]<inputData[i][j]) inputMax[i]=inputData[i][j];
				if(inputMin[i]>inputData[i][j]) inputMin[i]=inputData[i][j];
			}
		}
		for(i=0 ; i<outputSort ; i++){
			outputMax[i] = outputData[i][0];
			outputMin[i] = outputData[i][0];
			for(j=1 ; j<size ; j++){
				if(outputMax[i]<outputData[i][j]) outputMax[i]=outputData[i][j];
				if(outputMin[i]>outputData[i][j]) outputMin[i]=outputData[i][j];
			}
		}
		
		//rearrange max and min
		for(i=0 ; i<inputSort ; i++){
			inputMax[i] *= 1.2;
			inputMin[i] *= 0.8;
		}
		for(i=0 ; i<outputSort ; i++){
			outputMax[i] *= 1.2;
			outputMin[i] *= 0.8;
		}
		
		//normalize data values
		for(i=0 ; i<inputSort ; i++){
			tmpSpan = inputMax[i] - inputMin[i]; 
			for(j=0 ; j<size ; j++)
				inputData[i][j] = (inputData[i][j]-inputMin[i])/tmpSpan;			
		}
		for(i=0 ; i<outputSort ; i++){
			tmpSpan = outputMax[i] - outputMin[i];
			for(j=0 ; j<size ; j++)				
				outputData[i][j] = (outputData[i][j]-outputMin[i])/tmpSpan;
		}	
		
		//check largest group's ,sort's , node number in input layer. 
		maxInputNode = inputNodeBySort[0];
		for(i=1 ; i<inputSort ; i++){
			if(maxInputNode < inputNodeBySort[i]) 
				maxInputNode = inputNodeBySort[i];
		}
		
		//set max, min
		this.inputLayerMax = inputMax.clone();
		this.inputLayerMin = inputMin.clone();
		this.outputLayerMax = outputMax.clone();
		this.outputLayerMin = outputMin.clone();
		
		//set input, output layer's node data		
		tmpSize = size - maxInputNode + 1;
		this.nodeValueNumber = tmpSize;
		this.inputNodeValue = new double[tmpSize][this.inputNodeNumber];
		this.observedValue = new double[tmpSize][this.outputNodeNumber];		
		for(i=0 ; i<tmpSize ; i++){
			tmpIndex = 0;
			for(j=0 ; j<inputSort ; j++){
				for(k=0 ; k<inputNodeBySort[j] ; k++){
					this.inputNodeValue[i][tmpIndex]
					                      = inputData[j][i+maxInputNode-1-k];
					tmpIndex++;
				}
			}
			for(j=0 ; j<outputSort ; j++){
				this.observedValue[i][j] = outputData[j][i+maxInputNode-1];
			}
		}	
	}		

	public void makeOutputFile(String outputfile){
		
		//temporary variables for brief expression
		int i,j;
		int nodeNumber = this.outputNodeNumber;
		int dataNumber = this.nodeValueNumber;
		double[] max = this.outputLayerMax.clone();
		double[] min = this.outputLayerMin.clone();		
		double[][] predict = this.outputNodeValue.clone();
		double[][] observed = this.observedValue.clone();
				
		//restore size of predict, observed data
		for(i=0 ; i<dataNumber ; i++){
			for(j=0 ; j<nodeNumber ; j++){
				predict[i][j] = (predict[i][j]*(max[j]-min[j]))+min[j] ;
				observed[i][j] = (observed[i][j]*(max[j]-min[j]))+min[j] ;
			}
		}

		//make output file
		try{
			File file = new File(outputfile);
			PrintWriter pw = new PrintWriter(file);
			
			//write condition
			for(i=0 ; i<this.inputNodeSort ; i++) 
				pw.println("N_in_"+(i+1)+": "+this.inputNodeNumberBySort[i]);
			pw.println("N_hidden: "+this.hiddenNodeNumber);
			pw.println("N_output: "+this.outputNodeNumber);
			pw.println("N_data: "+this.nodeValueNumber);
			pw.println();
			for(i=0 ; i<this.inputNodeSort ; i++){
				pw.print("Max_input_"+(i+1)+": "+this.inputLayerMax[i]+"\t");
				pw.println("Min_input_"+(i+1)+": "+this.inputLayerMin[i]);
			}
			for(i=0 ; i<this.outputNodeSort ; i++){
				pw.print("Max_output_"+(i+1)+": "+this.outputLayerMax[i]+"\t");
				pw.println("Min_output_"+(i+1)+": "+this.outputLayerMin[i]);
			}
			pw.println();
			
			//write predict,observed data
			pw.print("no.\t");
			for(i=0 ; i<nodeNumber ; i++){
				pw.print("predict_"+(i+1)+"\tobserved_"+(i+1)+"\t");
			}
			pw.println();
			
			for(i=0 ; i<dataNumber ; i++){
				pw.print(i+1+"\t");
				for(j=0 ; j<nodeNumber ; j++){
					pw.print(predict[i][j]+"\t"+observed[i][j]+"\t");
				}
				pw.println();
			}
			
			pw.close();
			System.out.println("result: "+outputfile);
		}catch(FileNotFoundException e) { 
			e.printStackTrace(); 
		}
		
	}
	
	public void readWeightFile(String weightFile){
	
		int i,j;
		int inputNode, hiddenNode, outputNode;
		
   		try{
			File file = new File(weightFile);
			Scanner scan = new Scanner(file);
			
			//read condition
			for(i=0 ; i<this.inputNodeSort ; i++){ 
				scan.next();
				this.inputNodeNumberBySort[i] = scan.nextInt();
			}
			scan.next();
			inputNode = scan.nextInt();
			scan.next();
			hiddenNode = scan.nextInt();
			scan.next();
			outputNode = scan.nextInt();
			scan.next();
			scan.next();	//skip data number
			scan.next();
			this.totalError = scan.nextDouble();
		
			//store the node numbers
			this.inputNodeNumber = inputNode;
			this.hiddenNodeNumber = hiddenNode;
			this.outputNodeNumber = outputNode;
			
			//initiate
			this.inputHiddenWeight = new double[(inputNode+1)*hiddenNode];
			this.hiddenOutputWeight = new double[(hiddenNode+1)*outputNode];
			
			//Read weighting factor from learning result
			for(i=0; i<=inputNode; i++)
				for(j=0; j<hiddenNode; j++)
					this.inputHiddenWeight[j+i*hiddenNode] = scan.nextDouble();
			for(i=0; i<=hiddenNode; i++)
				for(j=0; j<outputNode; j++)
					this.hiddenOutputWeight[j+i*outputNode] = scan.nextDouble();
			scan.close();
		}catch(FileNotFoundException e) { 
			e.printStackTrace(); 
		}
	}
	
	public void makeWeightFile(String weightFile){

		//temporary variables for brief expression
		int i,j;
		int inputNode = this.inputNodeNumber;
		int hiddenNode = this.hiddenNodeNumber;
		int outputNode = this.outputNodeNumber;
		
   		try{
			File file = new File(weightFile);
			PrintWriter pw = new PrintWriter(file);
			
			//write condition
			for(i=0 ; i<this.inputNodeSort ; i++) 
				pw.println("N_in_"+(i+1)+": "+this.inputNodeNumberBySort[i]);
			pw.println("N_in: "+this.inputNodeNumber);
			pw.println("N_hidden: "+this.hiddenNodeNumber);
			pw.println("N_output: "+this.outputNodeNumber);
			pw.println("N_data: "+this.nodeValueNumber);
			pw.println("Total_error: "+this.totalError);
			pw.println();

			//write weight value
			for(i=0; i<=inputNode; i++){
				for(j=0; j<hiddenNode; j++){
					pw.printf("%10.6f",this.inputHiddenWeight[j+i*hiddenNode]);
					pw.println();
				}
			}
			for(i=0; i<=hiddenNode; i++){
				for(j=0; j<outputNode; j++){
					pw.printf("%10.6f",this.hiddenOutputWeight[j+i*outputNode]);
					pw.println();
				}
			}
			
			//print weight file name
			System.out.println("weight file: "+weightFile);
			pw.close();
		}catch(FileNotFoundException e) {
			e.printStackTrace();
		}  
	}

	public void printNodeValues(){
		
	}
	
	public void printNodeValues(String nodeValueFile){
		
		int i,j;
		
  		try{
			File file = new File(nodeValueFile);
			PrintWriter pw = new PrintWriter(file);
			
			//write condition
			for(i=0 ; i<this.inputNodeSort ; i++) 
				pw.println("N_in_"+(i+1)+": "+this.inputNodeNumberBySort[i]);
			pw.println("N_hidden: "+this.hiddenNodeNumber);
			pw.println("N_output: "+this.outputNodeNumber);
			pw.println("N_data: "+this.nodeValueNumber);
			pw.println();
			for(i=0 ; i<this.inputNodeSort ; i++){
				pw.print("Max_input_"+(i+1)+": "+this.inputLayerMax[i]+"\t");
				pw.println("Min_input_"+(i+1)+": "+this.inputLayerMin[i]);
			}
			for(i=0 ; i<this.outputNodeSort ; i++){
				pw.print("Max_output_"+(i+1)+": "+this.outputLayerMax[i]+"\t");
				pw.println("Min_output_"+(i+1)+": "+this.outputLayerMin[i]);
			}
			pw.println();
			
			//write node values			
			for(i=0 ; i<this.nodeValueNumber ; i++){
				pw.print((i+1)+"\t");
				for(j=0 ; j<this.inputNodeNumber ; j++)
					pw.print(this.inputNodeValue[i][j]+"\t");
				for(j=0 ; j<this.outputNodeNumber ; j++)
					pw.print(this.observedValue[i][j]+"\t");
				pw.println();
			}
			
			//print node value file name
			System.out.println("node_value file: "+nodeValueFile);
			pw.close();
		}catch(FileNotFoundException e) {
			e.printStackTrace();
		}  
	}
	
	public static void main(String[] args){
		
		String path = "d:/ann/test/";
		String lemInput = path+"rawdata_lem_test2.txt";
		String simInput = path+"rawdata_sim_test2.txt";
		String weight = path+"weight2.txt";
		String output = path+"result2.txt";
		String lemNode = path+"lem_node_20_19_1.txt";
		String simNode = path+"sim_node_20_19_1.txt";
		
		int[] nin = {10,10};
		int nhid = 19;
		int nout = 1;		
		
		ArtificialNeuralNetwork ann = new ArtificialNeuralNetwork();
		
		
		//ann.setNodeValues(lemInput, nin, nhid, nout);
		//System.out.println("node value setting ok.");
		
		//ann.printNodeValues(lemNode);
		//System.out.println("node value printing ok.");
		
		//ann.training();
		//System.out.println("training ok.");
		
		//ann.makeWeightFile(weight);
		//System.out.println("weight file making ok.");
		
		
		ann.setNodeValues(simInput, nin, nhid, nout);
		System.out.println("node value setting ok.");
		
		//ann.printNodeValues(simNode);
		//System.out.println("node value printing ok.");
		
		ann.readWeightFile(weight);
		System.out.println("weight file reading ok.");
		
		ann.simulation();
		System.out.println("simulation ok.");
		
		ann.makeOutputFile(output);
		System.out.println("output file printting ok.");
		
		System.out.println("process complete");
	
	}
}