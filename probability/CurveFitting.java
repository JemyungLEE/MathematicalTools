package probability;

import java.io.File;
import java.util.Scanner;

import probability.data.*;
import probability.distribution.*;
import probability.function.FirstOrder;
import probability.function.SecondOrder;
import matrix.*;


public class CurveFitting {
	
	public CurveFitting(){}
	
	public double[] getMean(TwoDimensionData data){
		
		//[0]: x, [1]: y
		double dataNumber;
		double[] average = new double[2];

		average[0] = 0;
		average[1] = 0;
		dataNumber = data.getDataNumber();

		for(int i=0 ; i<dataNumber ; i++){
			average[0] += data.getValueX(i);
			average[1] += data.getValueY(i);
		}
		average[0] /= dataNumber;
		average[1] /= dataNumber;
		
		return average;
	}
	
	public double[] getVariance(TwoDimensionData data){
		
		//[0]: x, [1]: y
		double dataNumber;
		double[] average = getMean(data);
		double[] variance = new double[2];
		
		variance[0] = 0;
		variance[1] = 0;
		dataNumber = data.getDataNumber();

		for(int i=0 ; i<dataNumber ; i++){
			variance[0] += Math.pow(data.getValueX(i) - average[0], 2);
			variance[1] += Math.pow(data.getValueY(i) - average[1], 2);
		}
		variance[0] /= (dataNumber - 1);
		variance[1] /= (dataNumber - 1);
		
		return variance;
	}
	
	public double[] getVariance(TwoDimensionData data, double[] average){
		
		//[0]: x, [1]: y
		double dataNumber;
		double[] variance = new double[2];
		
		variance[0] = 0;
		variance[1] = 0;
		dataNumber = data.getDataNumber();

		for(int i=0 ; i<dataNumber ; i++){
			variance[0] += Math.pow(data.getValueX(i) - average[0], 2);
			variance[1] += Math.pow(data.getValueY(i) - average[1], 2);
		}
		variance[0] /= (dataNumber - 1);
		variance[1] /= (dataNumber - 1);
		
		return variance;
	}
	
	public void LeastSquaresFit(TwoDimensionData data, FirstOrder function){
		
		int i;
		double dataSize = (double)data.getDataNumber();
		double x, y;
		double meanx, meany;
		double sumx, sumy, sumx2, sumxy;
		double st, sr, r2;
		double a0, a1;
		
		//calculate sums
		sumx = 0;
		sumy = 0;
		sumx2 = 0;
		sumxy = 0;
		for(i=0 ; i<dataSize ; i++){
			x = data.getValueX(i);
			y = data.getValueY(i);
			sumx += x;
			sumy += y;
			sumx2 += x*x;
			sumxy += x*y;
		}		
		
		//calculate mean
		meanx = sumx / dataSize;
		meany = sumy / dataSize;
		
		//calculate a1
		a1 = (dataSize*sumxy-sumx*sumy)/(dataSize*sumx2-sumx*sumx);	
		
		//calculate a0
		a0 = meany - a1 * meanx;
			
		//calculate coefficient of determination, r2
		st = 0;
		sr = 0;
		for(i=0 ; i< dataSize ; i++){
			x = data.getValueX(i);
			y = data.getValueY(i);
			st += Math.pow((y-meany),2);
			sr += Math.pow((y-a1*x-a0),2);
		}
		r2 = (st-sr)/st;
		
		//set function
		function.setVariables(a0, a1);
		function.setCoefficientOfDetermination(r2);
	}
	
	public void LeastSquaresFit(TwoDimensionData data, SecondOrder function){
		
		int i;
		double dataSize = (double)data.getDataNumber();
		double x, y;
		double meanx, meany;
		double sumx, sumx2, sumx3, sumx4, sumy, sumxy, sumx2y;
		double[][] matrixA;
		double[] matrixB, matrixC;
		double st, sr, r2;
		double a0, a1, a2;
		MakeInverseMatrix inverse = new MakeInverseMatrix();
		
		//calculate sums
		sumx = 0;
		sumx2 = 0;
		sumx3 = 0;
		sumx4 = 0;
		sumy = 0;
		sumxy = 0;
		sumx2y = 0;
		for(i=0 ; i<dataSize ; i++){
			x = data.getValueX(i);
			y = data.getValueY(i);
			sumx += x;
			sumx2 += x*x;
			sumx3 += Math.pow(x, 3);
			sumx4 += Math.pow(x, 4);
			sumy += y;
			sumxy += x*y;
			sumx2y += x*x*y;
		}	
		
		//calculate mean
		meanx = sumx / dataSize;
		meany = sumy / dataSize;
		
		//calculate a0, a1, a2
		matrixA = new double[3][3];
		matrixB = new double[3];
		matrixC = new double[3];
		matrixA[0][0] = dataSize;
		matrixA[0][1] = sumx;
		matrixA[0][2] = sumx2;
		matrixA[1][0] = sumx;
		matrixA[1][1] = sumx2;
		matrixA[1][2] = sumx3;
		matrixA[2][0] = sumx2;
		matrixA[2][1] = sumx3;
		matrixA[2][2] = sumx4;
		matrixB[0] = sumy;
		matrixB[1] = sumxy;
		matrixB[2] = sumx2y;		
			
		inverse.GaussJordanMethod(3, matrixA, matrixB, matrixC);
		a0 = matrixC[0];
		a1 = matrixC[1];
		a2 = matrixC[2];
		
		//calculate coefficient of determination, r2
		st = 0;
		sr = 0;
		for(i=0 ; i< dataSize ; i++){
			x = data.getValueX(i);
			y = data.getValueY(i);
			st += Math.pow((y-meany),2);
			sr += Math.pow((y-a0-a1*x-a2*x*x),2);
		}
		r2 = (st-sr)/st;
		
		//set function
		function.setVariables(a0, a1, a2);
		function.setCoefficientOfDetermination(r2);
	}
	
	public void LeastSquaresFit(TwoDimensionData data, Normal dist){
		
		int i;
		int dataSize = data.getDataNumber();
		double mean, variance, std;
		double[] tmpVariables;
		double verification;
		SecondOrder tmpFunction = new SecondOrder();
		TwoDimensionData tmpData = new TwoDimensionData();
		
		//data transform 
		tmpData.setDataNumber(dataSize);
		for(i=0 ; i<dataSize ; i++){
			tmpData.addValueX(data.getValueX(i));
			tmpData.addValueY(Math.log(data.getValueY(i)));
		}

		//least squares fitting
		this.LeastSquaresFit(tmpData, tmpFunction);		
		tmpVariables = tmpFunction.getVariables();
		
		//calculate mean, variance, and standard deviation
		variance = -0.5 / tmpVariables[2];
		std = Math.sqrt(variance);
		mean = tmpVariables[1] * variance; 
		
		//verification
		verification = -0.5*mean*mean/variance
						-Math.log(std*Math.sqrt(2*Math.PI));
		if(Math.abs(tmpVariables[0]-verification)<Math.abs(0.03*verification))	
			System.out.println("verified.");
		else System.err.println("error: not verified");
		
		//set distribution
		dist.setMean(mean);
		dist.setSTD(std);
		dist.setVariance(variance);
	}
		
	public static void main(String[] args) {
		
		String source = "c:/test/5.txt";
		TwoDimensionData data = new TwoDimensionData(source);
		//FirstOrder function = new FirstOrder();
		//SecondOrder function = new SecondOrder();
		Normal function = new Normal();
		CurveFitting cf = new CurveFitting();
		double[] variables = new double[3];
		double mean, std;
		
		cf.LeastSquaresFit(data, function);
		//variables = function.getVariables();
		mean = function.getMean();
		std = function.getSTD();
		
		System.out.println("mean: "+mean);
		System.out.println("std: "+std);
		/*
		System.out.println("a0: "+variables[0]);
		System.out.println("a1: "+variables[1]);
		System.out.println("a2: "+variables[2]);
		System.out.println("functoin: "
				+variables[2]+" x^2 + "+variables[1]+" x + "+variables[0]);
		System.out.println("r2: "+function.getCoefficientOfDetermination());
		*/
	}

}
