package by.framework.nm.utils.sle;

import java.util.ArrayList;

public class Gaus2 {
	
	private boolean REZULT_FOUND = true;
	private Errors error;
	
	public enum Errors {DETERMINANT, INPUT_SIZE};
	
	private double[][] srcArrayA;
	private double[] srcArrayB;
	private int SIZE;
	
	private ArrayA arrayA;
	private ArrayB arrayB;
	
	private Double[][] inputArrayA;
	private Double[] inputArrayB;
	
	public Gaus2(Double[][] arrayA, Double[] arrayB) {
		if(testArraysSize(arrayA, arrayB)) {
			this.srcArrayA = new double[arrayA.length][arrayA.length];
			this.srcArrayB = new double[arrayB.length];
			
			this.inputArrayA = arrayA;
			this.inputArrayB = arrayB;
			
			for(int i = 0; i < arrayA.length; i++) {
				for(int j = 0; j < arrayA[i].length; j++) {
					this.srcArrayA[i][j] = arrayA[i][j];
				}
			}
			
			for(int i = 0; i < arrayB.length; i++) {
				this.srcArrayB[i] = arrayB[i];
			}
			this.SIZE = arrayB.length;
			makeReplaceArrays();
		}
	}
	
	public Gaus2(ArrayList<ArrayList<Double>> arrayA, ArrayList<Double> arrayB) {
		if(testArraysSize(arrayA, arrayB)) {
			this.srcArrayA = new double[arrayA.size()][arrayA.size()];
			this.srcArrayB = new double[arrayB.size()];
			
			this.inputArrayA = new Double[arrayA.size()][arrayA.size()];
			this.inputArrayB = new Double[arrayB.size()];
			
			for(int i = 0; i < arrayA.size(); i++) {
				for(int j = 0; j < arrayA.get(i).size(); j++) {
					this.srcArrayA[i][j] = arrayA.get(i).get(j);
					this.inputArrayA[i][j] = arrayA.get(i).get(j);
				}
			}
			
			for(int i = 0; i < arrayB.size(); i++) {
				this.srcArrayB[i] = arrayB.get(i);
				this.inputArrayB[i] = arrayB.get(i);
			}
			this.SIZE = arrayB.size();
			makeReplaceArrays();
		}
	}
	
	public void solve() {
		if(REZULT_FOUND) {
			firstRun();
			if(REZULT_FOUND) {
				secondRun();
			}
		}
	}
	
	public double[] getResult() {
		if(REZULT_FOUND) {
			double[] res = new double[SIZE];
			for(int i = 0; i < SIZE; i++) {
				res[i] = arrayB.getR(i);
			}
			return res;
		} else {
			return null;
		}
	}
	
	public Errors getError() {
		return error;
	}
	
	public double testResult() {
		if(!REZULT_FOUND) {
			return 0;
		}
		double[] AX = new double[SIZE];
		double[] X = getResult();
		for(int i = 0; i < SIZE; i++) {
			AX[i] = 0;
		}
		for(int i = 0; i < SIZE; i++) {
			for(int j = 0; j < SIZE; j++) {
				AX[i] = AX[i] + inputArrayA[i][j]*X[j];
			}
		}
		for(int i = 0; i < SIZE; i++) {
			AX[i] = AX[i] - inputArrayB[i];
		}
		double res = 0;
		for(int i = 0; i < SIZE; i++) {
			res += AX[i]*AX[i];
		}
		res = Math.sqrt(res);
		return res;
	}
	
	private boolean testArraysSize(Double[][] A, Double[] B) {
		int S = A.length;
		for(int i = 0; i < S; i++) {
			if(A[i].length != S) {
				REZULT_FOUND = false;
				error = Errors.INPUT_SIZE;
				return false;
			}
		}
		if(B.length != S) {
			REZULT_FOUND = false;
			error = Errors.INPUT_SIZE;
			return false;
		}
		return true;
	}
	
	private boolean testArraysSize(ArrayList<ArrayList<Double>> A, ArrayList<Double> B) {
		int S = A.size();
		for(int i = 0; i < S; i++) {
			if(A.get(i).size() != S) {
				REZULT_FOUND = false;
				error = Errors.INPUT_SIZE;
				return false;
			}
		}
		if(B.size() != S) {
			REZULT_FOUND = false;
			error = Errors.INPUT_SIZE;
			return false;
		}
		return true;
	}
	
	private void makeReplaceArrays() {
		arrayA = new ArrayA(srcArrayA);
		arrayB = new ArrayB(srcArrayB);
	}
	
	private void firstRun() {
		for(int i = 0; i < SIZE; i++) {
			sortArrays(i);
			double main = arrayA.get(i, i);
			if(Math.abs((main)) < 10e-15 ) {
				REZULT_FOUND = false;
				error = Errors.DETERMINANT;
				return;
			}
			for(int j = i; j < SIZE; j++) {
				arrayA.set(i, j, arrayA.get(i, j)/main);
			}
			arrayB.set(i, arrayB.get(i)/main);
			
			for(int k = i+1; k < SIZE; k++) {
				double z = arrayA.get(k, i);
				arrayA.set(k, i, 0);
				for(int n = i+1; n < SIZE; n++) {
					double valueA = arrayA.get(k, n) - arrayA.get(i, n)*z;
					arrayA.set(k, n, valueA);
				}
				double valueB = arrayB.get(k) - arrayB.get(i)*z;
				arrayB.set(k, valueB);
			}
		}
	}
	
	private void secondRun() {
		for(int i = SIZE-1; i > 0; i--) {
			for(int j = i - 1; j >= 0; j--) {
				double value = arrayB.get(j) - arrayB.get(i)*arrayA.get(j, i);
				arrayB.set(j, value);
			}
		}
		arrayB.normolize(arrayA.getReplaceJ());
	}
	
	private void sortArrays(int start) {
		int maxI = start, maxJ = start;
		int k = start + 1;
		double max = Math.abs(arrayA.get(maxI, maxJ));
		
		for(int i = start; i < SIZE; i++) {
			for(int j = (k>start)?k--:k; j < SIZE; j++) {
				double fc = Math.abs(arrayA.get(i, j));
				if(Double.compare(fc, max) >= 0) {
					max = Math.abs(arrayA.get(i, j));
					maxI = i;
					maxJ = j;
				}
			}
		}
		
		if(start != maxI) {
			arrayA.replaceI(maxI, start);
			arrayB.replaceI(maxI, start);
		}
		
		if(start != maxJ) {
			arrayA.replaceJ(maxJ, start);
		}
	}
	
	private class ArrayA {
		
		private double[][] arrayA;
		
		private int[] replaceI;
		private int[] replaceJ;
		
		public ArrayA(double[][] arrayA) {
			this.arrayA = arrayA;
			replaceI = new int[arrayA.length];
			replaceJ = new int[arrayA.length];
			for(int i = 0; i < arrayA.length; i++) {
				replaceI[i] = i;
				replaceJ[i] = i;
			}
		}
		
		public double get(int i, int j) {
			return arrayA[replaceI[i]][replaceJ[j]];
		}
		
		public void set(int i, int j, double value) {
			arrayA[replaceI[i]][replaceJ[j]] = value;
		}
		
		public void replaceI(int max, int dest) {
			int tmp = replaceI[dest];
			replaceI[dest] = replaceI[max];
			replaceI[max] = tmp;
		}
		
		public void replaceJ(int max, int dest) {	
			
			int tmp = replaceJ[dest];
			replaceJ[dest] = replaceJ[max];
			replaceJ[max] = tmp;
		}
		
		public int[] getReplaceJ() {
			return replaceJ;
		}
		
	}
	
	private class ArrayB {
		private double[] arrayB;
		
		private int[] replaceI;
		
		public ArrayB(double[] arrayB) {
			this.arrayB = arrayB;
			replaceI = new int[arrayB.length];
			for(int i = 0; i < arrayB.length; i++) {
				replaceI[i] = i;
			}
		}
		
		public double get(int i) {
			return arrayB[replaceI[i]];
		}
		
		public void set(int i, double value) {
			arrayB[replaceI[i]] = value;
		}
		
		public void replaceI(int max, int dest) {	
			int tmp = replaceI[dest];
			replaceI[dest] = replaceI[max];
			replaceI[max] = tmp;
		}
		
		public double getR(int i) {
			return arrayB[i];
		}
		
		public void normolize(int[] replaceJ) {
			
			double[] tmp = new double[arrayB.length];
			for(int i = 0; i < arrayB.length; i++) {
				tmp[i] = arrayB[replaceI[i]];
			}
			for(int i = 0; i < arrayB.length; i++) {
				arrayB[replaceJ[i]] = tmp[i];
			}
			
		}
	}

}
