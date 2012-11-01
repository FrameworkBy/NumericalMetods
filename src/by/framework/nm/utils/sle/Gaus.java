package by.framework.nm.utils.sle;

import java.util.ArrayList;

public class Gaus {
	private ArrayList<ArrayList<Double>> arrayA;
	private ArrayList<Double> arrayB;
	
	private ArrayList<ArrayList<Double>> inputArrayA;
	private ArrayList<Double> inputArrayB;
	
	private int SIZE;
	private boolean REZULT_FOUND = true;
	
	public enum Errors {DETERMINANT, INPUT_SIZE};
	
	private class Replace {
		private int src;
		private int dest;
		public Replace(int src, int dest) {
			this.src = src;
			this.dest = dest;
		}
		public int getSrc() {
			return src;
		}
		public int getDest() {
			return dest;
		}
	}
	
	private ArrayList<Replace> replaces = new ArrayList<Replace>();
	
	private Errors error;
	
	public Gaus(Double[][] A, Double[] B) {
		arrayB = new ArrayList<Double>();
		for(int i = 0; i < B.length; i++) {
			arrayB.add(B[i]);
		}
		
		arrayA = new ArrayList<ArrayList<Double>>();
		for(int i = 0; i < A.length; i++) {
			ArrayList<Double> list = new ArrayList<Double>();
			for(int j = 0; j < A[i].length; j++) {
				list.add(A[i][j]);
			}
			arrayA.add(list);
		}
		SIZE = arrayA.size();
		saveInput();
	}
	
	public Gaus(ArrayList<ArrayList<Double>> A, ArrayList<Double> B) {
		arrayA = new ArrayList<ArrayList<Double>>();
		for(ArrayList<Double> list : A) {
			ArrayList<Double> newList = new ArrayList<Double>();
			for(Double a : list) {
				newList.add(a);
			}
			arrayA.add(newList);
		}
		
		arrayB = new ArrayList<Double>();
		for(Double b : B) {
			arrayB.add(b);
		}
		SIZE = arrayA.size();
		saveInput();
	}
	
	private void saveInput() {
		this.inputArrayA = new ArrayList<ArrayList<Double>>(SIZE);
		this.inputArrayB = new ArrayList<Double>(SIZE);
		for(ArrayList<Double> list : arrayA) {
			ArrayList<Double> newList = new ArrayList<Double>();
			for(Double a : list) {
				newList.add(a.doubleValue());
			}
			inputArrayA.add(newList);
		}
		for(Double b : arrayB) {
			inputArrayB.add(b);
		}
	}
	
	public void solve() {
		REZULT_FOUND = testSize();
		if(REZULT_FOUND) {
			firstRun();
			if(REZULT_FOUND)
				secondRun();
		}
	}
	
	public ArrayList<Double> getArrayResult() {
		if(!REZULT_FOUND) {
			return null;
		} else {
			return arrayB;
		}
	}
	
	public double[] getResult() {
		double[] res = new double[SIZE];
		for(int i = 0; i < SIZE; i++) {
			res[i] = arrayB.get(i);
		}
		return res;
	}
	
	public Errors getError() {
		return error;
	}
	
	private void sortArrays(int start) {
		int maxI = start, maxJ = start;
		int k = start + 1;
		Double max = Math.abs(arrayA.get(maxI).get(maxJ));
		for(int i = start; i < SIZE; i++) {
			for(int j = (k>start)?k--:k; j < SIZE; j++) {
				if(Double.compare(Math.abs(arrayA.get(i).get(j)), max) >= 0) {
					max = Math.abs(arrayA.get(i).get(j));
					maxI = i;
					maxJ = j;
				}
			}
		}
		if(start != maxI) {
			ArrayList<Double> tmpA = arrayA.get(start);
			arrayA.set(start, arrayA.get(maxI));
			arrayA.set(maxI, tmpA);
			
			Double tmpB = arrayB.get(start);
			arrayB.set(start, arrayB.get(maxI));
			arrayB.set(maxI, tmpB);
		}
		
		if(start != maxJ) {
			for(ArrayList<Double> list : arrayA) {
				Double tmpA2 = list.get(start);
				list.set(start, list.get(maxJ));
				list.set(maxJ, tmpA2);
			}
			replaces.add(new Replace(start, maxJ));
		}
	}
	
	private void firstRun() {
		for(int i = 0; i < SIZE; i++) {
			sortArrays(i);
			Double main = arrayA.get(i).get(i);
			
			if(Math.abs((main)) < 10e-15 ) {
				REZULT_FOUND = false;
				error = Errors.DETERMINANT;
				return;
			}
			for(int j = i; j < SIZE; j++) {
				arrayA.get(i).set(j, arrayA.get(i).get(j)/main);
			}
			arrayB.set(i, arrayB.get(i)/main);
			
			for(int k = i+1; k < SIZE; k++) {
				Double z = arrayA.get(k).get(i);
				arrayA.get(k).set(i, 0.);
				for(int n = i+1; n < SIZE; n++) {
					arrayA.get(k).set(n, arrayA.get(k).get(n)-arrayA.get(i).get(n)*z);
				}
				arrayB.set(k, arrayB.get(k)-arrayB.get(i)*z);
			}
		}
	}
	
	private void secondRun() {
		for(int i = SIZE-1; i > 0; i--) {
			for(int j = i - 1; j >= 0; j--) {
				arrayB.set(j, arrayB.get(j) - arrayB.get(i)*arrayA.get(j).get(i));
			}
		}
		
		for(int i = replaces.size()-1; i >= 0; i--) {
			Replace rep = replaces.get(i);
			Double tmp = arrayB.get(rep.getSrc());
			arrayB.set(rep.getSrc(), arrayB.get(rep.getDest()));
			arrayB.set(rep.getDest(), tmp);
		}
	}
	
	private boolean testSize() {
		for(ArrayList<Double> list : arrayA) {
			if(list.size() != SIZE) {
				error = Errors.INPUT_SIZE;
				return false;
			}
		}
		if(arrayB.size() != SIZE) {
			error = Errors.INPUT_SIZE;
			return false;
		}
		return true;
	}
	
	public Double testResult() {
		ArrayList<Double> AX = new ArrayList<Double>();
		for(int i = 0; i < SIZE; i++)
			AX.add(0.);
		for(int i = 0; i < SIZE; i++) {
			for(int j = 0; j < SIZE; j++) {
				AX.set(i, AX.get(i) + inputArrayA.get(i).get(j)*arrayB.get(j));
			}
		}
		for(int i = 0; i < AX.size(); i++) {
			AX.set(i, AX.get(i)-inputArrayB.get(i));
		}
		Double res = new Double(0.);
		for(Double a : AX) {
			res += a*a;
		}
		return Math.sqrt(res);
	}
}
