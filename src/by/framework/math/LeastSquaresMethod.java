package by.framework.math;

import static java.lang.System.out;

import java.util.List;


public class LeastSquaresMethod {

	private int N;
	private int m;
	private List<Pair> pairs;
	private double q = 0;

	public LeastSquaresMethod(List<Pair> pairs, int N, int m) {
		this.pairs = pairs;
		this.N = N;
		this.m = m;
	}

	public double[] solve() throws Exception {
		double[] POWERX = new double[2*m+1];
		double[] XP = new double[N];
		for(int i = 0; i < N; i++) {
			XP[i] = 1;
		}
		for(int k = 0; k <= 2*m; k++) {
			POWERX[k] = 0;
			for(int i = 0; i < N; i++) {
				POWERX[k] += XP[i];
				XP[i] *= pairs.get(i).getX();
			}
		}

		double[][] A = new double[m+1][m+1];
		for(int i = 0; i < m+1; i++) {
			for(int j = 0; j < m+1; j++) {
				A[i][j] = POWERX[j + i];
			}
		}
		
		double[] B = new double[m+1];
		for(int i = 0; i < N; i++) {
			XP[i] = 1;
		}
		for(int i = 0; i < m+1; i++) {
			for(int j = 0; j < N; j++) {
				B[i] += XP[j] * pairs.get(j).getY();
				XP[j] *= pairs.get(j).getX();
			}
		}

		Gaus gaus = new Gaus(A, B);
		gaus.solve();
		if(gaus.getError() != null) {
			out.println(gaus.getError());
			throw new Exception();
		}
		double[] res = gaus.getResult();

		double summ = 0;
		for(int j = 0; j < N; j++) {
			double result = res[res.length - 1];
			double x = pairs.get(j).getX(), y = pairs.get(j).getY();
			for(int i = res.length - 2; i >= 0; i--) {
				result = result*x + res[i];
			}
			summ += (result-y)*(result-y);
		}
		this.q = Math.sqrt(summ/(N-m-1));

		return res;
	}

	public double getQ() {
		return this.q;
	}

}
