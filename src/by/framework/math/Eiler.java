package by.framework.math;

import static java.lang.System.out;

public class Eiler {

	protected DiffFunction[] funcs;
	protected final int SIZE;
	protected double t;
	protected final double T;

	protected final double Edop = 1E-6;
	protected final double dTmax = 1E-6;

	public Eiler(double t0, double T, DiffFunction... funcs) {
		this.funcs = funcs;
		this.SIZE = funcs.length;
		this.t = t0;
		this.T = T;
	}

	public void solve(double... invars) throws Exception {

		double[] dU = new double[SIZE];
		double[] U = new double[SIZE];
		for(int i = 0; i < SIZE; i++) {
			U[i] = invars[i];
		}
		double dTk;

		int c = 0;
		do { 
			++c;
			dTk = dTmax;
			for(int i = 0; i < SIZE; i++) {
				funcs[i].setT(t);
				dU[i] = funcs[i].calculate(U);
				dTk = Math.min(dTk, Edop/(Math.abs(dU[i])));
			}
			if (t + dTk > T)
				dTk = T - t;
			
			for(int i = 0; i < SIZE; i++) {
				U[i] += dTk * dU[i];
			}

			t += dTk;

			/*
			System.out.print(c + " : " + t + " : ");
			out.print("[ ");
			for(double r : dU) {
				out.printf("%f ", r);
			}
			out.println("]");
			//*/

		} while (t < T);

		System.out.print(c + " : " + t + " : ");
		out.print("[ ");
		for(double r : U) {
			out.printf("%f ", r);
		}
		out.println("]");
	}

}
