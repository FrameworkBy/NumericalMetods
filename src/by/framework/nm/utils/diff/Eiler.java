package by.framework.nm.utils.diff;

import static java.lang.System.out;

public class Eiler {

	private DiffFunction[] funcs;
	private int SIZE;
	private double t;
	private double T;

	private double Eidop = 10E-4;
	private double dTmax = 10E-6;

	public Eiler(double t0, double T, DiffFunction... funcs) {
		this.funcs = funcs;
		this.SIZE = funcs.length;
		this.t = t0;
		this.T = T;
	}

	public void solve(double... invars) throws Exception {

		double[] xvars = new double[SIZE];
		double[] Tik = new double[SIZE];
		double[] tmpvars = new double[SIZE];
		
		int c = 0;

		do { 
			++c;
			
			for(int i = 0; i < SIZE; i++) {
				funcs[i].setT(t);
			}
			for(int i = 0; i < SIZE; i++) {
				xvars[i] = funcs[i].calculate(invars);
			}		

			for(int i = 0; i < SIZE; i++) {
				Tik[i] = Eidop/Math.abs(xvars[i] + Eidop/dTmax);
			}
			double Tk = Tik[0];
			for(int i = 1; i < SIZE; i++) {
				if(Tk>Tik[i]) Tk = Tik[i];
			}

			for(int i = 0; i < SIZE; i++) {
				tmpvars[i] = xvars[i];
				xvars[i] = invars[i] + Tk*xvars[i];
				invars[i] = tmpvars[i];
			}

			t += Tk;

			System.out.print(c + " : " + t + " : ");
			out.print("[ ");
			for(double r : xvars) {
				out.printf("%f ", r);
			}
			out.println("]");
			//out.println(); out.println();

		} while (t < T);
	}

}
