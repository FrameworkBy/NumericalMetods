package by.framework.math;

import static java.lang.System.out;

public class NEiler extends Eiler {

	protected final double dTmin = 1E-15;

	public NEiler(double t0, double T, DiffFunction... funcs) {
		super(t0, T, funcs);
	}

	@Override
	public void solve(double... invars) {
		double[] Yprev = new double[SIZE]; // Y k-1
		double[] Ycur = new double[SIZE]; // Y k
		double[] Ynext = new double[SIZE]; // Y k+1
		
		Function[] functions = new Function[SIZE];

		for(int i = 0; i < SIZE; i++) {
			Yprev[i] = invars[i];
			Ycur[i] = invars[i];
			Ynext[i] = invars[i];
		}

		double dTk = dTmin;
		double dTkprev = dTmin, dTknext = dTmin;
		double tnext;
		double maxEk;

		int c = 0;

		do {

			if (t + dTk > T)
				dTk = T - t;

			tnext = t + dTk; 

			for(int i = 0; i < SIZE; i++) {
				funcs[i].setT(tnext);
			}

			for(int i = 0; i < SIZE; i++) {
				final double yk = Ycur[i];
				final double dt = dTk;
				final Function f = funcs[i];
				final int j = i;
				functions[i] = new Function() {
					public double calculate(double... vars) throws Exception {
						return vars[j] - yk - dt*f.calculate(vars);
					}
				};
			}
			try {
				Newton newton = new Newton(functions);
				Ynext = newton.solve(Ynext);
			} catch (Exception e) {
				break;
			}

			maxEk = Double.MIN_VALUE;
			for(int i = 0; i < SIZE; i++) {
				maxEk = Math.max(maxEk, 
						Math.abs((dTk/(dTk-dTkprev)) 
								* (Ynext[i] - Ycur[i] - (dTk/dTkprev)*(Ycur[i]-Yprev[i]))));
			}

			if(maxEk>Edop) {
				dTk = dTk/2.;
				tnext = t;
				for(int i = 0; i < SIZE; i++) {
					Ynext[i] = Ycur[i];
				}
				continue;
			}

			if(maxEk<=(Edop/4.)) dTknext = 2.*dTk;
			else dTknext = dTk;

			for(int i = 0; i < SIZE; i++) {
				Yprev[i] = Ycur[i];
				Ycur[i] = Ynext[i];
			}

			dTkprev = dTk;
			dTk = dTknext;
			t = tnext;
			
			++c;

		} while(t < T);
		System.out.print(c + " : " + tnext + " : ");
		out.print("[ ");
		for(double r : Ynext) {
			out.printf("%f ", r);
		}
		out.println("]");
	}

}
