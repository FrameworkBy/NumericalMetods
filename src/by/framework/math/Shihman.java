package by.framework.math;

import static java.lang.System.out;

public class Shihman {

	private DiffFunction[] funcs;
	private double t;
	private final double T;
	private final int SIZE;

	private final double dTmin = 1E-15;
	private final double Edop = 1E-8;

	public Shihman(double t, double T, DiffFunction... funcs) {
		this.funcs = funcs;
		this.SIZE = funcs.length;
		this.t = t;
		this.T = T;
	}

	public void solve(double... invars) {

		double[] Yprpr = new double[SIZE]; // Y k-2
		double[] Ypr = new double[SIZE]; // Y k-1
		double[] Yc = new double[SIZE]; // Y k
		double[] Yn = new double[SIZE]; // Y k+1

		Function[] functions = new Function[SIZE];

		double maxEk;

		for(int i = 0; i < SIZE; i++) {
			Yprpr[i] = invars[i];
			Ypr[i] = invars[i];
			Yc[i] = invars[i];
			Yn[i] = invars[i];
		}

		double a0, a1, b0;

		double dTn = dTmin;
		double dTc = dTmin;
		double dTp = dTmin;
		double dTpp = dTmin;

		double tnext = 0;

		int iter = 0;
		boolean flag = true;

		while(iter < 2) {
			tnext = t + dTc;
			a0 = 1;
			a1 = 0;
			b0 = dTc;

			for(int i = 0; i < SIZE; i++) {
				funcs[i].setT(tnext);
			}
			for(int i = 0; i < SIZE; i++) {
				functions[i] = new ShFunction(Ypr[i], Yc[i], funcs[i], i, a0, a1, b0);
			}

			try {
				Newton newton = new Newton(functions);
				Yn = newton.solve(Yn);
			} catch (Exception e) {
				flag = false;
				break;
			}
			
			maxEk = Double.MIN_VALUE;
			for(int i = 0; i < SIZE; i++) {
				maxEk = Math.max(maxEk, 
						Math.abs((dTc/(dTc-dTp)) * (Yn[i] - Yc[i] - (dTc/dTp)*(Yc[i]-Ypr[i]))));
			}

			if(maxEk>Edop) {
				dTc = dTc/2.;
				tnext = t;
				for(int i = 0; i < SIZE; i++) {
					Yn[i] = Yc[i];
				}
				continue;
			}

			/*
			System.out.print(iter + " : " + tnext + " : ");
			out.print("[ ");
			for(double r : Yn) {
				out.printf("%f ", r);
			}
			out.println("]");
			//*/

			if(maxEk<=(Edop/4.)) dTn = 2.*dTc;
			else dTn = dTc;

			for(int i = 0; i < SIZE; i++) {
				Yprpr[i] = Ypr[i];
				Ypr[i] = Yc[i];
				Yc[i] = Yn[i];
			}

			dTpp = dTp;
			dTp = dTc;
			dTc = dTn;
			t = tnext;

			++iter;
		}
		
		double Rk, kpp, kp, kc, kn; 
		if(flag)
		while(t < T) {
			
			if (t + dTc > T)
				dTc = T - t;
			
			tnext = t + dTc;

			a0 = ((dTc+dTp) * (dTc + dTp)) / (dTp * (2*dTc+dTp));
			a1 = -(dTc*dTc) / (dTp * (2*dTc+dTp));
			b0 = (dTc * (dTc+dTp)) / (2*dTc+dTp);

			for(int i = 0; i < SIZE; i++) {
				funcs[i].setT(tnext);
			}
			for(int i = 0; i < SIZE; i++) {
				functions[i] = new ShFunction(Ypr[i], Yc[i], funcs[i], i, a0, a1, b0);
			}

			try {
				Newton newton = new Newton(functions);
				Yn = newton.solve(Yn);
			} catch (Exception e) {
				break;
			}

			Rk = (dTc*dTc*dTc * (dTc + dTp)) / (2*dTc+dTp);
			kpp = dTpp*(dTp+dTpp)*(dTc+dTp+dTpp);
			kp = dTp*dTpp*(dTc+dTp);
			kc = dTc*dTp*(dTp+dTpp);
			kn = dTc*(dTc+dTp)*(dTc+dTp+dTpp);

			maxEk = Double.MIN_VALUE;
			for(int i = 0; i < SIZE; i++) {
				maxEk = Math.max(maxEk, 
						Math.abs(Rk*(Yn[i]/kn - Yc[i]/kc + Ypr[i]/kp - Yprpr[i]/kpp)));
			}

			if(maxEk>Edop) {
				dTc = dTc/2.;
				tnext = t;
				for(int i = 0; i < SIZE; i++) {
					Yn[i] = Yc[i];
				}
				continue;
			}

			if(maxEk<=(Edop/4.)) dTn = 2.*dTc;
			else dTn = dTc;

			for(int i = 0; i < SIZE; i++) {
				Yprpr[i] = Ypr[i];
				Ypr[i] = Yc[i];
				Yc[i] = Yn[i];
			}

			dTpp = dTp;
			dTp = dTc;
			dTc = dTn;
			t = tnext;

			++iter;
		}
		
		System.out.print(iter + " : " + tnext + " : ");
		out.print("[ ");
		for(double r : Yn) {
			out.printf("%f ", r);
		}
		out.println("]");

	}

	private class ShFunction extends Function {

		final double yp;
		final double yc;
		final Function f;
		final int j;
		final double a0;
		final double a1;
		final double b0;

		public ShFunction(double Ypri, double Yci, Function f, int pos, double a0, double a1, double b0) {
			this.yp = Ypri;
			this.yc = Yci;
			this.f = f;
			this.j = pos;
			this.a0 = a0;
			this.a1 = a1;
			this.b0 = b0;
		}

		@Override
		public double calculate(double... vars) throws Exception {
			return vars[j] - a1*yp - a0*yc - b0*f.calculate(vars);
		}

	}

}
