package by.framework.app.diff;

import by.framework.nm.utils.diff.DiffFunction;
import by.framework.nm.utils.diff.Eiler;
import by.framework.nm.utils.diff.NEiler;
import by.framework.nm.utils.diff.Shihman;

public class DiffApplication {

	public static void main(String[] args) {
		double w = 25;
		DiffFunction f1 =  new Fun1(w);
		DiffFunction f2 =  new Fun2(w);
		double t = 1E-15;
		double T = 1;
		Eiler eiler = new Eiler(t, T, f1, f2);
		try {
			eiler.solve(0, -0.412);
		} catch(Exception ex) {}
		NEiler neiler = new NEiler(t, T, f1, f2);
		try {
			neiler.solve(0, -0.412);
		} catch(Exception ex) { }
		Shihman shihman = new Shihman(t, T, f1, f2);
		try {
			shihman.solve(0, -0.412);
		} catch(Exception ex) { }
	}

}

class Fun1 extends DiffFunction {

	public Fun1(double w) {
		super(w);
	}

	@Override
	public double calculate(double... vars) throws Exception {
		return -(vars[0]*vars[1])+Math.sin(t)/t;
	}
}

class Fun2 extends DiffFunction {

	private double a;

	public Fun2(double w) {
		super(w);
		a = 2.5 + w/40.;
	}

	@Override
	public double calculate(double... vars) throws Exception {
		return -(vars[1]*vars[1]) + (a*t)/(1+t*t);
	}
}