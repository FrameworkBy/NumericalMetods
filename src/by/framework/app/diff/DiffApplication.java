package by.framework.app.diff;

import by.framework.nm.utils.diff.DiffFunction;
import by.framework.nm.utils.diff.Eiler;

public class DiffApplication {
	
	public static void main(String[] args) throws Exception {
		double w = 25;
		Eiler eiler = new Eiler(0 + 10E-6, 1, new Fun1(w), new Fun2(w));
		eiler.solve(0, -0.412);
	}
	
}

class Fun1 extends DiffFunction {

	public Fun1(double w) {
		super(w);
	}

	@Override
	public double calculate(double... vars) throws Exception {
		return -vars[0]*vars[1]+Math.sin(t)/t;
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
		return -vars[1]*vars[1] + (a*t)/(1+t*t);
	}
}