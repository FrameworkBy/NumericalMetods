package by.framework.math;

public class IntegralTrapeze {

	protected double a;
	protected double b;
	protected double h;
	protected Function func;
	protected int steps = 1;
	protected double Eps = 1E-6;
	protected int iter = 0;

	public IntegralTrapeze(double start, double finish, Function func) {
		this.a = start;
		this.b = finish;
		this.func = func;
	}

	public double solve() throws Exception {
		
		double R = Double.MAX_VALUE;
		steps = 1;
		iter = 0;
		h = (b-a)/(double)steps;
		double I = (h/2.)*(func.calculate(a) + func.calculate(b));
		double I2;
		double x;
		double summ;
		do {
			++iter;
			steps *= 2;
			h /= 2.;
			summ = 0;
			x = a + h;
			for(int i = 1; i < steps; i++) {
				summ += func.calculate(x);
				x += h;
			}
			summ = (h/2.)*(func.calculate(a) + 2*summ + func.calculate(b));
			I2 = summ;
			R = Math.abs(I2 - I);
			I = I2;
		} while(R>=3*Eps);
		return I;
	}
	
	public int getIter() {
		return iter;
	}
	
	public int getSteps() {
		return steps;
	}

}
