package by.framework.math;

public class IntegralSimpson extends IntegralTrapeze {

	public IntegralSimpson(double start, double finish, Function func) {
		super(start, finish, func);
	}

	public double solve() throws Exception {
		double R = Double.MAX_VALUE;
		iter = 0;
		steps = 2;
		h = (b-a)/(double)steps;
		double I = (h/3.)*(func.calculate(a)+func.calculate(b) + 4.*func.calculate((a+b)/2.)); 
		double I2;
		double x;
		double summ2;
		double summ4;
		do {
			++iter;
			steps *= 2;
			h /= 2.;
			summ2 = 0;
			summ4 = 0;
			x = a + h;
			
			for(int i = 1; i <= steps/2; i++) {
				summ4 += func.calculate(x);
				x += h+h;
			}
			
			x = a + h + h;
			for(int i = 1; i < steps/2; i++) {
				summ2 += func.calculate(x);
				x += h+h;
			}

			I2 = (summ4*4. + summ2 * 2. + func.calculate(a) + func.calculate(b)) * (h/3.);
			R = Math.abs(I2 - I);
			I = I2;
		} while(R>=15*Eps);
		return I;
	}

}
