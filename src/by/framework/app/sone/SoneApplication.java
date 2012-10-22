package by.framework.app.sone;
import static java.lang.System.out;
import by.framework.nm.utils.sone.Function;
import by.framework.nm.utils.sone.Newton;

public class SoneApplication {
	
	public static void main(String[] args) {
		try {
			Newton newton = new Newton(new Fun1(), new Fun2());
			double[] res = newton.solve(0.5, 1.);
			out.print("[ ");
			for(double r : res) {
				out.printf("%f ", r);
			}
			out.println("]");
			out.println(); out.println();
		} catch (Exception e) {
			out.println("Can't solve!");
			e.printStackTrace();
		}
	}
	
	private static class Fun2 extends Function {

		@Override
		public double calculate(double... vars) throws Exception {
			return vars[0]*vars[0] - vars[1]*vars[1] - 1.;
		}
		
	}
	
	private static class Fun1 extends Function {

		@Override
		public double calculate(double... vars) throws Exception {
			return 2.*vars[0]*vars[0] + 2.*vars[1]*vars[1] + 1.2 
					- Math.exp(Math.sqrt(vars[0]*vars[0] + vars[1]*vars[1]));
		}
		
	}
}
