package by.framework.math;
import static java.lang.System.out;

public class Newton {
	
	private Function[] functions;
	private double[] arrayB;
	private double[][] arrayA;
	
	private double[] variables;
	
	private int SIZE = 0;
	
	private static final double eps = 10E-6;
	private static final double dx = 10E-8;
	private static final int MAX_ITER = 1000;
	
	public Newton(Function... functions) throws Exception {
		this.SIZE = functions.length;
		this.functions = functions; 
	}
	
	public double[] solve(double... inX) throws Exception {
		if(inX.length != SIZE) throw new Exception();
		
		double[] x = new double[SIZE];
		for(int i = 0; i < SIZE; i++) {
			x[i] = inX[i];
		}
		
		this.arrayA = new double[SIZE][SIZE];
		this.arrayB = new double[SIZE];
		
		variables = new double[SIZE];
		
		int iter = 0;
		double Rf, Rx;
		do {
			iter++;
			if(iter > MAX_ITER) throw new Exception();
			for(int i = 0; i < SIZE; i++) {
				variables[i] = x[i];
			}
			
			createJ();
			Gaus gaus = new Gaus(arrayA, arrayB);
			gaus.solve();
			if(gaus.getError() != null) {
				out.println(gaus.getError());
				throw new Exception();
			}
			double[] res = gaus.getResult();
			
			Rf = 0;
			Rx = 0;
			for(int j = 0; j < SIZE; j++) {
				x[j] = variables[j] + res[j];
			}
			
			for(int j = 0; j < SIZE; j++) {
				double b = functions[j].calculate(x);
				Rf += b*b;
				Rx += res[j]*res[j];
			}
			Rf = Math.sqrt(Rf);
			Rx = Math.sqrt(Rx);
			
		} while(Rx >= eps && Rf >= eps);
		
		return x;
	}
	
	private void createJ() throws Exception {
		for(int i = 0; i < SIZE; i++) {
			double[] tmp_vars1 = new double[SIZE];
			double[] tmp_vars2 = new double[SIZE];
			for(int vi = 0; vi < SIZE; vi++) {
				tmp_vars1[vi] = variables[vi];
				tmp_vars2[vi] = variables[vi];
			}
			tmp_vars1[i] = tmp_vars1[i] + dx;
			tmp_vars2[i] = tmp_vars2[i] - dx;
			
			for(int j = 0; j < SIZE; j++) {
				double df = calcDiff(j, tmp_vars1, tmp_vars2);
				arrayA[j][i] = df;
			}
			arrayB[i] = -functions[i].calculate(variables);
		}
	}
	
	private double calcDiff(int fun, double[] vars1, double[] vars2) 
			throws Exception {
		double f1 = functions[fun].calculate(vars1);
		double f2 = functions[fun].calculate(vars2);
		double t = f1 - f2;
		double d = t/(2*dx);
		return d;
	}
}
