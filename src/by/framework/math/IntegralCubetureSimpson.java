package by.framework.math;

public class IntegralCubetureSimpson {

	protected double a;
	protected double b;
	protected double c;
	protected double d;
	protected double hX;
	protected double hY;
	protected Function func;
	protected Function gfunc;
	protected int stepsX = 1;
	protected int stepsY = 1;
	protected double Eps = 1E-12;
	protected int iter = 0;
	protected int stepsMax = 1000;

	public IntegralCubetureSimpson(double start1, double finish1, 
			double start2, double finish2,
			Function func) {
		this.a = start1;
		this.b = finish1;
		this.c = start2;
		this.d = finish2;
		this.func = func;
		this.gfunc = new GFunction(start1, finish1, start2, finish2, func);
	}

	public double solve() throws Exception {
		double R = Double.MAX_VALUE;
		double I = 0, I2;
		double summ;
		do {
			++iter;
			stepsX *= 2;
			stepsY *= 2;
			hX = (b-a)/(double)(stepsX*2);
			hY = (d-c)/(double)(stepsY*2);
			summ = 0;
			Cache cache = new Cache(a, c, hX, hY, stepsX, stepsY, gfunc);
			for(int i = 0; i < stepsX; i++) {
				for(int j = 0; j < stepsY; j++) {
					summ += cache.get(2*i, 2*j);
					summ += 4*cache.get(2*i+1, 2*j);
					summ += cache.get(2*i+2, 2*j);
					summ += 4*cache.get(2*i, 2*j+1);
					summ += 16*cache.get(2*i+1, 2*j+1);
					summ += 4*cache.get(2*i+2, 2*j+1);
					summ += cache.get(2*i, 2*j+2);
					summ += 4*cache.get(2*i+1, 2*j+2);
					summ += cache.get(2*i+2, 2*j+2);
				}
			}
			summ = ((hX*hY)/9.)*summ;
			I2 = summ;
			R = Math.abs(I2 - I);
			I = I2;
		} while(R>=15*Eps);
		return I;
	}

	public int getIter() {
		return iter;
	}

	public int getSteps() {
		return stepsX;
	}
	
	private class Cache {
		Double[] cacheX;
		Double[] cacheY;
		
		private double stepX;
		private double stepY;
		
		private double startX;
		private double startY;
		
		private Double[][] cache;
		
		private Function func;
		
		public Cache(double startX, double startY, 
				double stepX, double stepY, 
				int n, int m,
				Function func) {
			cacheX = new Double[2*n+1];
			cacheY = new Double[2*m+1];
			cache = new Double[2*n+1][2*n+1];
			this.stepX = stepX;
			this.stepY = stepY;
			this.startX = startX;
			this.startY = startY;
			this.func = func;
		}
		
		public double get(int i, int j) throws Exception {	
			
			boolean t = true;
			
			if(t)
			return func.calculate(startX + i*stepX, startY + j*stepY);
			
			Double ans = cache[i][j];
			if(ans == null) {
				Double X = cacheX[i];
				if(X == null) {
					X = startX + i*stepX;
					cacheX[i] = X;
				}
				Double Y = cacheY[j];
				if(Y == null) {
					Y = startY + j*stepY;
					cacheY[j] = Y;
				}
				ans = func.calculate(X,Y);
				cache[i][j] = ans;
			}
			return ans;
		}
	}
	
	private class GFunction extends Function {

		private double a;
		private double b;
		private double c;
		private double d;
		private Function func;

		public GFunction(double a, double b, double c, double d, Function func) {
			this.a = a;
			this.b = b; 
			this.c = c;
			this.d = d;
			this.func = func;
		}

		@Override
		public double calculate(double... vars) throws Exception {
			if(vars[0] < a || vars[0] > b || vars[1] < c || vars[1] > d) return 0;
			return func.calculate(vars);
		}

	}

}
