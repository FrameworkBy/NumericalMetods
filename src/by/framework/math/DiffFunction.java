package by.framework.math;


public abstract class DiffFunction extends Function {
	
	protected double w;
	protected double t = 0;
	
	public DiffFunction(double w) {
		this.w = w;
	}
	
	public void setT(double T) {
		this.t = T;
	}
}
