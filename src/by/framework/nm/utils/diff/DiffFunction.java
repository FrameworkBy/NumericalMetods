package by.framework.nm.utils.diff;

import by.framework.nm.utils.sone.Function;

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
