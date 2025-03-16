package gammaspike.distribution;


/**
 * Provides a means to calculate the mean number of stubs on a branch, under a Poisson distribution
 */
public interface StubExpectation {
	/*
	 * The mean number of stubs along a lineage between two heights (parentHeight > height)
	 */
	public double getMeanStubNumber(double height, double parentHeight);
}
