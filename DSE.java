import java.io.*;

public class DSE {

	public static double [] means = new double[]{ 269.568837,	141.985532, 99.508883, 79.006262, 66.454092, 58.649068, 52.804705, 48.960141, 45.73406, 43.455542, 41.386889, 40.046129, 38.695398, 37.750183, 36.840822, 36.41925, 35.715777, 35.350802, 34.767695, 34.493228, 34.033093, 33.944396, 33.676262, 33.871417, 33.654013, 33.737325, 33.825452, 34.021759, 34.223725, 34.479067, 34.694037, 35.157394}; 
	public static double [] stddevs = new double [] {0.676473946, 0.335842225, 0.208458629, 0.319806191, 0.313754681, 0.340024999, 0.340440891, 0.441634464, 0.424892928, 0.45393722, 0.470011702, 0.479161768, 0.496199557, 0.493329504, 0.524060111, 0.523248507, 0.564837145, 0.590188106, 0.584813646, 0.602546264, 0.597015913, 0.600979201, 0.59552246, 0.595876665, 0.587283577, 0.645952011, 0.666749578, 0.732030737, 0.836999403, 0.942281274, 1.046115194, 1.103043517}; 

	public static int maxCores = 16;
	public static double deadline = 300;


	public static double [] [] [] [] [] Utilization = new double [maxCores] [maxCores] [maxCores] [maxCores] [maxCores];
	public static double [] [] [] [] [] Risk = new double [maxCores] [maxCores] [maxCores] [maxCores] [maxCores];
	public static boolean [] [] [] [] [] Valid = new boolean [maxCores] [maxCores] [maxCores] [maxCores] [maxCores];
	
	
	
	public static void main(String[] args) throws IOException {
		
		File logFile=new File("Data.csv");

	    BufferedWriter writer = new BufferedWriter(new FileWriter(logFile));
	    
		
		
		
		//Test for 5 Phases
		
		int counter = 0;
		for (int a = 0; a <maxCores; a++) {
			for (int b = 0; b<maxCores;b++) {
				for (int c = 0; c < maxCores; c++) {
					for (int d = 0; d < maxCores;d++) {
						for (int e = 0; e < maxCores; e++) {
							Utilization [a] [b] [c] [d] [e] = (a * means [a] + b * means [b] + c * means [c] + d * means [d] + e * means [e])/(deadline * maxCores) * 100;
							Risk [a] [b] [c] [d] [e] = (1 - normp(ZScore(deadline, means[a] + means [b] + means [c] + means [d] + means [e], stddevs[a] + stddevs[b] + stddevs[c] + stddevs[d] + stddevs[e])))*100;
							Valid [a] [b] [c] [d] [e] = true;
							
							if (Risk [a] [b] [c] [d] [e] > 1 && Risk [a] [b] [c] [d] [e] < 99) {
								counter++;
								//System.out.println(counter + ": "  + "(" + a + "," + b  + "," + c  + "," + d + "," + e + ")" + ": Utilization  = " + Utilization [a] [b] [c] [d] [e] + "%" + " : Risk = " + Risk [a] [b] [c] [d] [e] + "%") ;
								writer.write (counter+","+Utilization [a] [b] [c] [d] [e] +","+Risk [a] [b] [c] [d] [e]+"\n");
							}
						}
					}
				}
			}
		}
		
		writer.close();
		System.out.println ("Valid Points: " +  calculateValidPoints () );

	}
	
public static void DSE () {
	
	int numberOfSteps = 0;
	//int minUniformPoint 
	
	for (int i = 0; i<maxCores;i++) {
		numberOfSteps++;
		
	}
	
	
	
}
	
		
public static int calculateValidPoints () {
	
	int count = 0;
	for (int a = 0; a <maxCores; a++) {
		for (int b = 0; b<maxCores;b++) {
			for (int c = 0; c < maxCores; c++) {
				for (int d = 0; d < maxCores;d++) {
					for (int e = 0; e < maxCores; e++) {
						if (Valid[a][b][c][d][e]) count++; 
					}
				}
			}
		}
	}
	return count;
}
	
public static double valueFromZScore (double ZScore, double mean, double standardDeviation) {
		
		return (ZScore * standardDeviation) + mean;
	}
	
	public static double ZScore (double value, double mean, double standardDeviation) {
		
		return (value - mean)/standardDeviation;
		
	}
	
	
	 public static double normp(double z) {

	      double zabs;
	      double p;
	      double expntl,pdf;

	      final double p0 = 220.2068679123761;
	      final double p1 = 221.2135961699311;
	      final double p2 = 112.0792914978709;
	      final double p3 = 33.91286607838300;
	      final double p4 = 6.373962203531650;
	      final double p5 = .7003830644436881;
	      final double p6 = .3526249659989109E-01;

	      final double q0 = 440.4137358247522;
	      final double q1 = 793.8265125199484;
	      final double q2 = 637.3336333788311;
	      final double q3 = 296.5642487796737;
	      final double q4 = 86.78073220294608;
	      final double q5 = 16.06417757920695;
	      final double q6 = 1.755667163182642;
	      final double q7 = .8838834764831844E-1;

	      final double cutoff = 7.071;
	      final double root2pi = 2.506628274631001;

	      zabs = Math.abs(z);

	//  |z| > 37

	      if (z > 37.0) {

	         p = 1.0;

	         return p;

	      }

	      if (z < -37.0) {

	         p = 0.0;

	         return p;

	      }

	//  |z| <= 37.

	      expntl = Math.exp(-.5*zabs*zabs);

	      pdf = expntl/root2pi;

	//  |z| < cutoff = 10/sqrt(2).

	      if (zabs < cutoff) {

	         p = expntl*((((((p6*zabs + p5)*zabs + p4)*zabs + p3)*zabs +
	             p2)*zabs + p1)*zabs + p0)/(((((((q7*zabs + q6)*zabs +
	             q5)*zabs + q4)*zabs + q3)*zabs + q2)*zabs + q1)*zabs +
	             q0);

	      } else {

	         p = pdf/(zabs + 1.0/(zabs + 2.0/(zabs + 3.0/(zabs + 4.0/
	             (zabs + 0.65)))));

	      }

	      if (z < 0.0) {

	         return p;

	      } else {

	         p = 1.0 - p;

	         return p;

	      }

	   }
	 
	 public static double xnormi(double p) {

	      double arg,t,t2,t3,xnum,xden,qinvp,x,pc;

	      final double c[] = {2.515517, 
	      .802853,
	      .010328};

	      final double d[] = {1.432788,
	      .189269,
	      .001308};

	      if (p <= .5) {

	         arg = -2.0*Math.log(p);
	         t = Math.sqrt(arg);
	         t2 = t*t;
	         t3 = t2*t;

	         xnum = c[0] + c[1]*t + c[2]*t2;
	         xden = 1.0 + d[0]*t + d[1]*t2 + d[2]*t3;
	         qinvp = t - xnum/xden;
	         x = -qinvp;

	         return x;

	      }

	      else {

	         pc = 1.0 - p;
	         arg = -2.0*Math.log(pc);
	         t = Math.sqrt(arg);
	         t2 = t*t;
	         t3 = t2*t;

	         xnum = c[0] + c[1]*t + c[2]*t2;
	         xden = 1.0 + d[0]*t + d[1]*t2 + d[2]*t3;
	         x = t - xnum/xden;

	         return x;

	      }

	   }

}
