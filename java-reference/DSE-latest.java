import java.io.*;
import java.util.*;

public class DSE {

	public static int maxCores = 16;
	public static double [] means = new double[]{ 269.568837,	141.985532, 99.508883, 79.006262, 66.454092, 58.649068, 52.804705, 48.960141, 45.73406, 43.455542, 41.386889, 40.046129, 38.695398, 37.750183, 36.840822, 36.41925, 35.715777, 35.350802, 34.767695, 34.493228, 34.033093, 33.944396, 33.676262, 33.871417, 33.654013, 33.737325, 33.825452, 34.021759, 34.223725, 34.479067, 34.694037, 35.157394}; 
	public static double [] stddevs = new double [] {0.676473946, 0.335842225, 0.208458629, 0.319806191, 0.313754681, 0.340024999, 0.340440891, 0.441634464, 0.424892928, 0.45393722, 0.470011702, 0.479161768, 0.496199557, 0.493329504, 0.524060111, 0.523248507, 0.564837145, 0.590188106, 0.584813646, 0.602546264, 0.597015913, 0.600979201, 0.59552246, 0.595876665, 0.587283577, 0.645952011, 0.666749578, 0.732030737, 0.836999403, 0.942281274, 1.046115194, 1.103043517}; 

	public static double deadline = 216;
	public static double DMR  = 30;

	public static DesignPoint [][][][][] DesignPoints = new DesignPoint [maxCores] [maxCores] [maxCores] [maxCores] [maxCores];
	
	public static Stack <DesignPoint> DesignStack = new Stack  <DesignPoint> ();
	public static ArrayList <DesignPoint> RandomList = new ArrayList<>();
	public static int numberOfSteps = 0;
	
	
	public static double getRisk (DesignPoint Point) {
		
		if (!Point.operated) {
			numberOfSteps++;
			Point.operated = true;
		}
		return 100 * (1 - normp(ZScore(deadline,  means[Point.a] + means [Point.b] + means [Point.c] + means [Point.d] + means [Point.e], Math.sqrt(stddevs[Point.a]*stddevs[Point.a] + stddevs[Point.b]*stddevs[Point.b] + stddevs[Point.c]*stddevs[Point.c] + stddevs[Point.d]*stddevs[Point.d] + stddevs[Point.e]*stddevs[Point.e]))));
	}

	public static double getExecutionTime (DesignPoint Point) {
		/*if (!Point.operated) {
			numberOfSteps++;
			Point.operated = true;
		}*/
		return means[Point.a] + means [Point.b] + means [Point.c] + means [Point.d] + means [Point.e];
	}

	public static double getUtilization (DesignPoint Point) {
		/*if (!Point.operated) {
		numberOfSteps++;
		Point.operated = true;
		}*/
		return ((Point.a+1) * means [Point.a] + (Point.b+1) * means [Point.b] + (Point.c+1) * means [Point.c] + (Point.d+1) * means [Point.d] + (Point.e+1) * means [Point.e])/(deadline * maxCores) * 100;
	}

		
	
	
	public static void main(String[] args) throws IOException {
		
		double lowestUtilization = 100;
		DesignPoint optimalPoint = null;
			
		
		for (int a = 0; a <maxCores; a++) {
			for (int b = 0; b<maxCores;b++) {
				for (int c = 0; c < maxCores; c++) {
					for (int d = 0; d < maxCores;d++) {
						for (int e = 0; e < maxCores; e++) {
							DesignPoints [a] [b] [c] [d] [e] = new DesignPoint(a, b, c, d, e);
							RandomList.add(DesignPoints [a] [b] [c] [d] [e]);
							if (getRisk(DesignPoints [a] [b] [c] [d] [e]) < DMR && getUtilization(DesignPoints [a] [b] [c] [d] [e]) < lowestUtilization) {
								lowestUtilization = getUtilization(DesignPoints [a] [b] [c] [d] [e]);
								optimalPoint = DesignPoints [a] [b] [c] [d] [e];
							}

							DesignPoints [a] [b] [c] [d] [e].operated = false;
						}
					}
				}
			}
		}
		
		numberOfSteps = 0;
		
		System.out.println("Optimal Point: " + optimalPoint.print() + " With Risk :" + getRisk(optimalPoint) + " With Utilization" + getUtilization(optimalPoint));
		
		
		

		System.out.println ("Valid Points at Start: " +  calculateValidPoints () );

		Random randomNumberGenerator = new Random(0);
		
		int count = 0;
		while (!RandomList.isEmpty()) {
			DesignPoint randomDesignPoint = RandomList.get(randomNumberGenerator.nextInt(RandomList.size()));
			if (!randomDesignPoint.checked) {
				count++;
				//System.out.println( count + " Checking " + randomDesignPoint.print() +  " Remaining: " + RandomList.size());
				//System.out.println( count + "\t randomDesignPoint.print() +  " Remaining: " + RandomList.size());
				evaluate (randomDesignPoint);
				while (!DesignStack.empty()) {
					evaluateBlind (DesignStack.pop());
				}
			}
			
			RandomList.remove(randomDesignPoint);
			
		}
		
		
		
		
		
		
		System.out.println ("Valid Points at End: " +  calculateValidPoints () );

	}
	
	public static void evaluate (DesignPoint Point) {
		
		if (getRisk(Point) > DMR) {
			Point.MarkInvalid();
			Rule1 (Point);
			Rule3 (Point);
		}
		if (getRisk(Point) < DMR) {
			Rule2 (Point);
			Rule4 (Point);
		}
		
	}
	
	public static void evaluateBlind (DesignPoint Point) {
		
		if (Point.guessDMRRejection) {
			Point.MarkInvalid();
			Rule1 (Point);
			Rule3 (Point);
		}
		if (Point.guessUtilizationRejection) {
			Point.MarkInvalid();
			Rule2 (Point);
			Rule4 (Point);
		}
		
		
	}
	
	public static void Rule4 (DesignPoint Point) {
		
		if (Point.a>0) 
			if (Point.b<maxCores-1) {
				DesignPoint NewPoint = DesignPoints[Point.a-1][Point.b+1][Point.c][Point.d][Point.e];
				if (getUtilization(NewPoint)> getUtilization(Point)) {
					NewPoint.guessUtilizationRejection= true;
					PushToDesignStack(NewPoint,Point);
				}				
			}
		
		if (Point.a>0) 
			if (Point.c<maxCores-1) {
				DesignPoint NewPoint = DesignPoints[Point.a-1][Point.b][Point.c+1][Point.d][Point.e];
				if (getUtilization(NewPoint)> getUtilization(Point)) {
					NewPoint.guessUtilizationRejection= true;
					PushToDesignStack(NewPoint,Point);
				}				
			}
		
		if (Point.a>0) 
			if (Point.d<maxCores-1) {
				DesignPoint NewPoint = DesignPoints[Point.a-1][Point.b][Point.c][Point.d+1][Point.e];
				if (getUtilization(NewPoint)> getUtilization(Point)) {
					NewPoint.guessUtilizationRejection= true;
					PushToDesignStack(NewPoint,Point);
				}				
			}
		
		if (Point.a>0) 
			if (Point.e<maxCores-1) {
				DesignPoint NewPoint = DesignPoints[Point.a-1][Point.b][Point.c][Point.d][Point.e+1];
				if (getUtilization(NewPoint)> getUtilization(Point)) {
					NewPoint.guessUtilizationRejection= true;
					PushToDesignStack(NewPoint,Point);
				}				
			}
				
				
		
		
		
		if (Point.b>0) 
			if (Point.a<maxCores-1) {
				DesignPoint NewPoint = DesignPoints[Point.a+1][Point.b-1][Point.c][Point.d][Point.e];
				if (getUtilization(NewPoint)> getUtilization(Point)) {
					NewPoint.guessUtilizationRejection= true;
					PushToDesignStack(NewPoint,Point);
				}				
			}
		
		if (Point.b>0) 
			if (Point.c<maxCores-1) {
				DesignPoint NewPoint = DesignPoints[Point.a][Point.b-1][Point.c+1][Point.d][Point.e];
				if (getUtilization(NewPoint)> getUtilization(Point)) {
					NewPoint.guessUtilizationRejection= true;
					PushToDesignStack(NewPoint,Point);
				}				
			}
		
		if (Point.b>0) 
			if (Point.d<maxCores-1) {
				DesignPoint NewPoint = DesignPoints[Point.a][Point.b-1][Point.c][Point.d+1][Point.e];
				if (getUtilization(NewPoint)> getUtilization(Point)) {
					NewPoint.guessUtilizationRejection= true;
					PushToDesignStack(NewPoint,Point);
				}				
			}
		
		if (Point.b>0) 
			if (Point.e<maxCores-1) {
				DesignPoint NewPoint = DesignPoints[Point.a][Point.b-1][Point.c][Point.d][Point.e];
				if (getUtilization(NewPoint)> getUtilization(Point)) {
					NewPoint.guessUtilizationRejection= true;
					PushToDesignStack(NewPoint,Point);
				}				
			}
		
		
		
		if (Point.c>0) 
			if (Point.a<maxCores-1) {
				DesignPoint NewPoint = DesignPoints[Point.a+1][Point.b][Point.c-1][Point.d][Point.e];
				if (getUtilization(NewPoint)> getUtilization(Point)) {
					NewPoint.guessUtilizationRejection= true;
					PushToDesignStack(NewPoint,Point);
				}				
			}
		
		if (Point.c>0) 
			if (Point.b<maxCores-1) {
				DesignPoint NewPoint = DesignPoints[Point.a][Point.b+1][Point.c-1][Point.d][Point.e];
				if (getUtilization(NewPoint)> getUtilization(Point)) {
					NewPoint.guessUtilizationRejection= true;
					PushToDesignStack(NewPoint,Point);
				}				
			}
		
		if (Point.c>0) 
			if (Point.d<maxCores-1) {
				DesignPoint NewPoint = DesignPoints[Point.a][Point.b][Point.c-1][Point.d+1][Point.e];
				if (getUtilization(NewPoint)> getUtilization(Point)) {
					NewPoint.guessUtilizationRejection= true;
					PushToDesignStack(NewPoint,Point);
				}				
			}
		
		if (Point.c>0) 
			if (Point.e<maxCores-1) {
				DesignPoint NewPoint = DesignPoints[Point.a][Point.b][Point.c-1][Point.d][Point.e+1];
				if (getUtilization(NewPoint)> getUtilization(Point)) {
					NewPoint.guessUtilizationRejection= true;
					PushToDesignStack(NewPoint,Point);
				}				
			}
		
		
		if (Point.d>0) 
			if (Point.a<maxCores-1) {
				DesignPoint NewPoint = DesignPoints[Point.a+1][Point.b][Point.c][Point.d-1][Point.e];
				if (getUtilization(NewPoint)> getUtilization(Point)) {
					NewPoint.guessUtilizationRejection= true;
					PushToDesignStack(NewPoint,Point);
				}				
			}
		
		if (Point.d>0) 
			if (Point.b<maxCores-1) {
				DesignPoint NewPoint = DesignPoints[Point.a][Point.b+1][Point.c][Point.d-1][Point.e];
				if (getUtilization(NewPoint)> getUtilization(Point)) {
					NewPoint.guessUtilizationRejection= true;
					PushToDesignStack(NewPoint,Point);
				}				
			}
		
		if (Point.d>0) 
			if (Point.c<maxCores-1) {
				DesignPoint NewPoint = DesignPoints[Point.a][Point.b][Point.c+1][Point.d-1][Point.e];
				if (getUtilization(NewPoint)> getUtilization(Point)) {
					NewPoint.guessUtilizationRejection= true;
					PushToDesignStack(NewPoint,Point);
				}				
			}
		
		if (Point.d>0) 
			if (Point.e<maxCores-1) {
				DesignPoint NewPoint = DesignPoints[Point.a][Point.b][Point.c][Point.d-1][Point.e+1];
				if (getUtilization(NewPoint)> getUtilization(Point)) {
					NewPoint.guessUtilizationRejection= true;
					PushToDesignStack(NewPoint,Point);
				}				
			}
		
		if (Point.e>0) 
			if (Point.a<maxCores-1) {
				DesignPoint NewPoint = DesignPoints[Point.a+1][Point.b][Point.c][Point.d][Point.e-1];
				if (getUtilization(NewPoint)> getUtilization(Point)) {
					NewPoint.guessUtilizationRejection= true;
					PushToDesignStack(NewPoint,Point);
				}				
			}
		
		if (Point.e>0) 
			if (Point.b<maxCores-1) {
				DesignPoint NewPoint = DesignPoints[Point.a][Point.b+1][Point.c][Point.d][Point.e-1];
				if (getUtilization(NewPoint)> getUtilization(Point)) {
					NewPoint.guessUtilizationRejection= true;
					PushToDesignStack(NewPoint,Point);
				}				
			}
		
		if (Point.e>0) 
			if (Point.c<maxCores-1) {
				DesignPoint NewPoint = DesignPoints[Point.a][Point.b][Point.c+1][Point.d][Point.e-1];
				if (getUtilization(NewPoint)> getUtilization(Point)) {
					NewPoint.guessUtilizationRejection= true;
					PushToDesignStack(NewPoint,Point);
				}				
			}
		
		if (Point.e>0) 
			if (Point.d<maxCores-1) {
				DesignPoint NewPoint = DesignPoints[Point.a][Point.b][Point.c][Point.d+1][Point.e-1];
				if (getUtilization(NewPoint)> getUtilization(Point)) {
					NewPoint.guessUtilizationRejection= true;
					PushToDesignStack(NewPoint,Point);
				}				
			}
		
		
	}
	
	public static void Rule3 (DesignPoint Point) {
		
		
		if (Point.a>0) 
			if (Point.b<maxCores-1) {
				DesignPoint NewPoint = DesignPoints[Point.a-1][Point.b+1][Point.c][Point.d][Point.e];
				if (getExecutionTime(NewPoint)> getExecutionTime(Point)) {
					NewPoint.guessDMRRejection = true;
					PushToDesignStack(NewPoint,Point);
				}				
			}
		
		if (Point.a>0) 
			if (Point.c<maxCores-1) {
				DesignPoint NewPoint = DesignPoints[Point.a-1][Point.b][Point.c+1][Point.d][Point.e];
				if (getExecutionTime(NewPoint)> getExecutionTime(Point)) {
					NewPoint.guessDMRRejection = true;
					PushToDesignStack(NewPoint,Point);
				}				
			}
		
		if (Point.a>0) 
			if (Point.d<maxCores-1) {
				DesignPoint NewPoint = DesignPoints[Point.a-1][Point.b][Point.c][Point.d+1][Point.e];
				if (getExecutionTime(NewPoint)> getExecutionTime(Point)) {
					NewPoint.guessDMRRejection = true;
					PushToDesignStack(NewPoint,Point);
				}				
			}
		
		if (Point.a>0) 
			if (Point.e<maxCores-1) {
				DesignPoint NewPoint = DesignPoints[Point.a-1][Point.b][Point.c][Point.d][Point.e+1];
				if (getExecutionTime(NewPoint)> getExecutionTime(Point)) {
					NewPoint.guessDMRRejection = true;
					PushToDesignStack(NewPoint,Point);
				}				
			}
				
				
		
		
		
		if (Point.b>0) 
			if (Point.a<maxCores-1) {
				DesignPoint NewPoint = DesignPoints[Point.a+1][Point.b-1][Point.c][Point.d][Point.e];
				if (getExecutionTime(NewPoint)> getExecutionTime(Point)) {
					NewPoint.guessDMRRejection = true;
					PushToDesignStack(NewPoint,Point);
				}				
			}
		
		if (Point.b>0) 
			if (Point.c<maxCores-1) {
				DesignPoint NewPoint = DesignPoints[Point.a][Point.b-1][Point.c+1][Point.d][Point.e];
				if (getExecutionTime(NewPoint)> getExecutionTime(Point)) {
					NewPoint.guessDMRRejection = true;
					PushToDesignStack(NewPoint,Point);
				}				
			}
		
		if (Point.b>0) 
			if (Point.d<maxCores-1) {
				DesignPoint NewPoint = DesignPoints[Point.a][Point.b-1][Point.c][Point.d+1][Point.e];
				if (getExecutionTime(NewPoint)> getExecutionTime(Point)) {
					NewPoint.guessDMRRejection = true;
					PushToDesignStack(NewPoint,Point);
				}				
			}
		
		if (Point.b>0) 
			if (Point.e<maxCores-1) {
				DesignPoint NewPoint = DesignPoints[Point.a][Point.b-1][Point.c][Point.d][Point.e];
				if (getExecutionTime(NewPoint)> getExecutionTime(Point)) {
					NewPoint.guessDMRRejection = true;
					PushToDesignStack(NewPoint,Point);
				}				
			}
		
		
		
		if (Point.c>0) 
			if (Point.a<maxCores-1) {
				DesignPoint NewPoint = DesignPoints[Point.a+1][Point.b][Point.c-1][Point.d][Point.e];
				if (getExecutionTime(NewPoint)> getExecutionTime(Point)) {
					NewPoint.guessDMRRejection = true;
					PushToDesignStack(NewPoint,Point);
				}				
			}
		
		if (Point.c>0) 
			if (Point.b<maxCores-1) {
				DesignPoint NewPoint = DesignPoints[Point.a][Point.b+1][Point.c-1][Point.d][Point.e];
				if (getExecutionTime(NewPoint)> getExecutionTime(Point)) {
					NewPoint.guessDMRRejection = true;
					PushToDesignStack(NewPoint,Point);
				}				
			}
		
		if (Point.c>0) 
			if (Point.d<maxCores-1) {
				DesignPoint NewPoint = DesignPoints[Point.a][Point.b][Point.c-1][Point.d+1][Point.e];
				if (getExecutionTime(NewPoint)> getExecutionTime(Point)) {
					NewPoint.guessDMRRejection = true;
					PushToDesignStack(NewPoint,Point);
				}				
			}
		
		if (Point.c>0) 
			if (Point.e<maxCores-1) {
				DesignPoint NewPoint = DesignPoints[Point.a][Point.b][Point.c-1][Point.d][Point.e+1];
				if (getExecutionTime(NewPoint)> getExecutionTime(Point)) {
					NewPoint.guessDMRRejection = true;
					PushToDesignStack(NewPoint,Point);
				}				
			}
		
		
		if (Point.d>0) 
			if (Point.a<maxCores-1) {
				DesignPoint NewPoint = DesignPoints[Point.a+1][Point.b][Point.c][Point.d-1][Point.e];
				if (getExecutionTime(NewPoint)> getExecutionTime(Point)) {
					NewPoint.guessDMRRejection = true;
					PushToDesignStack(NewPoint,Point);
				}				
			}
		
		if (Point.d>0) 
			if (Point.b<maxCores-1) {
				DesignPoint NewPoint = DesignPoints[Point.a][Point.b+1][Point.c][Point.d-1][Point.e];
				if (getExecutionTime(NewPoint)> getExecutionTime(Point)) {
					NewPoint.guessDMRRejection = true;
					PushToDesignStack(NewPoint,Point);
				}				
			}
		
		if (Point.d>0) 
			if (Point.c<maxCores-1) {
				DesignPoint NewPoint = DesignPoints[Point.a][Point.b][Point.c+1][Point.d-1][Point.e];
				if (getExecutionTime(NewPoint)> getExecutionTime(Point)) {
					NewPoint.guessDMRRejection = true;
					PushToDesignStack(NewPoint,Point);
				}				
			}
		
		if (Point.d>0) 
			if (Point.e<maxCores-1) {
				DesignPoint NewPoint = DesignPoints[Point.a][Point.b][Point.c][Point.d-1][Point.e+1];
				if (getExecutionTime(NewPoint)> getExecutionTime(Point)) {
					NewPoint.guessDMRRejection = true;
					PushToDesignStack(NewPoint,Point);
				}				
			}
		
		if (Point.e>0) 
			if (Point.a<maxCores-1) {
				DesignPoint NewPoint = DesignPoints[Point.a+1][Point.b][Point.c][Point.d][Point.e-1];
				if (getExecutionTime(NewPoint)> getExecutionTime(Point)) {
					NewPoint.guessDMRRejection = true;
					PushToDesignStack(NewPoint,Point);
				}				
			}
		
		if (Point.e>0) 
			if (Point.b<maxCores-1) {
				DesignPoint NewPoint = DesignPoints[Point.a][Point.b+1][Point.c][Point.d][Point.e-1];
				if (getExecutionTime(NewPoint)> getExecutionTime(Point)) {
					NewPoint.guessDMRRejection = true;
					PushToDesignStack(NewPoint,Point);
				}				
			}
		
		if (Point.e>0) 
			if (Point.c<maxCores-1) {
				DesignPoint NewPoint = DesignPoints[Point.a][Point.b][Point.c+1][Point.d][Point.e-1];
				if (getExecutionTime(NewPoint)> getExecutionTime(Point)) {
					NewPoint.guessDMRRejection = true;
					PushToDesignStack(NewPoint,Point);
				}				
			}
		
		if (Point.e>0) 
			if (Point.d<maxCores-1) {
				DesignPoint NewPoint = DesignPoints[Point.a][Point.b][Point.c][Point.d+1][Point.e-1];
				if (getExecutionTime(NewPoint)> getExecutionTime(Point)) {
					NewPoint.guessDMRRejection = true;
					PushToDesignStack(NewPoint,Point);
				}				
			}
			
		
	}
	
	public static void Rule1 (DesignPoint Point) {
		
		if (Point.a > 0) {
			DesignPoints[Point.a-1][Point.b][Point.c][Point.d][Point.e].guessDMRRejection = true;
			PushToDesignStack(DesignPoints[Point.a-1][Point.b][Point.c][Point.d][Point.e],Point);
		}
		if (Point.b > 0) {
			DesignPoints[Point.a][Point.b-1][Point.c][Point.d][Point.e].guessDMRRejection = true;
			PushToDesignStack(DesignPoints[Point.a][Point.b-1][Point.c][Point.d][Point.e],Point);
		}
		if (Point.c > 0) {
			DesignPoints[Point.a][Point.b][Point.c-1][Point.d][Point.e].guessDMRRejection = true;
			PushToDesignStack(DesignPoints[Point.a][Point.b][Point.c-1][Point.d][Point.e],Point);
		}
		if (Point.d > 0) {
			DesignPoints[Point.a][Point.b][Point.c][Point.d-1][Point.e].guessDMRRejection = true;
			PushToDesignStack(DesignPoints[Point.a][Point.b][Point.c][Point.d-1][Point.e],Point);
		}
		if (Point.e > 0) {
			DesignPoints[Point.a][Point.b][Point.c][Point.d][Point.e-1].guessDMRRejection = true;
			PushToDesignStack(DesignPoints[Point.a][Point.b][Point.c][Point.d][Point.e-1],Point);
		}
		
	}
	
	public static void Rule2 (DesignPoint Point) {
		
		if (Point.a < maxCores-1) {
			DesignPoints[Point.a+1][Point.b][Point.c][Point.d][Point.e].guessUtilizationRejection = true;
			PushToDesignStack(DesignPoints[Point.a+1][Point.b][Point.c][Point.d][Point.e],Point);
		}
		if (Point.b < maxCores-1) {
			DesignPoints[Point.a][Point.b+1][Point.c][Point.d][Point.e].guessUtilizationRejection = true;
			PushToDesignStack(DesignPoints[Point.a][Point.b+1][Point.c][Point.d][Point.e],Point);
		}
		if (Point.c < maxCores-1) {
			DesignPoints[Point.a][Point.b][Point.c+1][Point.d][Point.e].guessUtilizationRejection = true;
			PushToDesignStack(DesignPoints[Point.a][Point.b][Point.c+1][Point.d][Point.e],Point);
		}
		if (Point.d < maxCores-1) {
			DesignPoints[Point.a][Point.b][Point.c][Point.d+1][Point.e].guessUtilizationRejection = true;
			PushToDesignStack(DesignPoints[Point.a][Point.b][Point.c][Point.d+1][Point.e],Point);
		}
		if (Point.e < maxCores-1) {
			DesignPoints[Point.a][Point.b][Point.c][Point.d][Point.e+1].guessUtilizationRejection = true;
			PushToDesignStack(DesignPoints[Point.a][Point.b][Point.c][Point.d][Point.e+1],Point);
		}
		
	}
	
	
	public static void PushToDesignStack (DesignPoint Point, DesignPoint Parent) {
		
			//System.out.println ("Trying to Push Point: " + Point.print() );
			if (!Point.checked && !DesignStack.contains(Point)) {
				Point.PushedBy = Parent;
				Point.checked = true;
				DesignStack.push(Point);
			}
		
	}
	
	
	public static int calculateValidPoints () {
		
		int validCount = 0;
		
		for (int a = 0; a <maxCores; a++) {
			for (int b = 0; b<maxCores;b++) {
				for (int c = 0; c < maxCores; c++) {
					for (int d = 0; d < maxCores;d++) {
						for (int e = 0; e < maxCores; e++) {
							if (DesignPoints [a] [b] [c] [d] [e].valid) validCount++;
						}
					}
				}
			}
		}
		
		return validCount;
		
		
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

	


