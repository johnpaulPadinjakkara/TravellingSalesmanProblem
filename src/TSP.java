import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.StringReader;
import java.util.ArrayList;
import java.util.Random;
import java.util.Scanner;
import java.util.stream.Collectors;
import java.util.stream.Stream;

public class TSP {

	int mutationRate;
	int crossoverRate;
	static int pointCount;
	static String dataSet = null;
	static int fileOption = 0;

	static ArrayList<ArrayList<Integer>> initialPopulation = new ArrayList<ArrayList<Integer>>();
	static ArrayList<Integer> resultOfCyclicCX = new ArrayList<Integer>();

	ArrayList<ArrayList<Integer>> initializePopulation(int populationSize, int pointCount) {

		int count = 0;
		int repetationCheck;

		while (count < populationSize) {
			Random rand = new Random();
			ArrayList<Integer> routes = new ArrayList<Integer>();
			for (int i = 1; i <= pointCount; i++) {
				repetationCheck = 0;
				int n = rand.nextInt(pointCount) + 1;
				for (int j = 0; j < routes.size(); j++) {
					if (n == routes.get(j)) {
						repetationCheck++;
					}
				}
				if (repetationCheck == 0) {
					routes.add(n);
				} else {
					i--;
				}
			}
			initialPopulation.add(routes);
			count++;
		}
		return initialPopulation;
	}

	public static String readFile() throws IOException {
		String everything;
		
		BufferedReader br = null;
		
		if(TSP.fileOption == 0) {
		 br = new BufferedReader(new FileReader("A2_TSP_WesternSahara_29.txt"));
		
		}
		else {
			 br = new BufferedReader(new FileReader("A2_TSP_Uruguay_734.txt"));
		}
			try {
				StringBuilder sb = new StringBuilder();
				String line = br.readLine();
				while (line != null) {
					sb.append(line);
					sb.append(System.lineSeparator());
					line = br.readLine();
					pointCount++;
				}
				everything = sb.toString();
			} finally {
				br.close();
		}
	

		return everything;
	}

	public static double fitnessCalc(ArrayList<Integer> chromosome) {
		double fitness = 0;
		for (int i = 0; i < pointCount - 1; i++) {
			fitness = fitness + TSP.distanceCalc(chromosome.get(i), chromosome.get(i + 1));
		}

		fitness = fitness + TSP.distanceCalc(chromosome.get(chromosome.size() - 1), chromosome.get(0));
		return fitness;
	}

	public static double distanceCalc(int point1, int point2) {
		double distance = 0;
		double latitude_1;
		double longitude_1;
		double latitude_2;
		double longitude_2;
		String[] linesInDataSet = TSP.dataSet.split("\\n+");
		String[] point1DataSet = linesInDataSet[point1 - 1].split("\\s+");
		String[] point2DataSet = linesInDataSet[point2 - 1].split("\\s+");
		latitude_1 = Double.parseDouble(point1DataSet[1]);
		longitude_1 = Double.parseDouble(point1DataSet[2]);
		latitude_2 = Double.parseDouble(point2DataSet[1]);
		longitude_2 = Double.parseDouble(point2DataSet[2]);
		distance = TSP.distanceInKmBetweenEarthCoordinates(latitude_1, longitude_1, latitude_2, longitude_2);
		return distance;
	}

	public static double distanceInKmBetweenEarthCoordinates(double lat1, double lon1, double lat2, double lon2) {
		double dLat = lat2 - lat1;
		double dLon = lon2 - lon1;
		dLat = dLat * dLat;
		dLon = dLon * dLon;
		double distance;
		double temp;
		temp = Math.abs(dLat + dLon);
		distance = Math.sqrt(temp);
		return distance;
	}

	public static ArrayList<Integer> selectionTournament(ArrayList<ArrayList<Integer>> population, int rand,
			int popSize) {
		ArrayList<Integer> best = new ArrayList<Integer>();
		best.clear();
		best = population.get(rand);
		ArrayList<Integer> temp = new ArrayList<Integer>();
		temp.clear();
		
		Random r = new Random();
		for (int i = 1; i <= rand; i++) {
			int rr = r.nextInt(popSize);
			temp = population.get(rr);
			if (best == null) {
				best = temp;
			} else if (TSP.fitnessCalc(temp) < TSP.fitnessCalc(best)) {
				best = temp;
			} else {
			}
		}
		return best;
	}
	//---------------
	public static ArrayList<Integer> Cyclic_Crossover(ArrayList<Integer> p1, ArrayList<Integer> p2,ArrayList<Integer> o1, ArrayList<Integer> o2, int initialParentSize) {
		ArrayList<Integer> offspring1 = new ArrayList<Integer>();
		ArrayList<Integer> offspring2 = new ArrayList<Integer>();
		ArrayList<Integer> Finaloffspring1 = new ArrayList<Integer>();
		ArrayList<Integer> Finaloffspring2 = new ArrayList<Integer>();
		int initialSize  = initialParentSize;


		offspring1.clear();
		offspring2.clear();
		if(p1.size()==p2.size()) {
		
		Finaloffspring1 = (ArrayList<Integer>) Stream.of(o1, offspring1) .flatMap(x -> x.stream()).collect(Collectors.toList());
		Finaloffspring2 = (ArrayList<Integer>) Stream.of(o2, offspring2) .flatMap(x -> x.stream()).collect(Collectors.toList());

		
		System.out.println("P1  "+p1);
		System.out.println("P2  "+p2);
		
		
		offspring1.add(p2.get(0));
		offspring2.add(p2.get(p1.indexOf(p2.get(p1.indexOf(p2.get(0))))));
		
		for(int i = 0 ; i < (p1.size()-1);i++) {
			
		offspring1.add(p2.get(p1.indexOf(offspring2.get(i))));
		offspring2.add(p2.get(p1.indexOf(p2.get(p1.indexOf(p2.get(p1.indexOf(offspring2.get(i))))))));
		
		
		if(offspring2.contains(p1.get(0))) {
			
			if(offspring1.size() == p1.size()) {
				System.out.println(" End of Crossover operation   ");
				Finaloffspring1 = (ArrayList<Integer>) Stream.of(o1, offspring1) .flatMap(x -> x.stream()).collect(Collectors.toList());
				Finaloffspring2 = (ArrayList<Integer>) Stream.of(o2, offspring2) .flatMap(x -> x.stream()).collect(Collectors.toList());
				System.out.println("Final o1  "+Finaloffspring1);
				System.out.println("Final o2  "+Finaloffspring2);
				break;

			}else {
				System.out.println("Case 2 detected  ");

				Finaloffspring1 = (ArrayList<Integer>) Stream.of(o1, offspring1) .flatMap(x -> x.stream()).collect(Collectors.toList());
				Finaloffspring2 = (ArrayList<Integer>) Stream.of(o2, offspring2) .flatMap(x -> x.stream()).collect(Collectors.toList());
				System.out.println("Final o1  "+Finaloffspring1);
				System.out.println("Final o2  "+Finaloffspring2);

				ArrayList<Integer> parent1ForRecrusion = new ArrayList<Integer>();
				ArrayList<Integer> parent2ForRecrusion = new ArrayList<Integer>();
				for(int l : p1) {
					if(!offspring2.contains(l)) {
						parent1ForRecrusion.add(l);
					}
				}
				for(int l : p2) {
					if(!offspring1.contains(l)) {
						parent2ForRecrusion.add(l);
					}
				}

				TSP.Cyclic_Crossover(parent1ForRecrusion, parent2ForRecrusion, Finaloffspring1, Finaloffspring2,initialSize);

				break;

			}
		}
		}
		

	

		System.out.println("Final Return  "+Finaloffspring1);

		System.out.println("Init Size  "+initialSize);
		System.out.println("Final Size  "+Finaloffspring1.size());

		Finaloffspring1 = (ArrayList<Integer>) Stream.of(o1, offspring1) .flatMap(x -> x.stream()).collect(Collectors.toList());
		Finaloffspring2 = (ArrayList<Integer>) Stream.of(o2, offspring2) .flatMap(x -> x.stream()).collect(Collectors.toList());
		if(initialSize == Finaloffspring1.size()) {
			
			resultOfCyclicCX = Finaloffspring1;
		}
		return Finaloffspring1;
		}
		else {
			return p1;
		}
	}
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	//-----------

	// Ordered crossover
	public static ArrayList<Integer> OX_Crossover(ArrayList<Integer> c1, ArrayList<Integer> c2) {
		ArrayList<Integer> offspring = new ArrayList<Integer>();
		offspring.clear();

		Random r = new Random();
		// int lowerLimit = 2;
		// int upperLimit = 6;
		int lowerLimit = r.nextInt(c1.size() / 2) + 1;
		int upperLimit = lowerLimit + Math.abs(c1.size() / 4);
		for (int x = 0; x < c1.size(); x++) {
			offspring.add(0);
		}
		System.out.println("Lower = " + lowerLimit);
		System.out.println("Upper = " + upperLimit);

		for (int i = lowerLimit; i < upperLimit; i++) {
			offspring.set(i, c2.get(i));
		}
		System.out.println("\nOrdered Crossover step 1 " + offspring + "");

		int k = upperLimit;
		int i = upperLimit;
		while (k < c1.size()) {

			if (offspring.contains(c1.get(i % c1.size()))) {
				i++;
			} else {
				offspring.set(k, c1.get(i % c1.size()));
				i++;
				k++;
			}

			// System.out.println("k = "+k +"i = "+i);

		}
		System.out.println("\nOrdered Crossover step 2 " + offspring + "");

		k = 0;
		i = 0;
		while (k < lowerLimit) {

			if (offspring.contains(c1.get(i % c1.size()))) {
				i++;
			} else {
				offspring.set(k, c1.get(i % c1.size()));
				i++;
				k++;
			}

			// System.out.println("k = "+k +"i = "+i);

		}
		System.out.println("\nOrdered Crossover step 3 " + offspring + "");

	
	

		return offspring;

	}

	public static ArrayList<Integer> PMX_Crossover(ArrayList<Integer> c1, ArrayList<Integer> c2) {
		ArrayList<Integer> offspring = new ArrayList<Integer>();
		offspring.clear();
		Random r = new Random();
		int count = 0;
		int index = 0;
		int index1 = 0;

		int lowerLimit = r.nextInt(c1.size() / 2) + 1;
		int upperLimit = Math.abs(c1.size() / 4);
		for (int x = 0; x < c1.size(); x++) {
			offspring.add(0);
		}
		for (int i = 1; i <= upperLimit; i++) {
			offspring.set(lowerLimit + (i - 1), c1.get(lowerLimit + (i - 1)));
		}
		System.out.println("\nPMX Crossover step 1 " + offspring + "");

		for (int i = lowerLimit; i < lowerLimit + upperLimit; i++) {
			if (offspring.contains(c2.get(i))) {
				count++;
			} else {
				for (int k = 0; k < c1.size(); k++) {
					if (c1.get(i) == c2.get(k)) {
						index = k;
					}
				}
				if (offspring.get(index) == 0) {
					offspring.set(index, c2.get(i));
				}

				else {

					int element_conflict = offspring.get(index);

					index1 = c2.indexOf(element_conflict);

					int index_of_conflict = c2.indexOf(c1.get(index1));
					offspring.set(index_of_conflict, c2.get(i));

				}
			}

		}
		System.out.println("PMX Crossover step 2 " + offspring + "");
		for (int l = 0; l < c1.size(); l++) {
			if (offspring.get(l) == 0) {
				offspring.set(l, c2.get(l));
			}
		}
		System.out.println("PMX Crossover step 3 " + offspring);
		return offspring;
	}

	public static ArrayList<Integer> Mutation(ArrayList<Integer> c) {
		ArrayList<Integer> result = new ArrayList<Integer>();
		result.clear();
		for (int x = 0; x < c.size(); x++) {
			result.add(c.get(x));
		}
		int swap_loc1;
		int swap_loc2;
		int temp;
		Random r = new Random();
		swap_loc1 = r.nextInt(c.size());
		swap_loc2 = r.nextInt(c.size());
		temp = c.get(swap_loc1);
		result.set(swap_loc1, c.get(swap_loc2));
		result.set(swap_loc2, temp);
		System.out.println("Mutated              " + result + "Fitness " + TSP.fitnessCalc(result));
		return result;
	}

	public static void main(String[] args) throws IOException {
		Random r = new Random();
		Random rRate = new Random();
		Scanner readerr = new Scanner(System.in);
		
		
		//Cyclic Cx test
				ArrayList<Integer> p1 = new ArrayList<Integer>();
				ArrayList<Integer> p2 = new ArrayList<Integer>();
				ArrayList<Integer> p3 = new ArrayList<Integer>();
				ArrayList<Integer> p4 = new ArrayList<Integer>();
				ArrayList<Integer> o1 = new ArrayList<Integer>();
				ArrayList<Integer> o2 = new ArrayList<Integer>();
				p1.add(3);
				p1.add(4);
				p1.add(8);
				p1.add(2);
				p1.add(7);
				p1.add(1);
				p1.add(6);
				p1.add(5);
				
				p2.add(4);
				p2.add(2);
				p2.add(5);
				p2.add(1);
				p2.add(6);
				p2.add(8);
				p2.add(3);
				p2.add(7);

				p3.add(1);
				p3.add(2);
				p3.add(3);
				p3.add(4);
				p3.add(5);
				p3.add(6);
				p3.add(7);
				p3.add(8);
				
				p4.add(2);
				p4.add(7);
				p4.add(5);
				p4.add(8);
				p4.add(4);
				p4.add(1);
				p4.add(6);
				p4.add(3);
				
				
				
//				for(int i=0;i<p1.size();i++) {
//					o1.add(0);
//					o2.add(0);
//				}
				
				TSP.Cyclic_Crossover(p1, p2,o1,o2,p1.size());
				System.out.println("Result of Cyclic CX 1 "+ resultOfCyclicCX);
				 TSP.Cyclic_Crossover(p3, p4,o1,o2,p1.size());
				System.out.println("Result of Cyclic CX 2  "+ resultOfCyclicCX);




				
				//-----

		System.out.println("PLEASE SELECT THE FILE TO RUN \n");
		System.out.println("Press 0 for A2_TSP_WesternSahara_29.txt \nPress 1 for A2_TSP_Uruguay_734.txt");
		TSP.fileOption = readerr.nextInt();
		System.out.println("Data in the selected file \n");
		TSP.dataSet = readFile();
		System.out.println(TSP.dataSet);

		ArrayList<Integer> selectedParent1 = new ArrayList<Integer>();
		ArrayList<Integer> selectedParent2 = new ArrayList<Integer>();
		ArrayList<Integer> offspring = new ArrayList<Integer>();
		ArrayList<Integer> WorstOffspring_Replace = new ArrayList<Integer>();
		ArrayList<Integer> BestOffspring = new ArrayList<Integer>();

		ArrayList<Integer> bestRouteAtLast = new ArrayList<Integer>();
		ArrayList<Integer> crossoverRateList = new ArrayList<Integer>();
		ArrayList<Integer> mutationRateList = new ArrayList<Integer>();
		
		ArrayList<Double> totalFitness = new ArrayList<Double>();
		ArrayList<Double> averageFitnessValues = new ArrayList<Double>();
		ArrayList<Double> bestFitness = new ArrayList<Double>();
		ArrayList<Double> worstFitness = new ArrayList<Double>();





		int genCount = 200;
		int worst_fitness_at_loc = 0;
		int best_fitness_at_loc = 0;

		double crossoverRate = 1;
		double mutationRate = 1;
		double averageFitness = 0;
		double currentChromosomeFitness = 0;
		
		

		System.out.println("Number of points read from file = " + TSP.pointCount + "\n");
		int popSize = 0;
		TSP tsp = new TSP();
		
		System.out.print("Enter the initial population size ");
		popSize = readerr.nextInt();
		
		System.out.print("PLEASE ENTER THE CROSSOVER RATE (enter a value between 0.0 to 1.0 ) \n ");
		crossoverRate = readerr.nextDouble();
		
		System.out.print("PLEASE ENTER THE MUTATION RATE (enter a value between 0.0 to 1.0 ) \n ");
		mutationRate = readerr.nextDouble();
		readerr.close();

		for (int i = 0; i < pointCount; i++) {
		}

		// Create initial population
		tsp.initializePopulation(popSize, TSP.pointCount);

		for (int i = 0; i < initialPopulation.size(); i++) {
			System.out.println(
					TSP.initialPopulation.get(i) + "   Fitness = " + TSP.fitnessCalc(initialPopulation.get(i)));
		}
		System.out.println("\n");
		crossoverRateList.add(0);
		int m = 0;
		// int temp;
		int limit = (int) Math.round(Math.abs((genCount - (crossoverRate * genCount))));
		System.out.println("L " + limit);
		while (m < limit) {
			int temp = rRate.nextInt(genCount + 1);
			if (!(crossoverRateList.contains(temp))) {
				// System.out.println(rRate.nextInt(genCount+1));
				crossoverRateList.add(temp);
				m = m + 1;
				// System.out.println(temp);

			} else {

			}

		}
		crossoverRateList.remove(0);

		System.out.println("crx list = " + crossoverRateList);

		mutationRateList.add(0);

		int mut = 0;
		int limitMut = (int) Math.round(Math.abs((genCount - (mutationRate * genCount))));
		System.out.println("LM " + limitMut);
		while (mut < limitMut) {
			int tempMut = rRate.nextInt(genCount + 1);
			if (!(mutationRateList.contains(tempMut))) {
				// System.out.println(rRate.nextInt(genCount+1));
				mutationRateList.add(tempMut);
				mut = mut + 1;
				// System.out.println(temp);

			} else {

			}

		}
		mutationRateList.remove(0);

		System.out.println("mut list = " + mutationRateList);
		// Initial population Created
		// GA loop with Selection, PMX Crossover, Mutation(Swap), Replacement in
		// population.
		while (genCount >= 0) {
			// PARENT SELECTION BY TOURNAMENT SELECTION
			// rand should be only from 1 to (popSize - 1)
			selectedParent1 = TSP.selectionTournament(TSP.initialPopulation, r.nextInt(popSize), popSize);
			selectedParent2 = TSP.selectionTournament(TSP.initialPopulation, r.nextInt(popSize), popSize);

			System.out.println("parent 1          " + selectedParent1 + "Fitness " + TSP.fitnessCalc(selectedParent1));
			System.out.println("parent 2          " + selectedParent2 + "Fitness " + TSP.fitnessCalc(selectedParent2));

			if ((crossoverRateList.contains(genCount + 1))) {
				// Perform Swap MUTATION
				offspring = TSP.Mutation(selectedParent1);
				System.out.println(
						"----------------------CROSSOVER RATE :: NO CROSSOVER DONE ONLY MUTATION ------------------------------------------------------------------------------------");

			} else if (mutationRateList.contains(genCount + 1)) {

				System.out.println(
						"----------------------MUTATION RATE :: NO MUTATION DONE ONLY CROSSOVER ------------------------------------------------------------------------------------");
//				// Perform OX CROSSOVER
//				offspring = TSP.OX_Crossover(selectedParent1, selectedParent2);
				
				//Perform Cyclic Crossover
				TSP.Cyclic_Crossover(selectedParent1, selectedParent2,o1,o2,selectedParent1.size());
				offspring = resultOfCyclicCX;
				System.out
						.println("\nOffspring            " + offspring + "Fitness " + TSP.fitnessCalc(selectedParent2));

			} else {
				System.out.println(
						"----------------------CROSSOVER & MUTATION ------------------------------------------------------------------------------------");

				//Perform Cyclic Crossover
				TSP.Cyclic_Crossover(selectedParent1, selectedParent2,o1,o2,selectedParent1.size());
				offspring = resultOfCyclicCX;
				System.out
						.println("\nOffspring            " + offspring + "Fitness " + TSP.fitnessCalc(selectedParent2));
				// Perform Swap MUTATION
				offspring = TSP.Mutation(selectedParent1);
			}

			

			// CODE TO PERFORM REPLACEMENT
			WorstOffspring_Replace = initialPopulation.get(0);
			worst_fitness_at_loc = 0;
			for (int i = 0; i < initialPopulation.size(); i++) {
				if (TSP.fitnessCalc(WorstOffspring_Replace) < TSP.fitnessCalc(initialPopulation.get(i))) {
					WorstOffspring_Replace = initialPopulation.get(i);
					worst_fitness_at_loc = i;
				} else {

				}
			}
			System.out.println(" offspring  to be replaced" + WorstOffspring_Replace + "Fitness "
					+ TSP.fitnessCalc(WorstOffspring_Replace) + " Located at " + worst_fitness_at_loc);

				initialPopulation.set(worst_fitness_at_loc, offspring);

				
			System.out.println("\nGeneration " + genCount);
			for (int i = 0; i < initialPopulation.size(); i++) {
				System.out.println(
						TSP.initialPopulation.get(i) + "   Fitness = " + TSP.fitnessCalc(initialPopulation.get(i)));
			}
			BestOffspring = initialPopulation.get(0);
			best_fitness_at_loc = 0;
			for (int i = 0; i < initialPopulation.size(); i++) {
				if (TSP.fitnessCalc(BestOffspring) > TSP.fitnessCalc(initialPopulation.get(i))) {
					BestOffspring = initialPopulation.get(i);
					best_fitness_at_loc = i;
				} else {

				}
			}
			System.out.println("The best fitness from the current generation  = " + BestOffspring + " \nFitness = "+TSP.fitnessCalc(BestOffspring));
			bestFitness.add(TSP.fitnessCalc(BestOffspring));
			System.out.println("The Worst fitness from the current generation = " + WorstOffspring_Replace + " \nFitness = "+TSP.fitnessCalc(WorstOffspring_Replace));
			worstFitness.add(TSP.fitnessCalc(WorstOffspring_Replace));

			// Average Fitness of the Current generation
						currentChromosomeFitness = 0;
						averageFitness = 0;
						for (int i = 0; i < initialPopulation.size(); i++) {
							currentChromosomeFitness = currentChromosomeFitness + TSP.fitnessCalc(initialPopulation.get(i));

						}
						System.out.println("TOTAL FITNESS OF CURRENT GENERATION = " + currentChromosomeFitness);
						totalFitness.add(currentChromosomeFitness);

						averageFitness = currentChromosomeFitness / initialPopulation.size();

						System.out.println("AVERAGE FITNESS OF CURRENT GENERATION = " + averageFitness);
						averageFitnessValues.add(averageFitness);

			genCount--;

		}

		bestRouteAtLast = initialPopulation.get(0);
		System.out.println("THE BEST ROUTE FROM THE FINAL GENERATION ");
		for (int j = 0; j < initialPopulation.size(); j++) {
			if (TSP.fitnessCalc(bestRouteAtLast) > TSP.fitnessCalc(initialPopulation.get(j))) {
				bestRouteAtLast = initialPopulation.get(j);
			}
		}
		
		//Data generated for GRAPH
		System.out.println(" FINAL ROUTE\n" + bestRouteAtLast + "Fitness " + TSP.fitnessCalc(bestRouteAtLast));
		System.out.println("Total "+totalFitness);
		System.out.println("aver "+averageFitnessValues);
		System.out.println("best "+bestFitness);
		System.out.println("worst "+worstFitness);




		
	}

}
