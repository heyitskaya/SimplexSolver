
import java.util. *;
import org.ejml.simple.SimpleMatrix;
import org.ejml.data.DMatrixRMaj;
import org.ejml.dense.row.CommonOps_DDRM;
import java.util.Iterator;
/** This class uses double values to give greater precision, it assumes that the user will input the LP in 
 * standard form
 * @author kayani
 *
 */
public class SimplexSolver{
	int numRows=0;
	int numCols=0;
	DMatrixRMaj internalMatrix;
	HashMap<Integer,SimpleMatrix> colMap = new HashMap<Integer,SimpleMatrix>();
	HashMap<Integer, Integer> objectFuncCostMap;
	static List<Set<Integer>> possibleCombinations;
	int currCombo=0;
	SimpleMatrix B;
	SimpleMatrix BInverse;
	int[] basicVars;
	int[] nonBasicVars;
	int enteringVar;
	int leavingVar;
	double[][] bArr;
	SimpleMatrix b ;
	
	public SimplexSolver(double[][] arr, double[][] bArr, HashMap<Integer,Integer> objectFuncCostMap) {
		this.objectFuncCostMap = objectFuncCostMap;
		 internalMatrix= new DMatrixRMaj(arr);
		 this.bArr = bArr;
		 b = new SimpleMatrix(bArr);
		 this.numRows = arr.length;
		 this.numCols = arr[0].length;
		 nonBasicVars = new int[numCols-numRows];
		 basicVars = new int[numRows];
		 SimpleMatrix simple = SimpleMatrix.wrap(internalMatrix);
		 for(int i=1;i<=numCols;i++) {
			 colMap.put(i, simple.extractVector(false, i-1));
		 }
		
		 ArrayList<Integer> superSet = new ArrayList<Integer>();
		 for(int i=1;i<=numCols;i++) {
			 superSet.add(i);
		 }
		 getSubsets(superSet,numRows);
		 objectFuncCostMap.put(5,200);
		 objectFuncCostMap.put(6,300);
		 System.out.println("The initial BFS looks like");
		 findBFS();
		 double reducedCost = getReducedCost();
		
		 while(reducedCost>0) {
			 findEnteringVariable();
			 findLeavingVariable();
			 //update xB and B matrix and 
			 updateBasicVars(leavingVar, enteringVar);
			 updateNonBasicVars(leavingVar, enteringVar);
			 //once you've updated non basic variables and basic variables 
			 //so far B is only created in findBFS which is a problem because we won't call findBFS again
			 constructBMatrix(basicVars);
			 reducedCost = getReducedCost();
			 
		 }
		 findObjectFunctionValue();
	}
	
	/** This is the last step in the algorithm. Once we've exited our while loop that means there 
	 * is no good reduced cost, therefore we calculate the reduced cost using Cb * Cx
	 * @return the optimal value 
	 */
	private double findObjectFunctionValue() {
		double value = 0 ;
		SimpleMatrix bArrMatrix = new SimpleMatrix(bArr);
		SimpleMatrix Xb = BInverse.mult(bArrMatrix);
		for(int i=0;i<numRows;i++) {
			if(objectFuncCostMap.get(basicVars[i])!=null) {
				value = value + Xb.get(i, 0) * objectFuncCostMap.get(basicVars[i]);
			}
			
		}
		System.out.println("objectiveFunctionValue is " + value);
		return value;
		
		
	}
	/** For debugging reasons to print out what our current BFS looks like **/
	private void printBFS() {
		SimpleMatrix m1 = new SimpleMatrix(numRows,numRows);
		for(int i=0;i<basicVars.length;i++) {
			m1.insertIntoThis(0, i, colMap.get(basicVars[i]));
		}
		
		m1.print();
		
	}
	
	/**Helper method for finding our initial BFS. Worst case, we'll need to consider all subsets of size numRows so here we recursively generate
	 * them all and store it in an ArrayList. Once the ArrayList is calculated it will be stored as a global variable
	 * for us to iterate through to find our initial BFS. This is only called at the beginning of the algorithm.
	 * @param superSet
	 * @param k
	 * @param index
	 * @param current
	 * @param solution
	 */
	private static void getSubsets(ArrayList<Integer> superSet, int k, int index, Set<Integer> current,List<Set<Integer>> solution) {
	    if (current.size() == k) {
	        solution.add(new HashSet<>(current));
	        return;
	    }
	    if (index == superSet.size()) return;
	    Integer x = superSet.get(index);
	    current.add(x);
	    //"guess" x is in the subset
	    getSubsets(superSet, k, index+1, current, solution);
	    current.remove(x);
	    //"guess" x is not in the subset
	    getSubsets(superSet, k, index+1, current, solution);
	}

	/** Wrapper method that calls helper method to create ArrayList of all possible subsets of size k
	 * 
	 * @param superSet
	 * @param k
	 * @return a list of subsets of size k
	 */
	public static List<Set<Integer>> getSubsets(ArrayList<Integer> superSet, int k) {
	    possibleCombinations = new ArrayList<>();
	    getSubsets(superSet, k, 0, new HashSet<Integer>(), possibleCombinations);
	    return possibleCombinations;
	}

	/** This method finds our initial BFS.
	 * It iterates through the arraylist of subsets and constructs square matrices based on the variable in the subsets
	 * Then it finds the determinant through an api call to determine whether the vectors in the square matrix are linearly 
	 * independent. If so then we initialize our BFS.
	 */
	public void findBFS() {
		boolean found = false;
		while(!found && currCombo<=possibleCombinations.size()) {
		//iterate through possible combinations
		//the determinant will be one if all vectors are linearly independent
			Set<Integer> combo1 = possibleCombinations.get(currCombo); //1,2
			int[] arr = new int[combo1.size()];
			int index =0;
			for(Integer i:combo1) {
				arr[index]=i;
				index++;
			}
			Iterator<Integer> it =combo1.iterator();
			//construct our square matrix
			double[][] squareArr = new double[numRows][numRows];
		
			SimpleMatrix m1 = new SimpleMatrix(numRows,numRows);
			for(int i=0;i<arr.length;i++) {
				m1.insertIntoThis(0,i, colMap.get(arr[i]));
			}
			if(CommonOps_DDRM.det(m1.getMatrix())!=0) {
				//then they are linearly independent which means we've
				//currently found our BFS
				B=m1;
				found = true;
				currCombo++;
				B.print(); 
			
				int pos1 = 0;
				for(int i=1;i<=numCols;i++) {
					if(!combo1.contains(i)) {
						nonBasicVars[pos1] = i;
						pos1++;
					}
				}
				int pos2 = 0;
				for(int i=1;i<=numCols;i++) {
					if(combo1.contains(i)) {
						basicVars[pos2] = i;
						pos2++;
					}
				}
				return;
			}
			else {
				currCombo++;
			}
		}
		
	}
	/**Constructing a B matrix after we've update basic variables to include entering variable
	 * @param basicVars
	 * @return The new BFS
	 */
	public SimpleMatrix constructBMatrix(int[] basicVars) {
		SimpleMatrix m1 = new SimpleMatrix(numRows,numRows);
		
		for(int i=0;i<basicVars.length;i++) {
			
			m1.insertIntoThis(0,i, colMap.get(basicVars[i]));
		}
		B=m1;
		System.out.println("The new BFS looks like");
		B.print();
		return B;
	}
	
	/** To find entering variable we first compute the reduced cost of each non basic variable and pick the largest one
	 * greater than 0
	 * Stores entering variable in global variable enteringVar
	 * 
	 */
	public void findEnteringVariable() {
		//first compute reduced cost of each non basic variable and then pick the 
		//biggest one 
		double greatestCost = Integer.MIN_VALUE;
		for(int i=0;i<nonBasicVars.length;i++) {
			int currVar = nonBasicVars[i];
			int objectiveFunctionCost;
			if(objectFuncCostMap.get(currVar)==null) {
				objectiveFunctionCost = 0;
			}
			else {
				objectiveFunctionCost = objectFuncCostMap.get(currVar);
			}
			SimpleMatrix cbTranspose = new SimpleMatrix(1,numRows);
			int currCol = 0;
			int cost;
			SimpleMatrix sm;
			double[][] d = new double[1][1];
			//creating cbTranspose
			for(int j=0;j<basicVars.length;j++) {
				if(objectFuncCostMap.get(basicVars[j])==null){
					cost = 0;
				}
				else {
					cost = objectFuncCostMap.get(basicVars[j]);
				}
				d[0][0] = cost;
				sm = new SimpleMatrix(d);
				cbTranspose.insertIntoThis(0, currCol, sm);
				currCol++;
				
			}
	
			//create BInverse
			BInverse = B.invert();
			SimpleMatrix tempMatrix1 = cbTranspose.mult(BInverse);
			SimpleMatrix tempMatrix2 = tempMatrix1.mult(colMap.get(currVar));
			double currReducedCost = objectiveFunctionCost-tempMatrix2.get(0,0);
			if(currReducedCost>greatestCost && currReducedCost > 0) {
				greatestCost = currReducedCost;
				enteringVar = currVar;
			}
			
		}
		System.out.println("entering vector we found was "+ colMap.get(enteringVar));
	}
	
	/** Iterates through array of non basic variables and calculates the reduced cost for each one using formula
	 * and returns the greatest reduced cost
	 * Used for determining whether to stay in while loop or exit while loop to calculate optimal value
	 * 
	 * @return the greatest reduced cost
	 */
	public double getReducedCost() {
	//	System.out.println("In findEnteringVariable");
		//first compute reduced cost of each non basic variable and then pick the 
		//biggest one 
		double greatestCost = Integer.MIN_VALUE;
		for(int i=0;i<nonBasicVars.length;i++) {
			int currVar = nonBasicVars[i];
			int objectiveFunctionCost;
			if(objectFuncCostMap.get(currVar)==null) {
				objectiveFunctionCost = 0;
			}
			else {
				objectiveFunctionCost = objectFuncCostMap.get(currVar);
			}
			SimpleMatrix cbTranspose = new SimpleMatrix(1,numRows);
			
			int currCol = 0;
			int cost;
			SimpleMatrix sm;
			double[][] d = new double[1][1];
			//creating cbTranspose
			for(int j=0;j<basicVars.length;j++) {
				if(objectFuncCostMap.get(basicVars[j])==null){
					cost = 0;
				}
				else {
					cost = objectFuncCostMap.get(basicVars[j]);
				}
				d[0][0] = cost;
				sm = new SimpleMatrix(d);
				cbTranspose.insertIntoThis(0, currCol, sm);
				currCol++;
				
			}
			BInverse = B.invert();
			SimpleMatrix tempMatrix1 = cbTranspose.mult(BInverse);
			SimpleMatrix tempMatrix2 = tempMatrix1.mult(colMap.get(currVar));
			double currReducedCost = objectiveFunctionCost-tempMatrix2.get(0,0);
			if(currReducedCost>greatestCost ) {
				greatestCost = currReducedCost;
			}
			
		}
		if(greatestCost>0) {
			System.out.println("reducedCost is "+greatestCost);
			return greatestCost;
		}
		else {
			System.out.println("reducedCost is " + greatestCost);
			return greatestCost;
		}
		
	}
	
	/** Called after we've found entering and leaving variable, we update the basicVars array to more 
	 * easily construct our new BFS 
	 * @param leaving
	 * @param entering
	 */
	private void updateBasicVars(int leaving, int entering) {
		for(int i=0;i<numRows;i++) {
			if(basicVars[i] == leaving) {
				basicVars[i] = entering;
			}
		}
		System.out.println("Basic vars "+Arrays.toString(basicVars));
	}
	/** Called after we've found entering and leaving variable, we update the nonBasicVars array 
	 * to more easily find reduced cost in next iteration
	 * @param leaving
	 * @param entering
	 */
	private void updateNonBasicVars(int leaving, int entering) {
		for(int i=0;i<numCols-numRows;i++) {
			if(nonBasicVars[i] == entering) {
				nonBasicVars[i] = leaving;
			}
		}
		System.out.println("Non basic vars " +Arrays.toString(nonBasicVars));
	}
	/**Finds leaving variable after we've found entering variable 
	 * Calculates the most restrictive constraint and chooses the variable corresponding with that constraint
	 */
	public void findLeavingVariable() {
		SimpleMatrix m1 = BInverse.mult(b);
		SimpleMatrix m2 = BInverse.mult(colMap.get(enteringVar));
		//iterate through m1 and m2 col by col and get divide them
		//we want to find the smallest number
		double mostRestrictive = Integer.MAX_VALUE;
		for(int i=0;i<numRows;i++) {
			double numerator = m1.get(i, 0);
			double denominator = m2.get(i,0);
			double quotient = numerator/denominator;
			if(/**quotient == 0 &&**/ denominator<0 && numerator >= 0) { //this actually has no upper bound
				//therefore we don't consider this case
			}
			else if(quotient<mostRestrictive) {
				
				leavingVar = basicVars[i];
				mostRestrictive = quotient;
			}
		}
		
		System.out.println("leaving vector we found is "+ colMap.get(leavingVar));
	}
	public static void main(String args[]){
		double[][] bArr = {{100},{120},{45},{30}};
		double[][] arr = {{1,0,0,0,3,2}, {0,1,0,0,2,4},{0,0,1,0,1,1}, {0,0,0,1,0,1}};
		HashMap<Integer, Integer> objectFuncCostMap = new HashMap<Integer,Integer>();
		objectFuncCostMap.put(5,200);
		objectFuncCostMap.put(6,300);
		SimplexSolver ss= new SimplexSolver(arr, bArr, objectFuncCostMap);
		
	}
}
