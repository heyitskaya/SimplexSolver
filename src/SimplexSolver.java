
import java.util. *;
//import org.ejml.simple.SimpleBase<SimpleMatrix>;
import org.ejml.simple.SimpleMatrix;
import org.ejml.data.DMatrixRMaj;
import org.ejml.data.DMatrixD1;
import org.ejml.data.DMatrix1Row;
import org.ejml.simple.SimpleBase;
//import org.ejml.simple.SimpleBase<SimpleMatrix>;
import org.ejml.UtilEjml;
import org.ejml.simple.ops.SimpleOperations_DDRM;
import org.ejml.dense.row.CommonOps_DDRM;

import org.ejml.ops.*;

import org.ejml.data.DMatrix3;
import org.ejml.data.DMatrix3x3;
import org.ejml.data.DMatrixRMaj;
import org.ejml.dense.fixed.CommonOps_DDF3;
import org.ejml.ops.ConvertDMatrixStruct;
import org.ejml.simple.SimpleMatrix;
import org.ejml.data.DMatrix;
import org.ejml.simple.SimpleBase;
//import org.ejml.simple.SimpleBase<T>;
import java.util.Iterator;

public class SimplexSolver{
	/**Using EJML**/
	int numRows=0;
	int numCols=0;
	DMatrixRMaj internalMatrix;
	HashMap<Integer,SimpleMatrix> colMap = new HashMap<Integer,SimpleMatrix>();
	HashMap<Integer, Integer> objectFuncCostMap = new HashMap<Integer,Integer>();
	
	static List<Set<Integer>> possibleCombinations;
	int currCombo=0;
	SimpleMatrix B;
	SimpleMatrix BInverse;
	int[] basicVars;
	int[] nonBasicVars;
	int enteringVar;
	int leavingVar;
	
	double[][] bArr = {{100},{120},{45},{30}};
	SimpleMatrix b = new SimpleMatrix(bArr);

	/**use double for more precision**/
	/**for getting all subsets of size numRows to create our BFS**/
	
	public SimplexSolver(double[][] arr) {
		 internalMatrix= new DMatrixRMaj(arr);
		 this.numRows = arr.length;
		 this.numCols = arr[0].length;
		 nonBasicVars = new int[numCols-numRows];
		 basicVars = new int[numRows];
		 System.out.println("numRows " +numRows);
		 System.out.println("numCols " +numRows);
		 internalMatrix.print();
		 SimpleMatrix simple = SimpleMatrix.wrap(internalMatrix);
		 for(int i=1;i<=numCols;i++) {
			 colMap.put(i, simple.extractVector(false, i-1));
		 }
		 //public static List<Set<Integer>> getSubsets(ArrayList<Integer> superSet, int k) {
		 ArrayList<Integer> superSet = new ArrayList<Integer>();
		 for(int i=1;i<=numCols;i++) {
			 superSet.add(i);
		 }
		 getSubsets(superSet,numRows);
		 objectFuncCostMap.put(5,200);
		 objectFuncCostMap.put(6,300);
		 findBFS();
		 findEnteringVariable();
		 findLeavingVariable();
		 
		
	}
	private static void getSubsets(ArrayList<Integer> superSet, int k, int index, Set<Integer> current,List<Set<Integer>> solution) {
	    //successful stop clause
	    if (current.size() == k) {
	        solution.add(new HashSet<>(current));
	        return;
	    }
	    //unseccessful stop clause
	    if (index == superSet.size()) return;
	    Integer x = superSet.get(index);
	    current.add(x);
	    //"guess" x is in the subset
	    getSubsets(superSet, k, index+1, current, solution);
	    current.remove(x);
	    //"guess" x is not in the subset
	    getSubsets(superSet, k, index+1, current, solution);
	}

	public static List<Set<Integer>> getSubsets(ArrayList<Integer> superSet, int k) {
	    possibleCombinations = new ArrayList<>();
	    getSubsets(superSet, k, 0, new HashSet<Integer>(), possibleCombinations);
	    return possibleCombinations;
	}
	
	/**https://gist.github.com/Cellane/398372/23a3e321daa52d4c6b68795aae093bf773ce2940**/
	public static double matrixDeterminant (double[][] matrix) {
		double temporary[][];
		double result = 0;

		if (matrix.length == 1) {
			result = matrix[0][0];
			return (result);
		}

		if (matrix.length == 2) {
			result = ((matrix[0][0] * matrix[1][1]) - (matrix[0][1] * matrix[1][0]));
			return (result);
		}

		for (int i = 0; i < matrix[0].length; i++) {
			temporary = new double[matrix.length - 1][matrix[0].length - 1];

			for (int j = 1; j < matrix.length; j++) {
				for (int k = 0; k < matrix[0].length; k++) {
					if (k < i) {
						temporary[j - 1][k] = matrix[j][k];
					} else if (k > i) {
						temporary[j - 1][k - 1] = matrix[j][k];
					}
				}
			}

			result += matrix[0][i] * Math.pow (-1, (double) i) * matrixDeterminant (temporary);
		}
		return (result);
	}
	public void test() {
		double[][] dMatrix= {{1,0,9},{2,7,6}};
		SimpleMatrix sm = new SimpleMatrix(dMatrix);
	}
	
/**	public int findDeterminant() {
		DMatrixRMaj A = new DMatrixRMaj(2,3,true,1.1,2.34,3.35436,4345,59505,0.00001234);
		A.print();
	    System.out.println();
	    A.print("%e");
	    System.out.println();
	    A.print("%10.2f");
	    //System.out.println("determinant " + determinant(A));
	    System.out.println("Determinant = "+ CommonOps_DDRM.det(A));
	    return -1;
	} **/
	
	
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
			if(CommonOps_DDRM.det(m1.getMatrix())==1) {
				//then they are linearly independent which means we've
				//currently found our BFS
				B=m1;
				found = true;
				currCombo++;
				B.print(); //beautiful!
				ArrayList<Integer> nonBasic = new ArrayList<Integer>();
				
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
		
				for(int i=0;i<basicVars.length;i++) {
					System.out.println(basicVars[i]);
				} 
				return;
			}
			else {
				currCombo++;
			}
		}
		
	}
	
	public void findEnteringVariable() {
		System.out.println("In findEnteringVariable");
		//first compute reduced cost of each non basic variable and then pick the 
		//biggest one 
		int greatestCost = Integer.MIN_VALUE;
		for(int i=0;i<nonBasicVars.length;i++) {
			int currVar = nonBasicVars[i];
			int objectiveFunctionCost;
			if(objectFuncCostMap.get(currVar)==null) {
				objectiveFunctionCost = 0;
			}
			else {
				objectiveFunctionCost = objectFuncCostMap.get(currVar);
			}
			//pick up here
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
			System.out.println("cbTranspose looks like this ");
			cbTranspose.print(); //looks legit
			//create BInverse
			BInverse = B.invert();
			System.out.println("BInverse looks like");
			BInverse.print(); //also looks legit 
		//	int matrixProduct = 
		//	int reducedCost = objectiveFunctionCost
			SimpleMatrix tempMatrix1 = cbTranspose.mult(BInverse);
			SimpleMatrix tempMatrix2 = tempMatrix1.mult(colMap.get(currVar));
			tempMatrix2.print();
			double currReducedCost = objectiveFunctionCost-tempMatrix2.get(0,0);
			if(currReducedCost>greatestCost) {
				enteringVar = currVar;
			}
			
		}
		System.out.println("enteringVar we found was "+ enteringVar);
	}
	
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
			if(quotient<mostRestrictive) {
				leavingVar = basicVars[i];
			}
		}
		System.out.println("leaving var we found is "+ leavingVar);
	}
	public static void main(String args[]){
		ArrayList<Integer> superSet = new ArrayList<>();
		superSet.add(1);
		superSet.add(2);
		superSet.add(3);
		superSet.add(4);
		System.out.println(getSubsets(superSet,2));
		double[][] arr = {{1,0,0,0,3,2}, {0,1,0,0,2,4},{0,0,1,0,1,1}, {0,0,0,1,0,1}};
		SimplexSolver ss= new SimplexSolver(arr);
		//ss.findDeterminant();
		
	}
}
