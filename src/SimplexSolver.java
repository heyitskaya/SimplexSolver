
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
	static List<Set<Integer>> possibleCombinations;
	int currCombo=0;
	SimpleMatrix B;

	/**use double for more precision**/
	/**for getting all subsets of size numRows to create our BFS**/
	
	public SimplexSolver(double[][] arr) {
		 internalMatrix= new DMatrixRMaj(arr);
		 this.numRows = arr.length;
		 this.numCols = arr[0].length;
		 System.out.println("numRows " +numRows);
		 System.out.println("numCols " +numRows);
		 internalMatrix.print();
		 SimpleMatrix simple = SimpleMatrix.wrap(internalMatrix);
		 SimpleMatrix col1 = simple.extractVector(false, 0);
		 SimpleMatrix col2 = simple.extractVector(false, 1);
		 SimpleMatrix col3 = simple.extractVector(false, 2);
		 colMap.put(1,col1);
		 colMap.put(2, col2);
		 colMap.put(3, col3);
		 col1.print();
		 col2.print();
		 col3.print();
		 //public static List<Set<Integer>> getSubsets(ArrayList<Integer> superSet, int k) {
		 ArrayList<Integer> superSet = new ArrayList<Integer>();
		 superSet.add(1);
		 superSet.add(2);
		 superSet.add(3);
		 getSubsets(superSet,numRows);
		 findBFS();
		 
		
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
	
	public int findDeterminant() {
		DMatrixRMaj A = new DMatrixRMaj(2,3,true,1.1,2.34,3.35436,4345,59505,0.00001234);
		A.print();
	    System.out.println();
	    A.print("%e");
	    System.out.println();
	    A.print("%10.2f");
	    //System.out.println("determinant " + determinant(A));
	    System.out.println("Determinant = "+ CommonOps_DDRM.det(A));
	    return -1;
	}
	
	
	public void findBFS() {
		boolean found = false;
		while(!found) {
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
				return;
			}
			else {
				currCombo++;
			}
	}
		
	
		//m1.insertIntoThis(0,1,colMap.get(2));
		
		
		
		
	}
	
	public static void main(String args[]){
		ArrayList<Integer> superSet = new ArrayList<>();
		superSet.add(1);
		superSet.add(2);
		superSet.add(3);
		superSet.add(4);
		System.out.println(getSubsets(superSet,2));
		double[][] arr = {{1.0,2.0,3.0}, {4.0,5.0,6.0}};
		SimplexSolver ss= new SimplexSolver(arr);
		//ss.findDeterminant();
		
	}
}
