import java.util. *;
public class SimplexSolver {
	/**Using EJML**/
	int numRows=0;
	int numCols=0;

	/**use double for more precision**/
	/**for getting all subsets of size numRows to create our BFS**/
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
	    List<Set<Integer>> res = new ArrayList<>();
	    getSubsets(superSet, k, 0, new HashSet<Integer>(), res);
	    return res;
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
	
	
	
	
	public static void main(String args[]){
		ArrayList<Integer> superSet = new ArrayList<>();
		superSet.add(1);
		superSet.add(2);
		superSet.add(3);
		superSet.add(4);
		System.out.println(getSubsets(superSet,2));
		
	}
}
