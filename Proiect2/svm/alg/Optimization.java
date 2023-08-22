package alg;

import io.*;
import java.util.Random;

public class Optimization{
	// aici comentariu //
	
	public static int Wolfe(int n, int m, double c[][],double d[], double a[][], double b[], int maxiterations, double sol[]){
		int mplusn,m2n2,m2n3,m2n3plus1,big,unrestricted; 
		int i,j,iterations,temp,column,row=0; 
		int index[] = new int[n+m+1]; 
		double total,candidate,dividend; 
		double tableau[][] = new double[n+m+1][2*m+3*n+2]; 
		double price[] = new double[2*m+3*n+2]; 
		double change[] = new double[2*m+3*n+2]; 
		double net[] = new double[2*m+3*n+2]; 
		double gain[] = new double[2*m+3*n+2]; 
		double fraction[] = new double[n+m+1]; 
	
	
		//initializari
		for (i=0; i<=n; i++) 
			sol[i] = 0.; 
			
		for (i=0; i<=n; i++) 
			for (j=0; j<=n; j++) 
				if (i != j) c[i][j] /= 2.; 
		big = Integer.MAX_VALUE; 
		mplusn = m + n; 
		m2n2 = m + m + n + n; 
		m2n3 = m2n2 + n; 
		m2n3plus1 = m2n3 + 1; 
		for (i=1; i<=mplusn; i++) 
			for (j=1; j<=m2n3plus1; j++) 
				tableau[i][j] = 0.; 
		for (i=1; i<=m; i++) 
			tableau[i][1] = b[i]; 
		for (i=m+1; i<=mplusn; i++) 
			tableau[i][1] = -d[i-m]; 
		for (i=1; i<=m; i++) 
			for (j=1; j<=n; j++) 
				tableau[i][j+1] = a[i][j]; 
		for (i=1; i<=n; i++) 
			for (j=1; j<=n; j++) 
				tableau[i+m][j+1] = 2. * c[i][j]; 
		for (i=m+1; i<=mplusn; i++) 
			for (j=n+2; j<=mplusn+1; j++) 
				tableau[i][j] = a[j-n-1][i-m]; 
		for (i=1; i<=mplusn; i++) { 
			temp = i + mplusn + n + 1; 
			for (j=m2n2+2; j<=m2n3plus1; j++) 
				if (j == temp) tableau[i][j] = 1.; 
		} 
		for (i=m+1; i<=mplusn; i++) { 
			temp = i - m + mplusn + 1; 
			for (j=mplusn+2; j<=m2n3plus1; j++) 
				if (j == temp) tableau[i][j] = -1.; 
		} 
		for (j=1; j<=m2n3; j++) 
			price[j] = 0.; 
		for (i=1; i<=m; i++) 
			price[n+1+i] = tableau[i][1]; 
		for (j=m2n2+2; j<=m2n3plus1; j++) 
			price[j] = big - 1; 
		for (i=1; i<=mplusn; i++) 
			index[i] = m2n3 - mplusn + i; 
		iterations = 0;
		while (true) { 
			// iteration start 
			iterations++; 
			for (j=1; j<=m2n3plus1; j++) 
				gain[j] = 0.; 
			for (j=1; j<=m2n3plus1; j++) { 
				total = 0.; 
				for (i=1; i<=mplusn; i++) 
				total += price[index[i]+1] * tableau[i][j]; 
				gain[j] = total; 
				change[j] = price[j] - gain[j]; 
			} 
			// search for the pivot element 
			column = 0; 
			candidate = 0.; 
			// get the variable with largest gain 
			for (i=2; i<=m2n3plus1; i++) 
			if (change[i] < candidate) { 
				candidate = change[i]; 
				column = i; 
			} 
			if (column <= 0) break; 
			unrestricted = 0; 
			for (i=1; i<=mplusn; i++) { 
				if (tableau[i][column] > 0) 
					fraction[i] = tableau[i][1] / tableau[i][column]; 
				else { 
					unrestricted++; 
					if (unrestricted == mplusn) 
						// objective function is unbounded 
						return 1; 
					else 
						fraction[i] = Double.MAX_VALUE; 
				} 
			} 
			// remove limiting variable 
			for (i=1; i<=mplusn; i++) 
			if (fraction[i] >= 0) { 
				if (fraction[i] > big) fraction[i] = big; 
				candidate = fraction[i]; 
				row = i; 
				break; 
			} 
			for (i=1; i<=mplusn; i++) 
				if (candidate > fraction[i]) { 
					candidate = fraction[i]; 
					row = i; 
				}
			// perform pivoting and introduce new variable 
			dividend = tableau[row][column]; 
			for (j=1; j<=m2n3plus1; j++) 
				tableau[row][j] /= dividend; 
			for (i=1; i<=mplusn; i++) 
				if (i != row) { 
					for (j=1; j<=m2n3plus1; j++) 
						net[j] = tableau[row][j] * tableau[i][column] / tableau[row][column]; 
					for (j=1; j<=m2n3plus1; j++) 
					tableau[i][j] -= net[j]; 
				}
			price[row] = price[column]; 
			index[row] = column - 1; 
			// recompute the price 
			for (j=1; j<=m2n2+1; j++) 
				price[j] = 0.; 
			for (i=1; i<=mplusn; i++) { 
				if (index[i] <= mplusn) 
					temp = index[i] + mplusn + 1; 
				else { 
					if (index[i] > m2n2) continue; 
						temp = index[i] - (mplusn - 1);
				} 
				price[temp] = tableau[i][1]; 
			} 
			if (iterations >= maxiterations) 
			// maximum number of iterations exceeded 
			return 2; 
		} 
		// return the optimal solution 
		total = 0.; 
		for (i=1; i<=mplusn; i++) 
			if (index[i] <= n) total += d[index[i]] * tableau[i][1]; 
		sol[0] = total; 
		total =0.; 
		for (i=1; i<=mplusn; i++) 
			for (j=1; j<=mplusn; j++) { 
				if (index[i] > n) continue; 
				if (index[j] > n) continue; 
				total += c[index[i]][index[j]] * tableau[i][1] * tableau[j][1]; 
			} 
		sol[0] += total; 
		for (i=1; i<=mplusn; i++) 
			if ((tableau[i][1] != 0) && (index[i] <= n)) 
				sol[index[i]] = tableau[i][1]; 
		return 0; 
	} 	
}
