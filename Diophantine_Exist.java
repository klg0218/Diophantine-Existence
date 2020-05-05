import java.util.Scanner;

public class Diophantine_Exist {
	
	/*
	 This method determines whether the Diophantine equation 
	 ax + by = c has an integer solution with x,y>=0
	 */
	public static boolean Diophantine3(int a, int b, int c) {
		boolean sol = false; 
		
		for(int i=0; i<=c; i++)
		{
			for(int j=0; j<=c; j++)
			{
				if((a*i + b*j)==c)
				{
					sol = true;
				}
			}
		}
		
		return sol; 
	}
	
	/*
	 This method returns all integer solutions of the 
	 Diophantine equation (x^a)(y^b) = M. 
	 */
	public static int[][] Diophantine2(int a, int b, int M) {
		int count = 0; 
		
		for(int i=(-1)*M; i<=M; i++)
		{
			for(int j=(-1)*M; j<=M; j++)
			{
				int pow1 = (int) Math.pow((double) i, (double) a); 
				int pow2 = (int) Math.pow((double) j, (double) b); 
				
				if((pow1*pow2)==M)
				{
					count++;
				}
			}
		}
		
		int[][] output = new int[count][2]; 
		
		count = 0; 
		
		for(int i=(-1)*M; i<=M; i++)
		{
			for(int j=(-1)*M; j<=M; j++)
			{
				int pow1 = (int) Math.pow((double) i, (double) a); 
				int pow2 = (int) Math.pow((double) j, (double) b); 
				
				if((pow1*pow2)==M)
				{
					output[count][0] = i; 
					output[count][1] = j; 
					
					count++;
				}
			}
		}
		
		return output; 
	}
	
	public static boolean prime_check(int p) {
		boolean prime = true;
		
		for(int i=2; i<p; i++)
		{
			if(p%i == 0)
			{
				prime = false;
				break;
			}
		}
		
		return prime; 
	}
	
	public static int[][] prime_fac(int m) {
		int M = Math.abs(m); 
		
		int factor_count = 0; 
		for(int i=2; i<=M; i++)
		{
			if((prime_check(i)==true)&&(M%i==0))
			{
				factor_count++; 
			}
		}
		
		int[][] output = new int[2][factor_count];
		int exp = 0;
		int count = 0; 
		
		for(int i=2; i<=M; i++)
		{
			if((prime_check(i)==true)&&(M%i==0))
			{
				exp = 0; 
				
				int M1 = M;
				
				while(M1%i==0) {
					M1 = M1/i; 
					exp++; 
				}
				
				output[0][count] = i;
				output[1][count] = exp; 
				
				count++;
			}
		}
		
		return output; 
	}
	
	/*
	 This method determines whether an integer solution to the Diophantine
	 equation (x^a)(y^b) = M exists. 
	 */
	public static boolean Diophantine1(int a, int b, int M)
	{
		boolean exists=true; 
		
		if((M==0)||(M==1))
		{
			return true; 
		}
		
		if(M<0)
		{
			if((a%2==0)&&(b%2==0))
			{
				return false;  
			}
		}
		
		int[][] pf1 = prime_fac(M); 
		int[] pf = pf1[1];
		int l = pf.length; 
		
		for(int i=0; i<l; i++)
		{
			if(Diophantine3(a,b,pf[i])==false)
			{
				exists = false;
				break; 
			}
		}
		
		return exists; 
	}

	public static int[] double_convert(double x) {
		
		double c = 1;
		
		double y = c*x - Math.round(c*x); 
		
		while(Math.abs(y)>0.000001) {
			c++; 
			
			y = c*x - Math.round(c*x); 
		}
		
		int[] output = new int[2];
		output[1] = (int) Math.round(c); 
		output[0] = (int) Math.round(c*x);
		
		return output; 
	}
	
	/*
	 This method uses the rational root theorem to determine whether a univariate  
	 polynomial with integer coefficients has any integer roots. 
	 */
	public static boolean zero_check(double[] input) {
		double A0 = input[0]; 
		double a0 = Math.abs(A0); 
		int n = input.length; 
		
		boolean zero = false; 
		
		if(a0==0)
		{
			zero = true;
		}
		else {
			
			for(double i = ((-1)*a0); i <= a0; i++)
			{
				if(i!=0)
				{
					
				int I = (int) i;
				int A = (int) a0; 
				
				if(A%I==0)
				{
					double z = 0;  
					
					for(int j=0; j<n; j++)
					{
						double J = (double) j;
						z += input[j]*(Math.pow(I, J)); 
					}
					
					if(z==0)
					{
						zero = true;
						break;
					}
				}
			}
			}
		}
		
		return zero; 
	}
	
	public static double[][] metric_factor(double[][] poly) {
		int l = poly.length; 
		
		double[][] output = new double[l][l];
		for(int i=0; i<l; i++)
		{
			for(int j=0; j<l; j++)
			{
				output[i][j] = 0; 
			}
		}
		
		int x = 0;
		int y = 0; 
		
		for(int i=0; i<l; i++)
		{
			for(int j=0; j<l; j++)
			{
				if(poly[i][j]!=0)
				{
					x = i;
					y = j;
				}
			}
		}
		
		double a = poly[x][y]; 
		double N = poly[0][0];
		
		if(y!=0)
		{
			for(int i=0; i<l; i++)
			{
				for(int j=0; j<l; j++)
				{
					if((i!=x)||(j!=y))
					{
						double J = (double) j+1; 
						output[i][j] = poly[i][j]/(a*J); 
					}
				}
			}
			
			output[0][0] = N/a; 
			output[x][y] = 1/(y+1); 
			
			double K = 0; 
			
			for(int i=0; i<=x; i++)
			{
				K += output[i][0];
			}
			
			K *= -1; 
			
			output[x][0] += K; 
		}
		
		if(y==0)
		{
			for(int i=0; i<=x; i++)
			{
				double I = (double) i+1;
				output[i][0] = poly[i][0]/(a*I); 
			}
			
			double K = 0; 
			
			for(int i=0; i<=x; i++)
			{
				K += output[i][0];
			}
			
			K *= -1; 
			
			output[0][0] += K; 
		}
		
		return output;
	}
	
	public static double[][] reduce(double[][] poly) {
		int l = poly.length; 
		
		int x = 0;
		int y = 0; 
		
		for(int i=0; i<l; i++)
		{
			for(int j=0; j<l; j++)
			{
				if(poly[i][j]!=0)
				{
					x = i;
					y = j;
				}
			}
		}
		
		double[][] metric = metric_factor(poly); 
		
		double[][] output = new double[l][l]; 
		for(int i=0; i<l; i++)
		{
			for(int j=0; j<l; j++)
			{
				output[i][j] = 0;
			}
		}
		
		if(metric[x][y]!=0)
		{
		double c = poly[x][y]/metric[x][y]; 
		
		for(int i=0; i<=x; i++)
		{
			for(int j=0; j<=y; j++)
			{
				output[i][j] = poly[i][j] - c*metric[i][j]; 
			}
		}
		output[x][y] = 0; 
		
		return output; 
		}
		else {
			return metric; 
		}
	}
	
	public static void main(String[] args) {
		
		Scanner input = new Scanner(System.in); 
		System.out.println("Enter the maximum exponent which appears in the polynomial:");
		int m = input.nextInt();
		
		System.out.println(" "); 
		
		double[][] poly = new double[m+1][m+1]; 
		for(int i=0; i < (m+1); i++)
		{
			for(int j=0; j < (m+1); j++)
			{
				System.out.println("Input the coefficent of the x^" + i + " y^" + j + " term. (Integers only.)"); 
				poly[i][j] = input.nextDouble(); 
			}
		}
		input.close(); 
		
		double[][] temp1 = new double[m+1][m+1];
		for(int i=0; i < (m+1); i++)
		{
			for(int j=0; j < (m+1); j++)
			{
				temp1[i][j] = poly[i][j]; 
			}
		}
		
		System.out.println(" "); 
		
		if(m==0)
		{
			if(poly[0][0]==0)
			{
				System.out.println("The Diophantine equation has an integer solution.");
			}
			else {
				System.out.println("The Diophantine equation does not have an integer solution.");
			}
			
			System.exit(0);
		}
		
		if(poly[0][0]==0)
		{
			System.out.println("The Diophantine equation has an integer solution.");
			System.exit(0);
		}
		
		double[] xpoly = new double[m+1];
		for(int i=0; i < (m+1); i++)
		{
			xpoly[i] = poly[i][0]; 
		}
		
		double[] ypoly = new double[m+1]; 
		for(int i=0; i < (m+1); i++)
		{
			ypoly[i] = poly[0][i];
		}
		
		if((zero_check(xpoly)==true)||(zero_check(ypoly)==true))
		{
			System.out.println("The Diophantine equation has an integer solution."); 
			System.exit(0);
		}
		
		int x = 0;
		int y = 0; 
		
		for(int i=0; i < (m+1); i++)
		{
			for(int j=0; j < (m+1); j++)
			{
				if(poly[i][j]!=0)
				{
					x = i;
					y = j;
				}
			}
		}
		
		double[][] temp = new double[m+1][m+1];
		
		while((x!=0)||(y!=0)) {
			
			for(int i=0; i < (m+1); i++)
			{
				for(int j=0; j < (m+1); j++)
				{
					temp[i][j] = poly[i][j]; 
				}
			}
			
			poly = reduce(temp); 
			
			for(int i=0; i < (m+1); i++)
			{
				for(int j=0; j < (m+1); j++)
				{
					if(poly[i][j]!=0)
					{
						x = i;
						y = j;
					}
				}
			}
		}
		
		for(int i=0; i < (m+1); i++)
		{
			for(int j=0; j < (m+1); j++)
			{
				poly[i][j] = temp[i][j]; 
			}
		}
		
		int[][][] poly1 = new int[m+1][m+1][2]; 
		
		for(int i=0; i < (m+1); i++)
		{
			for(int j=0; j < (m+1); j++)
			{
				poly1[i][j] = double_convert(poly[i][j]); 
			}
		}
		
		int coeff_count = 0;
		for(int i=0; i < (m+1); i++)
		{
			for(int j=0; j < (m+1); j++)
			{
				if(poly1[i][j][0] != 0)
				{
					coeff_count++; 
				}
			}
		}
		
		if(coeff_count==2)
		{
			int x1 = 0;
			int y1 = 0;
			
			for(int i=0; i < (m+1); i++)
			{
				for(int j=0; j < (m+1); j++)
				{
					if((i!=0)||(j!=0))
					{
						if(poly1[i][j][0] != 0)
						{
							x1 = i;
							y1 = j; 
						}
					}
				}
			}
			
			int[] q = new int[2]; 
			q[0] = (-1)*poly1[0][0][0]*poly1[x1][y1][1]; 
			q[1] = poly1[0][0][1]*poly1[x1][y1][0]; 
			
			
			if(q[0]%q[1] != 0)
			{
				System.out.println("The Diophantine equation does not have an integer solution.");
				System.exit(0); 
			}
			
			int M = q[0]/q[1]; 
			
			if(x1==0)
			{
				double M1 = (double) M; 
				double Y1 = (double) y1; 
				
				double root = Math.pow(M1, 1/Y1); 
				double c = root - Math.round(root); 
				
				if(c!=0)
				{
					System.out.println("The Diophantine equation does not have an integer solution.");
					System.exit(0);
				}
				
				double[] poly_root = new double[m+1]; 
				
				for(int i=0; i < (m+1); i++)
				{
					poly_root[i] = 0;
					
					for(int j=0; j < (m+1); j++)
					{
						poly_root[i] += temp1[i][j]*Math.pow(root, j); 
					}
				}
				
				if(zero_check(poly_root)==true)
				{
					System.out.println("The Diophantine equation has an integer solution."); 
					System.exit(0);
				}
				
				if(y1%2==0)
				{
					root *= -1; 
					
					for(int i=0; i < (m+1); i++)
					{
						poly_root[i] = 0;
						
						for(int j=0; j < (m+1); j++)
						{
							poly_root[i] += temp1[i][j]*Math.pow(root, j); 
						}
					}
		
					if(zero_check(poly_root)==true)
					{
						System.out.println("The Diophantine equation has an integer solution.");
						System.exit(0);
					}
					else {
						System.out.println("The Diophantine equation does not have an integer solution."); 
						System.exit(0);
					}
				}
				else {
					System.out.println("The Diophantine equation does not have an integer solution.");
					System.exit(0);
				}
			}
			
			if(y1==0)
			{
				double M1 = (double) M; 
				double X1 = (double) x1; 
				
				double root = Math.pow(M1, 1/X1); 
				double c = root - Math.round(root); 
				
				if(c!=0)
				{
					System.out.println("The Diophantine equation does not have an integer solution.");
					System.exit(0);
				}
				
				double[] poly_root = new double[m+1]; 
				
				for(int i=0; i < (m+1); i++)
				{
					poly_root[i] = 0; 
					
					for(int j=0; j < (m+1); j++)
					{
						poly_root[i] += temp1[j][i]*Math.pow(root, j); 
					}
				}
				
				if(zero_check(poly_root)==true)
				{
					System.out.println("The Diophantine equation has an integer solution."); 
					System.exit(0); 
				}
				
				if(x1%2==0)
				{
					root *= -1; 
					
					for(int i=0; i < (m+1); i++)
					{
						poly_root[i] = 0; 
						
						for(int j=0; j < (m+1); j++)
						{
							poly_root[i] += temp1[j][i]*Math.pow(root, j); 
						}
					}
					
					if(zero_check(poly_root)==true)
					{
						System.out.println("The Diophantine equation has an integer solution.");
						System.exit(0);
					}
					else {
						System.out.println("The Diophantine equation does not have an integer solution.");
						System.exit(0);
					}
				}
				else {
					System.out.println("The Diophantine equation does not have an integer solution.");
					System.exit(0);
				}
			}
			
			if((x1!=0)&&(y1!=0))
			{
				if(Diophantine1(x1,y1,M)==true)
				{
					int[][] D = Diophantine2(x1,y1,M); 
					int l = D.length; 
					
					for(int i=0; i<l; i++)
					{
						double x2 = (double) D[i][0];
						double y2 = (double) D[i][1]; 
						
						double eval = 0; 
						for(int j=0; j < (m+1); j++)
						{
							for(int k=0; k < (m+1); k++)
							{
								eval += temp1[j][k]*Math.pow(x2, (double) j)*Math.pow(y2, (double) k); 
							}
							
							if(eval==0)
							{
								System.out.println("The Diophantine equation has an integer solution.");
								System.exit(0);
							}
						}
					}
					
					System.out.println("The Diophantine equation does not have an integer solution.");
					System.exit(0);
				}
				else {
					System.out.println("The Diophantine equation does not have an integer solution.");
					System.exit(0);
				}
			}
		}
		
		if(coeff_count>2)
		{
			double[][] metric = metric_factor(poly); 
			
			int[][][] metric1 = new int[m+1][m+1][2];
			for(int i=0; i < (m+1); i++)
			{
				for(int j=0; j < (m+1); j++)
				{
					metric1[i][j] = double_convert(metric[i][j]); 
				}
			}
			
			int[][] x_poly = new int[m+1][2]; 
			for(int i=0; i < (m+1); i++)
			{
				for(int j=0; j<2; j++)
				{
					x_poly[i][j] = metric1[i][0][j]; 
				}
			}
			
			int denom = 1;
			
			for(int i=0; i < (m+1); i++)
			{
				denom *= x_poly[i][1]; 
			}
			
			int[] x_poly1 = new int[m+1];
			for(int i=0; i < (m+1); i++)
			{
				x_poly1[i] = denom*x_poly[i][0]/x_poly[i][1];
			}
			
			int xa0 = Math.abs(x_poly1[0]); 
			
			for(int i = (-1)*xa0; i<=xa0; i++)
			{
				if(i!=0)
				{
					if(xa0%i==0)
					{
						int z = 0; 
						
						for(int j=0; j < (m+1); j++)
						{
							z += x_poly1[j]*((int) Math.pow((double) i, (double) j));
						}
						
						if(z==0)
						{
							double[] reduce_poly = new double[m+1];
							
							for(int k=0; k < (m+1); k++)
							{
								reduce_poly[k] = 0; 
								
								for(int s=0; s < (m+1); s++)
								{
									reduce_poly[k] += temp1[s][k]*Math.pow((double) i, (double) s);
								}
							}
							
							if(zero_check(reduce_poly)==true)
							{
								System.out.println("The Diophantine equation has an integer solution."); 
								System.exit(0); 
							}
						}
					}
				}
			}
			
			int[][] y_poly = new int[m+1][2]; 
			for(int i=0; i < (m+1); i++)
			{
				for(int j=0; j<2; j++)
				{
					y_poly[i][j] = metric1[0][i][j]; 
				}
			}
			
			int ydenom = 1;
			
			for(int i=0; i < (m+1); i++)
			{
				ydenom *= y_poly[i][1]; 
			}
			
			int[] y_poly1 = new int[m+1];
			for(int i=0; i < (m+1); i++)
			{
				y_poly1[i] = ydenom*y_poly[i][0]/y_poly[i][1];
			}
			
			int ya0 = Math.abs(y_poly1[0]); 
			
			for(int i = (-1)*ya0; i<=ya0; i++)
			{
				if(i!=0)
				{
					if(ya0%i==0)
					{
						int z = 0; 
						
						for(int j=0; j < (m+1); j++)
						{
							z += y_poly1[j]*((int) Math.pow((double) i, (double) j));
						}
						
						if(z==0)
						{
							double[] reduce_poly = new double[m+1];
							
							for(int k=0; k < (m+1); k++)
							{
								reduce_poly[k] = 0; 
								
								for(int s=0; s < (m+1); s++)
								{
									reduce_poly[k] += temp1[k][s]*Math.pow((double) i, (double) s);
								}
							}
							
							if(zero_check(reduce_poly)==true)
							{
								System.out.println("The Diophantine equation has an integer solution."); 
								System.exit(0); 
							}
						}
					}
				}
			}
			
			System.out.println("The Diophantine equation does not have an integer solution.");
		}
	}
}
