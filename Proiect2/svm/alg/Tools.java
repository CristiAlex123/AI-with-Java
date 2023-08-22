package alg;

public class Tools{

	public static float[] translate(float[] w, io.Vector[] V){
		int dim = w.length - 1;
		int N = V.length;
		float max = Float.MIN_VALUE, min;
		int imax = -1, imin;		
		for(int i = 0; i < N; i++){
			float d =  distFromHiperplanToVector(w, V[i]);
			if(d > max){
				max = d; imax = i;
			}
		}
		
		int y = V[imax].cl.Y;		
		float s = 0;
		for(int j = 0; j < dim; j++) s += w[j]*V[imax].X[j];
		w[dim] = -s;

		min = Float.MAX_VALUE; max = Float.MIN_VALUE;
		imax = -1; imin = -1;			

		for(int i = 0; i < N; i++){
			float d = distFromHiperplanToVector(w, V[i]);
			if(V[i].cl.Y == y){
				if(d > max){max = d; imax = i;}
			}else{
				if(d < min){min = d; imin = i;}
			}
		}
		
		float bmax=0, bmin=0;
		s = 0;
		for(int j = 0; j < dim; j++) s += w[j]*V[imax].X[j];
		bmax = s;
		s = 0;
		for(int j = 0; j < dim; j++) s += w[j]*V[imin].X[j];
		bmin = s;
		
		float[] b = new float[4];
		b[0] = bmin; 
		b[1] = bmax;
		b[2] = (bmin + bmax)/2;		
		
		w[dim] = -bmax;		
		b[3] = distFromHiperplanToVector(w, V[imin]); // the margin;
		
		return b;
	}	


	public static float[] translate(float[] w, io.Vector[] V, int y){
		int dim = w.length - 1;
		int N = V.length;
		float min = Float.MAX_VALUE;
		int imin = -1;			

		for(int i = 0; i < N; i++){	
			if(V[i].cl.Y == y){
				float d = distFromHiperplanToVector(w, V[i]);
				if(d < min){min = d; imin = i;}
			}
		}
		
		float b = 0;
		for(int j = 0; j < dim; j++) b += w[j]*V[imin].X[j];
		
		float[] B = new float[2];
		B[0] = b; 
		B[1] = (float)imin;
		return B;
	}	

	public static float[] translate_(float[] w, io.Vector[] V){
		float[] b = translate(w, V, 0);
		float b0 = b[0];
		int i0 = (int)b[1];

		b = translate(w, V, 1);
		float b1 = b[0];
		int i1 = (int)b[1];

		b = new float[4];
		b[0] = b0; 
		b[1] = b1;
		b[2] = (b0 + b1)/2;		
		
		w[w.length-1] = -b1;		
		b[3] = distFromHiperplanToVector(w, V[i0]); // the margin;		
		
		return b;
	}		

	
	public static float distFromHiperplanToVector(float[] w, io.Vector V){
		float dist = 0;
		float norm = 0;
		for(int j = 0; j < w.length-1; j++){
			dist += w[j]*V.X[j];
			norm += w[j]*w[j];
		}
		dist += w[w.length-1];
		dist = Math.abs(dist);
		norm = (float)Math.sqrt(norm);
		dist /= norm;
		return dist;
	}
	
	

}