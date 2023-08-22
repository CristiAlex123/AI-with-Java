package alg;

import svm.SVM;

public class MaxMarginRandom extends Algorithm{
	int R = 1000;
	
	public MaxMarginRandom(SVM svm){
		super(svm);
		if(svm.ind.V != null){
			name = "MaxMarginRandom";
			svm.outd.algorithm = name;
			svm.outd.showInputData();
		}
	}
	
	public void run(){
		t = System.currentTimeMillis();
		boolean flag = false;
		
		float margin = Float.MIN_VALUE;
		svm.design.calculates = false;
		
		
		float[] w = new float[dim+1];
		
		float[] M0 = new float[dim];
		float[] M1 = new float[dim];
		int k0 = 0, k1 = 0;		
		for(int i = 0; i < N; i++)
			if(svm.ind.V[i].cl.Y == 0) k0++; else k1++;
		for(int j = 0; j < dim; j++)
			for(int i = 0; i < N; i++)
				if(svm.ind.V[i].cl.Y == 0) 
					M0[j] += svm.ind.V[i].X[j];
				else 
					M1[j] += svm.ind.V[i].X[j];
		for(int j = 0; j < dim; j++){
			M0[j] /= k0;
			M1[j] /= k1;
		}
		float[] X0 = new float[dim];
		for(int j = 0; j < dim; j++){
			X0[j] = (M0[j] + M1[j])/2;
			w[j] =  M1[j] - M0[j];
			w[dim] -= w[j] * X0[j];
		}
		
		for(long r = 1;r <= R; r++){
			for(int j = 0; j < dim; j++) w[j] = -0.5f + (float)Math.random();
			for(int p = 1; p <= P; p++){
			boolean erori = false; 
			for(int i = 0; i < N; i++){
				float s = 0;
				for(int j = 0; j < dim; j++)
					s += w[j]*svm.ind.V[i].X[j];
				s += w[dim];
				int y = s < 0 ? 0 : 1;
				int e = svm.ind.V[i].cl.Y - y;
				if(e!=0){
					erori = true;
					for(int j = 0; j < dim; j++)
						w[j] += eta*e*svm.ind.V[i].X[j];
					w[dim] += eta*e;
				}
			}
			if(!erori){
				
				float[] b = translate(w, svm.ind.V);
				w[dim] = -b[2];
				float[] w0 = new float[dim+1];
				for(int j = 0; j < dim; j++) w0[j] = w[j];
				w0[dim] = -b[0];
				float[] w1 = new float[dim+1];
				for(int j = 0; j < dim; j++) w1[j] = w[j];
				w1[dim] = -b[1];
				
				if(margin < b[3]){
					margin = b[3];
					if(dim==2) svm.design.setPointsOfMaxLine(w,w0,w1);
					svm.outd.accuracy = getAccuracy(w);
					svm.outd.w = w;
					svm.outd.margin = b[3];
					svm.control.ta.append("Stage "+r+"/ Margin = "+ b[3] + "\n");
					String s = "";
					for(int j = 0; j < w.length; j++) s+= "w["+j+"] = "+ w[j] + "; ";
					svm.control.ta.append(s + "\n");
				}else
					svm.control.ta.append("Stage "+ r +"\n");
				
				break;
				}
		}
		
		svm.outd.stages_count = R;
		svm.outd.computing_time = System.currentTimeMillis() - t;
		svm.outd.showInputData();
		svm.outd.showOutputData();
		}
		
		svm.control.start.enable(true);
		
		
	}
	
	public static float[] translate(float[] w, io.Vector[] V){
		float[] b = translate(w,V,0);
		float b0 = b[0];
		int i0 = (int)b[1];
		
		b = translate(w, V, 1);
		float b1 = b[0];
		int i1 = (int)b[1];
		
		b = new float[4];
		b[0] = b0;
		b[1] = b1;
		b[2] = (b0 + b1)/2;
		
		w[w.length - 1] = -b1;
		b[3] = distFromHiperplanToVector(w, V[i0]); //the margin
		
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
				if(d < min){ min = d; imin = i;}
			}
		}
		
		float b = 0;
		for(int j = 0; j < dim; j++) b += w[j]*V[imin].X[j];
		
		float[] B = new float[2];
		B[0] = b;
		B[1] = (float)imin;
		return B;
	}
	
	public static float distFromHiperplanToVector(float[] w, io.Vector V){
		float dist = 0;
		float norm = 0;
		for(int j = 0; j < w.length - 1; j++)	
		{	dist += w[j] * V.X[j];
			norm += w[j] * w[j];	
		}
		dist += w[w.length-1];
		dist = Math.abs(dist);
		norm = (float)Math.sqrt(norm);
		dist /= norm;
		return dist;
	}
	
	
}