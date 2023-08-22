package alg;

import svm.SVM;
import io.*;
import java.util.*;

//*******Sequential Minimal Optimization (SMO) Algorithm*******
//***Reference***
//http://fourier.eng.hmc.edu/e176/lectures/ch9/node9.html

public class SMO extends Algorithm{
	SVM svm;
	
	public SMO(SVM svm){
		super(svm);
		this.svm=svm;
		if(svm.ind.V!=null){
			name="SMO";
			svm.outd.algorithm = name;
			svm.outd.max_stages_count =1;
			svm.outd.showInputData();
		}
	}
	
	public float dotProduct(float[] x,float[] y){
		float sum=0;
		for(int i=0;i<x.length;i++){
			sum+=x[i]*y[i];
		}
		return sum;
	}
	
	public void run(){
	t = System.currentTimeMillis();
		float C=10000;
		float[] alpha=new float[N];//alpha
		float bias=0;
		int i,j,k;
		float tol=0.01f;
		int changed_alphas;//number of changed alphas
		int it=0;//iteration index
		int maxit=100;
		float Ei=0;
		float Ej=0;
		int y_changed_i;
		int y_changed_j;
		int y_changed_k;
		float ai=0;
		float aj=0;
		float L = 0;
		float H = 0;
		float eta;
		float bi,bj;
		Random rand=new Random();
		float[] w3 = new float[dim];
		float[] w2 = new float[dim+1];
		// Initialize alpha's with 0
		for(i = 0; i < N; i++)
			alpha[i] = 0;
		
		while(it < maxit){//number of iterations less than maximum
			it++;
			changed_alphas = 0;
			
			for(i = 0; i < N; i++){
                Ei = 0;
				y_changed_i = svm.ind.V[i].cl.Y == 0? -1 : 1;
				for(k = 0; k < N; k++){
					y_changed_k = svm.ind.V[k].cl.Y == 0? -1 : 1;
					Ei += alpha[k] * y_changed_k * dotProduct(svm.ind.V[k].X, svm.ind.V[i].X);
				}
				
				Ei = Ei + bias - y_changed_i;
				
				if((alpha[i] < C&&(y_changed_i * Ei)< -tol) ||(alpha[i] > 0&&(y_changed_i * Ei > tol))){
					// Selecting a random j != i
					do{
						j = rand.nextInt(N);
					}while(j == i);
					
                    y_changed_j = svm.ind.V[j].cl.Y == 0? -1 : 1;
					
                    Ej = 0;
					
					for(k = 0; k < N; k++){
						y_changed_k = svm.ind.V[k].cl.Y == 0? -1 : 1;
						Ej += alpha[k] * y_changed_k * dotProduct(svm.ind.V[k].X, svm.ind.V[j].X);
					}
					Ej = Ej + bias - y_changed_j;
					ai=alpha[i];//alpha_i old
					aj=alpha[j];//alpha_j old
					if(y_changed_i != y_changed_j){
						L = Math.max(0, alpha[j] - alpha[i]);
						H = Math.min(C, C + alpha[j] - alpha[i]);
					}
					else if(y_changed_i == y_changed_j){
						L = Math.max(0, alpha[i] + alpha[j] - C);
						H = Math.min(C, alpha[i] + alpha[j]);
					}
					if(L == H)
						continue;
					eta = 2 * dotProduct(svm.ind.V[i].X, svm.ind.V[j].X) 
						- dotProduct(svm.ind.V[i].X, svm.ind.V[i].X) 
						- dotProduct(svm.ind.V[j].X, svm.ind.V[j].X);
					if (eta >= 0)
						continue;
					alpha[j] = alpha[j] + y_changed_j * (Ej - Ei) / eta;
					if(alpha[j] > H)
						alpha[j] = H;
					else if(alpha[j] < L)
						alpha[j] = L;
					if(Math.abs(alpha[j] - aj) < tol)
						continue;
					alpha[i] = alpha[i] - y_changed_i * y_changed_j * (alpha[j] - aj);
					bi = bias - Ei 
					- y_changed_i * (alpha[i] - ai) * dotProduct(svm.ind.V[i].X,svm.ind.V[i].X) 
					- y_changed_j * (alpha[j] - aj) * dotProduct(svm.ind.V[j].X,svm.ind.V[i].X);
					bj = bias - Ej 
					- y_changed_i * (alpha[i] - ai) * dotProduct(svm.ind.V[i].X,svm.ind.V[j].X) 
					- y_changed_j * (alpha[j] - aj) * dotProduct(svm.ind.V[j].X,svm.ind.V[j].X);
					if(0 < alpha[i] && alpha[i] < C)
						bias = bi;
					else if(0 < alpha[j] && alpha[j] < C)
						bias = bj;
					else
						bias = (bi + bj)/2;
					changed_alphas = changed_alphas + 1;
				}
			}
			if(changed_alphas == 0)
                it++;
            else it = 0;
		}
		
		// Representing the support vectors
		
		ArrayList<Integer> vS = new ArrayList<Integer>();
		for(i = 0; i < N; i++){
			if(alpha[i] == Float.NaN)
				continue;
			if(alpha[i] != 0)
				vS.add(i);
		}
		
		for(i = 0; i < vS.size(); i++){
			int index = vS.get(i);
			int y_i = svm.ind.V[index].cl.Y == 0 ? -1 : 1;
			for(k = 0; k < dim; k++)
				w3[k] += alpha[index] * svm.ind.V[index].X[k] * y_i;
		}
	
		for(i = 0; i < dim; i++)
			w2[i] = w3[i];
		
		 w2[dim] = bias;
		
		svm.outd.computing_time = System.currentTimeMillis() - t;
		
		float[] b = translate(w2,svm.ind.V);
		w2[dim] = -b[2];
		float[] w0 = new float[dim+1];
		for(j=0; j<dim; j++) w0[j] = w2[j];
		w0[dim] = -b[0];
		float[] w1 = new float[dim+1];
		for(j=0; j<dim; j++) w1[j] = w2[j];
		w1[dim] = -b[1];			
		if(dim==2) svm.design.setPointsOfMaxLine(w2,w0,w1);	
		
		svm.outd.w = w2;
		svm.outd.accuracy = getAccuracy(w2);
		svm.outd.margin = b[3];
		svm.outd.stages_count = 1;
		svm.outd.max_stages_count = 1;
		svm.outd.showInputData();
		svm.outd.showOutputData();
		svm.design.calculates = false;
		svm.design.repaint();
		svm.control.start.enable(false);		
	
	}
	
	 public static float[] translate(float[] w, io.Vector[] V){
		int dim = w.length - 1;
		int N = V.length;
		float max = Float.MIN_VALUE, min;
		int imax = -1, imin;
		for(int i = 0; i < N; i++){
			float d = distFromHiperplanToVector(w, V[i]);
			if(d > max){max = d; imax = i;}
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
		b[3] = distFromHiperplanToVector(w, V[imin]); 
		return b;
	}

	public static float distFromHiperplanToVector(float[] w, io.Vector V){
		float dist = 0, norm = 0;
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