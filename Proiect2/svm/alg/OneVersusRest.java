package alg;

import svm.SVM;
import io.*;
import java.util.*;
import io.Vector;

//*******Sequential Minimal Optimization (SMO) Algorithm*******
//***Reference***
//http://fourier.eng.hmc.edu/e176/lectures/ch9/node9.html

public class OneVersusRest extends Algorithm{
	SVM svm;
	ArrayList<Vector> v;//multimea de probe cu f(x) potrivit
	Clasa CP,CM;//f(x)
	int nC;//numarul de clase
	float C=10000;//costul din SMO
	ArrayList<Vector> vS;//vectorii suport
	int N;//numarul de vector
	float bias;//bias din SMO
	static float sigma=100;//sigma din calculul nucleului Gaussian
	float[] alpha;//alpha din SMO
	
	static class KernelGaussian{//stocam datele pentru care am aplicat SMO si datele de iesire (alpha,bias)
		public  ArrayList<Vector> vectorF;//vectorii 
		public  float[] alphaF;//alpha
		public  ArrayList<Vector> vSF;//vectorii support
		public  float biasF;//bias
		public int NF;//numarul de vectori
		
		public KernelGaussian(ArrayList<Vector> vector,float[] alpha,ArrayList<Vector> vS,float bias,int N){
			this.vectorF=vector;
			this.alphaF=alpha;
			this.vSF=vS;//alpha>0&&alpha<C
			this.biasF=bias;
			this.NF=N;
		}
		
		
	}
	
	ArrayList<KernelGaussian> classifiers;//date corespunzatoare fiecarei clase dupa rezolvarea SMO cu multimea de probe corespunzatoare
	
	public OneVersusRest(SVM svm){
		super(svm);
		this.svm=svm;
		nC=svm.ind.classes.length;
		CP=new Clasa("",1,null);
		CM=new Clasa("",-1,null);
		if(svm.ind.V!=null){
			name="OneVersusRest";
			svm.outd.algorithm = name;
			svm.outd.max_stages_count =1;
			svm.outd.showInputData();
		}
	}
	
	//K(x,y)=exp(-||x-y||^2/2*sigma^2)
	float KGaussian(float[] x,float[] y){
		return (float)(Math.exp(-(dotProduct(x,x)+dotProduct(y,y)-2*(float)Math.sqrt(dotProduct(x,y))))/(2*sigma*sigma));
	}
	
	float dotProduct(float[] x,float[] y){
		float sum=0;
		for(int i=0;i<x.length;i++){
			sum+=x[i]*y[i];
		}
		return sum;
	}
	
	float[] calculAlpha(float[] alpha){
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
		while(it < maxit){//number of iterations less than maximum
			it++;
			changed_alphas = 0;
			
			for(i = 0; i < N; i++){
				Ei = 0;
				Vector v_i=v.get(i); 
				y_changed_i = v_i.cl.Y == 0? -1 : 1;
				for(k = 0; k < N; k++){	
					Vector v_k=v.get(k);
					y_changed_k = v_k.cl.Y == 0? -1 : 1;
					Ei += alpha[k] * y_changed_k *KGaussian(v_i.X, v_k.X);//problema calcul kgaussian
				}
				
				Ei = Ei + bias - y_changed_i;
				
				if((alpha[i] < C&&(y_changed_i * Ei)< -tol) ||(alpha[i] > 0&&(y_changed_i * Ei > tol))){
					// Selecting a random j != i
					do{
						j = rand.nextInt(N);
					}while(j == i);
					Vector v_j=v.get(j);
					y_changed_j = v_j.cl.Y == 0? -1 : 1;
					
					Ej = 0;
					
					for(k = 0; k < N; k++){
						Vector v_k=v.get(k);
						y_changed_k = v_k.cl.Y == 0? -1 : 1;
						Ej += alpha[k] * y_changed_k * KGaussian(v_k.X, v_j.X);
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
					eta = 2 *KGaussian(v_i.X, v_j.X) 
						- KGaussian(v_i.X, v_i.X) 
						- KGaussian(v_j.X, v_j.X);
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
					- y_changed_i * (alpha[i] - ai) * KGaussian(v_i.X,v_i.X) 
					- y_changed_j * (alpha[j] - aj) * KGaussian(v_j.X,v_i.X);
					bj = bias - Ej 
					- y_changed_i * (alpha[i] - ai) * KGaussian(v_i.X,v_j.X) 
					- y_changed_j * (alpha[j] - aj) * KGaussian(v_j.X,v_j.X);
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
		return alpha;
	}
	
	public float evaluate(float[] inp,KernelGaussian kg){
			float res=0;
			for(int i=0;i<kg.vSF.size();i++){
				Vector x_i=kg.vSF.get(i);
				float y_changed_i=x_i.cl.Y<0?-1:1;
				res+=kg.alphaF[i]*y_changed_i*KGaussian(inp,x_i.X);
			}
			return res-kg.biasF;//calcul hk in x dat
		}
	
	int classify(float[] in){
		int max_i=-1;
		float v_max=Float.MIN_VALUE;
		for(int i=0;i<nC;i++){
			KernelGaussian kg=classifiers.get(i);
			float val=evaluate(in,kg);
			if(val>v_max){
				v_max=val;
				max_i=i;
			}
		}
		return max_i;
	}
	
	void param(ArrayList<Vector> vect){
		v=new ArrayList<Vector>();
		v=vect;
		N=vect.size();
		vS=new ArrayList<Vector>();
		bias=0;
		alpha=new float[N];
		alpha=calculAlpha(alpha); 
		for(int i=0;i<N;i++){
			if(alpha[i]>0&&alpha[i]<C){
				vS.add(vect.get(i));
			}
		}
	}
	
	public void run(){
		v=new ArrayList<Vector>();
		t=System.currentTimeMillis();
		classifiers=new ArrayList<KernelGaussian>();
		for(int i=0;i<nC;i++){
			v.clear();
			for(int k=0;k<svm.ind.V.length;k++){
				v.add(new Vector(svm.ind.V[k].X,svm.ind.V[k].cl.Y==i?CP:CM));//adaugam vectorii folosind
			}
			param(v);//calculam datele SMO
			classifiers.add(new KernelGaussian(v,alpha,vS,bias,N));
		}
		float[] input=new float[svm.ind.V[0].X.length];
		for(int i=0;i<input.length;i++){
			input[i]=0;
		}
		int arg=classify(input);
		System.out.println(arg);
        svm.outd.w=new float[dim+1];
        svm.outd.stages_count = 1;
        svm.outd.max_stages_count = 1;
        svm.outd.computing_time = System.currentTimeMillis() - t;
        svm.outd.showInputData();
        svm.outd.showOutputData();
        svm.design.calculates = false;
        svm.design.repaint();
        svm.control.start.enable(false);
	}
}