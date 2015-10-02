/*
Copyright (c) 2015 Mengxiao Lin

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

package com.sweetdum.hmm;

/**
 * HMMRepresentation stores the parameters of a hidden markov model.
 * With this representation object, user can use the model to create sequences.
 * Or modify the representation with data.
 *
 * Created by Mengxiao Lin on 2015/10/2.
 */
public class HMMRepresentation {
    /**
     * a: state transfer matrix
     * a[i][j] := P (s[i]->s[j] | s[i]);
     */
    double a[][];
    /**
     * b: observation symbol emit possibilities
     * b[i][j][k] := P(o[k] | s[i]->s[j])
     */
    double b[][][];
    /**
     * pi: the initial probabilities of all states
     */
    double pi[];
    /**
     * the count of state symbols
     */
    int stateCount;
    /**
     * the count of observation symbols;
     */
    int observationCount;

    public HMMRepresentation(int stateCount,int observationCount){
        this.stateCount=stateCount;
        this.observationCount=observationCount;
        a=new double[stateCount][stateCount];
        b=new double[stateCount][stateCount][observationCount];
        pi=new double[stateCount];
    }

    public int getStateCount() {
        return stateCount;
    }

    public int getObservationCount() {
        return observationCount;
    }

    public double getA(int i,int j) {
        return a[i][j];
    }

    public double getB(int i,int j,int k) {
        return b[i][j][k];
    }

    public double getPi(int i) {
        return pi[i];
    }

    public void setA(int i,int j, double a) {
        this.a[i][j] = a;
    }

    public void setB(int i,int j,int k,double b) {
        this.b[i][j][k] = b;
    }

    public void setPi(int i,double pi) {
        this.pi[i] = pi;
    }

    /**
     * produce an observations sequence randomly with the hmm
     * @param length the length of the observations
     * @return an array of observations sequence
     */
    public int[] produceObservations(int length){
        int nowState=0,nextState=0;
        int ret[]=new int[length];
        double random=Math.random();

        //produce the first state
        for (int i=0;i<stateCount;++i){
            random-=pi[i];
            if (random<=0){
                nowState=i;
                break;
            }
        }

        for (int pos=0;pos<length;++pos){
            //produce the next state
            random=Math.random();
            for (int i=0;i<stateCount;++i){
                random-=a[nowState][i];
                if (random<=0){
                    nextState=i;
                    break;
                }
            }
            //produce an observation
            random=Math.random();
            for (int i=0;i<observationCount;++i){
                random-=b[nowState][nextState][i];
                if (random<=0){
                    ret[pos]=i;
                    break;
                }
            }
            nowState=nextState;
        }
        return ret;
    }

    public int[] produceObservationFromState(int[] state){
        int ret[]=new int[state.length-1];
        double random;
        for (int i=1;i<state.length;++i){
            random=Math.random();
            for (int k=0;k<observationCount;++k){
                random-=b[state[i-1]][state[i]][k];
                if (random<=0){
                    ret[i-1]=k;
                    break;
                }
            }
        }
        return ret;
    }
    /**
     * init the parameter to uniform distribution
     */
    public void uniformParameter(){
        double count=0;
        for (int i=0;i<stateCount;++i){
            pi[i]=Math.random();count+=pi[i];
        }
        for (int i=0;i<stateCount;++i){
            pi[i]/=count;
        }
        for (int i=0;i<stateCount;++i){
            count=0;
            for (int j=0;j<stateCount;++j){
                a[i][j]=Math.random();
                count+=a[i][j];
            }
            for (int j=0;j<stateCount;++j) a[i][j]/=count;
        }
        for (int i=0;i<stateCount;++i){
            for (int j=0;j<stateCount;++j){
                count=0;
                for (int k=0;k<observationCount;++k){
                    b[i][j][k]=Math.random();
                    count+=b[i][j][k];
                }
                for (int k=0;k<observationCount;++k){
                    b[i][j][k]/=count;
                }
            }
        }
    }

    /**
     * output the parameter of the hmm representation
     */
    public void printParameter(){
        for (int i=0;i<2;++i){
            for (int j=0;j<2;++j){
                System.out.print("A[" + i + "][" + j + "]=" + a[i][j] + " ");
                for (int k=0;k<4;++k){
                    System.out.print("B["+i+"]["+j+"]["+k+"]="+b[i][j][k]+" ");
                }
                System.out.println();
            }
        }
        for (int i=0;i<2;++i){
            System.out.println("Pi["+i+"]="+pi[i]);
        }
    }

    /**
     * get the difference of parameters
     * @param hmm the HMM representation to compare
     * @return the square difference
     */
    public double diffWithOtherHMM(HMMRepresentation hmm){
        double error=0;
        for (int i=0;i<stateCount;++i){
            for (int j=0;j<stateCount;++j){
                error+=(hmm.a[i][j]-a[i][j])*(hmm.a[i][j]-a[i][j]);
            }
        }
        for (int i=0;i<stateCount;++i){
            for (int j=0;j<stateCount;++j){
                for (int k=0;k<observationCount;++k) {
                    error += (hmm.b[i][j][k] - b[i][j][k]) * (hmm.b[i][j][k] - b[i][j][k]);
                }
            }
        }
        for (int i=0;i<stateCount;++i){
            error+=(hmm.pi[i]-pi[i])*(hmm.pi[i]-pi[i]);
        }
        error/=stateCount+stateCount*stateCount+stateCount*stateCount*observationCount;
        return error;
    }

    @Override
    public Object clone() {
        HMMRepresentation ret=new HMMRepresentation(stateCount,observationCount);
        for (int i=0;i<stateCount;++i){
            ret.pi[i]=pi[i];
            for (int j=0;j<stateCount;++j){
                ret.a[i][j]=a[i][j];
                for (int k=0;k<observationCount;++k){
                    ret.b[i][j][k]=b[i][j][k];
                }
            }
        }
        return ret;
    }
}
