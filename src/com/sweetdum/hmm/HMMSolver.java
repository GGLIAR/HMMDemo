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
 * HMMSolver provider the method to solve the three center questions of HMM
 * It need a hmm representation and an observations sequence as input
 *
 * Created by Mengxiao Lin on 2015/10/2.
 */
public class HMMSolver {
    private final double EPS=1e-15;
    HMMRepresentation hmm;
    HMMForwardBackwardManipulator fbManipulator;
    int o[];
    public HMMSolver(HMMRepresentation hmm,int observations[]){
        this.hmm=hmm;
        o=observations;
        this.fbManipulator=new HMMForwardBackwardManipulator(hmm,observations);
    }

    /**
     * input the observation sequence into the solver
     * @param observations
     */
    public void readNewObservation(int observations[]){
        o=observations;
        this.fbManipulator=new HMMForwardBackwardManipulator(hmm,observations);
    }

    /**
     * get the probability of the observation sequence with the HMM
     * @return the probability
     */
    public double getObservationProbability(){
        double ret=0;
        int time=o.length;
        int N=hmm.stateCount;
        for (int t=0;t<=time;++t){
            for (int i=0;i<N;++i){
                ret+=fbManipulator.alpha[t][i]*fbManipulator.beta[t][i];
            }
        }
        ret/=time;
        return ret;
    }

    /**
     * Get the state sequence which explain the observations in the largest
     * probability way
     * Use Viterbi algorithm
     *
     * @return the state ids sequence
     */
    public int[] getBestStateSequence(){
        int endTime=o.length;
        int N=hmm.stateCount;
        double delta[][]=new double[endTime+1][N];
        int phi[][]=new int[endTime+1][N];

        for (int i=0;i<N;++i) delta[0][i]=hmm.pi[i];
        for (int t=1;t<=endTime;++t) {
            for (int i = 0; i < N; ++i) {
                for (int j = 0; j < N; ++j) {
                    //transfer from j to i
                    double tmp = delta[t - 1][j] * hmm.a[j][i] * hmm.b[j][i][o[t - 1]];
                    if (tmp > delta[t][i]) {
                        delta[t][i] = tmp;
                        phi[t][i] = j;
                    }
                }
            }
        }

        int ret[]=new int[endTime+1];
        for (int i=0;i<N;++i) {
            if (delta[endTime][i]>delta[endTime][ret[endTime]]) ret[endTime]= i;
        }
        for (int t=endTime-1;t>=0;--t){
            ret[t]=phi[t+1][ret[t+1]];
        }
        return ret;
    }

    /**
     * re-estimate the hmm parameter with forward-backward algorithm
     */
    public void estimateParameter(){
        int endTime=o.length;
        int N=hmm.stateCount;
        double p[][][]=new double[endTime][N][N];
        for (int t=0;t<endTime;++t){
            double count=EPS;
            for (int i=0;i<N;++i){
                for (int j=0;j<N;++j){
                    p[t][i][j]=fbManipulator.alpha[t][i]*hmm.a[i][j]*hmm.b[i][j][o[t]]*fbManipulator.beta[t+1][j];
                    count +=p[t][j][j];
                }
            }
            for (int i=0;i<N;++i){
                for (int j=0;j<N;++j){
                    p[t][i][j]/=count;
                }
            }
        }
        double pSumThroughTime[][]=new double[N][N];
        for (int i=0;i<N;++i){
            for (int j=0;j<N;++j){
                pSumThroughTime[i][j]=EPS;
            }
        }
        for (int t=0;t<endTime;++t){
            for (int i=0;i<N;++i){
                for (int j=0;j<N;++j){
                    pSumThroughTime[i][j]+=p[t][i][j];
                }
            }
        }

        //rebuild gamma
        for (int t=0;t<endTime;++t){
            for (int i=0;i<N;++i) fbManipulator.gamma[t][i]=0;
            for (int i=0;i<N;++i){
                for (int j=0;j<N;++j){
                    fbManipulator.gamma[t][i]+=p[t][i][j];
                }
            }
        }

        double gammaCount=EPS;
        //maximum parameter
        for (int i=0;i<N;++i){
            gammaCount+=fbManipulator.gamma[0][i];
            hmm.pi[i]=fbManipulator.gamma[0][i];
        }
        for (int i=0;i<N;++i){
            hmm.pi[i]/=gammaCount;
        }

        for (int i=0;i<N;++i) {
            gammaCount = EPS;
            for (int t = 0; t < endTime; ++t) gammaCount += fbManipulator.gamma[t][i];
            for (int j = 0; j < N; ++j) {
                hmm.a[i][j] = pSumThroughTime[i][j] / gammaCount;
            }
        }
        for (int i=0;i<N;i++){
            for (int j=0;j<N;++j){
                double tmp[]=new double[hmm.observationCount];
                for (int t=0;t<endTime;++t){
                    tmp[o[t]]+=p[t][i][j];
                }
                for (int k=0;k<hmm.observationCount;++k){
                    hmm.b[i][j][k]=(tmp[k]+1e-8)/pSumThroughTime[i][j];
                }
            }
        }
        //rebuild fbManipulator
        fbManipulator=new HMMForwardBackwardManipulator(hmm,o);
    }
}
